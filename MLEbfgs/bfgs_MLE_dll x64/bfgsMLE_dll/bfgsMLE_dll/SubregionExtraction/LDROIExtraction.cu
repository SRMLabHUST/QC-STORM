/*
This program is free software: you can redistribute it and/or modify
it under the terms of the GNU LESSER GENERAL PUBLIC LICENSE as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU LESSER GENERAL PUBLIC LICENSE for more details.

You should have received a copy of the GNU LESSER GENERAL PUBLIC LICENSE
along with this program.  If not, see <https://www.gnu.org/licenses/>.
*/

#include "LDROIExtraction.h"


#include "bfgsMLE_core.h"

#include <math.h>


#define Max(a,b)    (((a) > (b)) ? (a) : (b))
#define Min(a,b)    (((a) < (b)) ? (a) : (b))




void LDROIExtractData_TypeDef::ExtractMolecules(unsigned short *pImgData, int ImageSource, LocalizationPara & LocPara, int StartFrame_CurBatch, int BatchedImageNum, cudaStream_t cstream)
{
	const int BatchedImageHigh = BatchedImageNum * LocPara.ImageHigh;

	const int BatchImagePixelNum = BatchedImageNum * LocPara.ImageWidth * LocPara.ImageHigh;


	// first batch, localization may be perfomed after multi extraction batches
	if ((FirstFrame < 0) && (EndFrame < 0))
	{
		FirstFrame = StartFrame_CurBatch;
		EndFrame = FirstFrame + BatchedImageNum - 1;
	}
	else
	{
		EndFrame = EndFrame + BatchedImageNum;
	}

	// get raw image
	if (ImageSource == ImageSource_CPU_Pinned)
	{
		// image from CPU
		cudaMemcpyAsync(d_RawImg, pImgData, BatchImagePixelNum * sizeof(unsigned short), cudaMemcpyHostToDevice, cstream); // h_RawImg
	}
	else if (ImageSource == ImageSource_CPU_Normal)
	{
		// image from CPU
		cudaMemcpy(d_RawImg, pImgData, BatchImagePixelNum * sizeof(unsigned short), cudaMemcpyHostToDevice); // h_RawImg
	}
	else if (ImageSource == ImageSource_GPU)
	{
		// always cause problem
		//		cudaMemcpy(d_RawImg, pImgData, BatchImagePixelNum * sizeof(unsigned short), cudaMemcpyDeviceToDevice); // h_RawImg
		//		cudaMemcpyAsync(d_RawImg, pImgData, BatchImagePixelNum * sizeof(unsigned short), cudaMemcpyDeviceToDevice, cstream);

		printf("data error\n");
		return;
	}
	else
	{
		printf("data error\n");
		// image data error
		return;
	}

	// ImageFiltering
	ImageFiltering(LocPara.LocType, LocPara.MultiEmitterFitEn, LocPara.ImageWidth, LocPara.ImageHigh, BatchedImageNum, cstream);

	ROIExtraction(LocPara.ROISize, LocPara.LocType, LocPara.MultiEmitterFitEn, LocPara.ImageWidth, LocPara.ImageHigh, BatchedImageNum, StartFrame_CurBatch, cstream);

}


void LDROIExtractData_TypeDef::ROIExtraction(int ROISize, int LocType, int MultiEmitterFitEn, int ImageWidth, int ImageHigh, int BatchedImageNum, int StartFrame, cudaStream_t cstream)
{
	const int BatchedImageHigh = BatchedImageNum*ImageHigh;
	const int ROIWholeSize = ROISize*(ROISize + 1);


	cudaMemsetAsync(d_ValidROINum, 0, sizeof(int), cstream);

	cudaMemsetAsync(d_MoleculePosImage, 0, ImageWidth * BatchedImageHigh * sizeof(short), cstream);

	// find possible ROI position

	ROIFindingLD(d_RawImg_Smoothed, d_MoleculePosImage, d_ValidROINum, d_ROIMarkInfArray, ImageVariance, LocType, MultiEmitterFitEn, ROISize, ImageWidth, BatchedImageHigh, ImageHigh, cstream);


	cudaMemcpyAsync(h_ValidROINum, d_ValidROINum, sizeof(int), cudaMemcpyDeviceToHost, cstream);
	cudaStreamSynchronize(cstream);

	int ROINumber = *h_ValidROINum;


	// pre-clasify molecule by single or multi by ROI molecule distance

	DetectedROIClasify(d_MoleculePosImage, d_ROIMarkInfArray, ROINumber, MultiEmitterFitEn, ROISize, ImageWidth, cstream);

	cudaStreamSynchronize(cstream);
	
	int ROINum_CurBatch = *h_ValidROINum;


	// ROI extraction

	SubregionExtractionLD(d_RawImg, d_ImageROI, d_ROIMarkInfArray, ROISize, ROINum_CurBatch, ImageWidth, BatchedImageHigh, ImageHigh, StartFrame, cstream);


	int ROIAddrOffset = TotalROINumber*ROIWholeSize;

	cudaMemcpyAsync(&h_ImageROI[ROIAddrOffset], d_ImageROI, ROINum_CurBatch * ROIWholeSize * sizeof(short), cudaMemcpyDeviceToHost, cstream);
	cudaStreamSynchronize(cstream); // wait task of this stream finish


#if(WLE_ENABLE == 1)

	// estimate WLE parameter
	int WLEParaAddrOffset = TotalROINumber*WLE_ParaNumber;

	WLEParameterEstimator->WLEParameterEstimate(d_ImageROI, LocType, MultiEmitterFitEn, ROISize, ROINum_CurBatch, cstream);// &h_ImageROI[ROIAddrOffset]

	cudaMemcpyAsync(&WLEParameterEstimator->h_WLEPara[WLEParaAddrOffset], WLEParameterEstimator->d_WLEPara, ROINum_CurBatch * WLE_ParaNumber * sizeof(float), cudaMemcpyDeviceToHost, cstream);
	cudaStreamSynchronize(cstream); // wait task of this stream finish

#endif // WLE_ENABLE

	TotalROINumber += ROINum_CurBatch;

}

void LDROIExtractData_TypeDef::ImageVarianceCalc(unsigned short *d_iRawImg, int ImageWidth, int ImageHigh, cudaStream_t cstream)
{
	int ImageCenterX = ImageWidth / 2;
	int ImageCenterY = ImageHigh / 2;

	int StdImageWidth = ImageWidth / 2 - 5;
	int StdImageHigh = ImageHigh / 2 - 5;


	StdImageWidth = Min(StdImageWidth, 256);
	StdImageHigh = Min(StdImageHigh, 256);


	StdImageWidth = Max(StdImageWidth, 1);
	StdImageHigh = Max(StdImageHigh, 1);


	int StdImageSize = Min(StdImageWidth, StdImageHigh);

	cudaMemsetAsync(d_MeanDataX, 0, sizeof(float), cstream);
	cudaMemsetAsync(d_MeanDataX2, 0, sizeof(float), cstream);

	ImageVarianceCalcOnSelectRegion(d_iRawImg, d_MeanDataX, d_MeanDataX2, ImageCenterX, ImageCenterY, StdImageSize, ImageWidth, ImageHigh, cstream);


	cudaMemcpyAsync(h_MeanDataX, d_MeanDataX, sizeof(float), cudaMemcpyDeviceToHost, cstream);
	cudaMemcpyAsync(h_MeanDataX2, d_MeanDataX2, sizeof(float), cudaMemcpyDeviceToHost, cstream);

	cudaStreamSynchronize(cstream);


	ImageVariance = *h_MeanDataX2 - (*h_MeanDataX)*(*h_MeanDataX);

	if (ImageVariance < 100.0f) ImageVariance = 100.0f;


	//	printf("ImageVariance:%f\n", ImageVariance);
}



int LDROIExtractData_TypeDef::GetAccumulatedROINum()
{
	return TotalROINumber;
}


void LDROIExtractData_TypeDef::ResetROINum()
{
	TotalROINumber = 0;

	FirstFrame = -1;
	EndFrame = -1;

}

float * LDROIExtractData_TypeDef::Get_h_WLEPara()
{
	return WLEParameterEstimator->h_WLEPara;
}

float * LDROIExtractData_TypeDef::Get_d_WLEPara()
{
	return WLEParameterEstimator->d_WLEPara;

}

int LDROIExtractData_TypeDef::GetMaxBatchedNumForCurrentImageSize(int ImageWidth, int ImageHigh)
{
	int MaxBatchImgNum = (2048 * 2 / (ImageWidth + ImageHigh));

	MaxBatchImgNum = MaxBatchImgNum / 2 * 2;

	if (MaxBatchImgNum < 1)MaxBatchImgNum = 1;

	return MaxBatchImgNum;
}


void LDROIExtractData_TypeDef::Init(LocalizationPara & LocPara)
{
	cudaError_t err;


	int MaxBatchImgNum = GetMaxBatchedNumForCurrentImageSize(LocPara.ImageWidth, LocPara.ImageHigh);

	WLEParameterEstimator = new WLEParameterEstimation_TypeDef();
	WLEParameterEstimator->Init(LocPara);


	// raw image
	err = cudaMallocHost((void **)&h_RawImg, MaxBatchedImageSize * sizeof(short));
	HandleErr(err, "cudaMalloc h_RawImg");

	err = cudaMalloc((void **)&d_RawImg, MaxBatchedImageSize * sizeof(short));
	HandleErr(err, "cudaMalloc d_RawImg");

	err = cudaMalloc((void **)&d_RawImg_Smoothed, MaxBatchedImageSize * sizeof(short));
	HandleErr(err, "cudaMalloc d_RawImg_Smoothed");


	err = cudaMalloc((void **)&d_BackgroundImage, MaxBatchedImageSize * sizeof(short));
	HandleErr(err, "cudaMalloc d_BackgroundImage");

	err = cudaMalloc((void **)&d_LineFilterImage_t, MaxBatchedImageSize * sizeof(short));
	err = cudaMalloc((void **)&d_LineFilterImage_t1, MaxBatchedImageSize * sizeof(short));


	err = cudaMalloc((void **)&d_MoleculePosImage, MaxBatchedImageSize * sizeof(short));

	// image std use

	err = cudaMallocHost((void **)&h_MeanDataX, sizeof(float));
	err = cudaMallocHost((void **)&h_MeanDataX2, sizeof(float));

	err = cudaMalloc((void **)&d_MeanDataX, sizeof(float));
	err = cudaMalloc((void **)&d_MeanDataX2, sizeof(float));


	// image filter

	//	err = cudaMallocHost((void **)&h_LineFilterH_Bkg, LineFilterSize * sizeof(float));
	//	err = cudaMalloc((void **)&d_LineFilterH_Bkg, LineFilterSize * sizeof(float));

	err = cudaMallocHost((void **)&h_LineFilterH_Signal, LineFilterSize * sizeof(float));
//	err = cudaMalloc((void **)&d_LineFilterH_Signal, LineFilterSize * sizeof(float));


	// extracted molecular ROI data
	const int ROIWholeSize = LocPara.ROISize*(LocPara.ROISize + 1);

	err = cudaMallocHost((void **)&h_ImageROI, MaxPointNum * ROIWholeSize * sizeof(unsigned short));
	err = cudaMalloc((void **)&d_ImageROI, MaxPointNum * ROIWholeSize * sizeof(unsigned short));

	err = cudaMallocHost((void **)&h_ValidROINum, sizeof(int));
	err = cudaMalloc((void **)&d_ValidROINum, sizeof(int));


	err = cudaMalloc((void **)&d_ROIMarkInfArray, MaxPointNum * ROIMarkInfNum * sizeof(int));


	//
	FilterInit(LocPara.ROISize, LocPara.LocType, LocPara.MultiEmitterFitEn);

	ResetROINum();

}


void LDROIExtractData_TypeDef::Deinit()
{
	cudaError_t err;

	WLEParameterEstimator->Deinit();

	delete WLEParameterEstimator;


	// raw image filtering
	err = cudaFreeHost(h_RawImg);
	HandleErr(err, "cudaFreeHost h_RawImg");

	err = cudaFree(d_RawImg);
	HandleErr(err, "cudaFree d_RawImg");

	err = cudaFree(d_RawImg_Smoothed);


	err = cudaFree(d_BackgroundImage);
	HandleErr(err, "cudaFree d_BackgroundImage");

	err = cudaFree(d_LineFilterImage_t);
	err = cudaFree(d_LineFilterImage_t1);

	// image std use

	err = cudaFreeHost(h_MeanDataX);
	err = cudaFreeHost(h_MeanDataX2);

	err = cudaFree(d_MeanDataX);
	err = cudaFree(d_MeanDataX2);


	// image filter
	//	err = cudaFreeHost(h_LineFilterH_Bkg);
	//	err = cudaFree(d_LineFilterH_Bkg);

	err = cudaFreeHost(h_LineFilterH_Signal);
//	err = cudaFree(d_LineFilterH_Signal);


	// extracted molecular ROI data
	err = cudaFreeHost(h_ImageROI);
	err = cudaFree(d_ImageROI);

	err = cudaFreeHost(h_ValidROINum);
	err = cudaFree(d_ValidROINum);

	err = cudaFree(d_ROIMarkInfArray);

}



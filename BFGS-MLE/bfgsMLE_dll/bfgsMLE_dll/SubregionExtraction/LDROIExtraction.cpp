#include "stdafx.h"

#include "LDROIExtraction.h"


#include "bfgsMLE_core.h"


#define Max(a,b)    (((a) > (b)) ? (a) : (b))
#define Min(a,b)    (((a) < (b)) ? (a) : (b))




void LDROIExtractData_TypeDef::ExtractMolecules(unsigned short *pImgData, int ImageSource, LocalizationPara & LocPara, int StartFrame_CurBatch, int BatchedImageNum, cudaStream_t cstream)
{
	const int BatchedImageHigh = BatchedImageNum* LocPara.ImageHigh;

	const int BatchImagePixelNum = BatchedImageNum * LocPara.ImageWidth * LocPara.ImageHigh;


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

	/*
	cudaError_t cudaStatus;

	cudaStatus = cudaGetLastError();
	HandleErr(cudaStatus, "memory copy");
	*/

	// ImageFiltering
	ImageFiltering(LocPara.ImageWidth, LocPara.ImageHigh, BatchedImageNum, cstream);

	ROIExtraction(LocPara.ROISize, LocPara.ImageWidth, LocPara.ImageHigh, BatchedImageNum, StartFrame_CurBatch, cstream);

}


void LDROIExtractData_TypeDef::ImageFiltering(int ImageWidth, int ImageHigh, int BatchedImageNum, cudaStream_t cstream)
{
	const int BatchedImageHigh = BatchedImageNum*ImageHigh;


	// calculate the background intensity, then smooth the background
	BackgroundIntensityCalc(d_RawImg, d_BackgroundImage, ImageWidth, BatchedImageHigh, cstream);

	ImageConv_Seperable_H(d_BackgroundImage, d_LineFilterImage_t, d_LineFilterH_Bkg, ImageWidth, BatchedImageHigh, cstream);
	ImageConv_Seperable_V(d_LineFilterImage_t, d_BackgroundImage, d_LineFilterH_Bkg, ImageWidth, BatchedImageHigh, cstream);


	// ImageSubstract d_RawImg by d_BackgroundImage
	ImageSubstract(d_RawImg, d_BackgroundImage, d_RawImg_Smoothed, ImageWidth, BatchedImageHigh, cstream);


	ImageVarianceCalc(d_RawImg_Smoothed, ImageWidth, BatchedImageHigh, cstream);


	// smooth d_RawImg
	ImageConv_Seperable_H(d_RawImg_Smoothed, d_LineFilterImage_t, d_LineFilterH_Signal, ImageWidth, BatchedImageHigh, cstream);
	ImageConv_Seperable_V(d_LineFilterImage_t, d_RawImg_Smoothed, d_LineFilterH_Signal, ImageWidth, BatchedImageHigh, cstream);

	// debug
	//	cudaMemcpyAsync(h_BackgroundImage, d_RawImg, BatchImagePixelNum * sizeof(short), cudaMemcpyDeviceToHost, cstream);

	cudaStreamSynchronize(cstream);
//	printf("var:%f\n", ImageVariance);

}



void LDROIExtractData_TypeDef::ROIExtraction(int ROISize, int ImageWidth, int ImageHigh, int BatchedImageNum, int StartFrame, cudaStream_t cstream)
{
	const int BatchedImageHigh = BatchedImageNum*ImageHigh;
	const int ROIWholeSize = ROISize*(ROISize + 1);

	const int MaxFluoNumPerImage = GetMaxFluoNumPerImage(BatchedImageNum);


	cudaMemsetAsync(d_ROINumPerImage, 0, BatchedImageNum * sizeof(int), cstream);


	ROIFindingLD(d_RawImg_Smoothed, d_ROIPosArray, d_ROINumPerImage, ImageVariance, ROISize, ImageWidth, BatchedImageHigh, ImageHigh, cstream);

	cudaMemcpyAsync(h_ROINumPerImage, d_ROINumPerImage, BatchedImageNum * sizeof(int), cudaMemcpyDeviceToHost, cstream);
	cudaStreamSynchronize(cstream);


	// extract molecules per image in batched image
	int RegionNum_CurImage = 0;
	int RegionNum_CurBatch = 0;


	for (int cnt = 0; cnt < BatchedImageNum; cnt++)
	{
		// since whole memory is seperated for batch frames, and there are zero gaps, see gpuSubregionFindingLD
		RegionNum_CurImage = h_ROINumPerImage[cnt];


		unsigned short* d_ROIMem_CurImage = &d_ROIMem[RegionNum_CurBatch*ROIWholeSize];

		int *d_ROIPosArray_CurImage = &d_ROIPosArray[cnt*MaxFluoNumPerImage*RegionPosInfNum];


		SubregionExtractionLD(d_RawImg, d_ROIMem_CurImage, d_ROIPosArray_CurImage, ROISize, RegionNum_CurImage, ImageWidth, BatchedImageHigh, ImageHigh, StartFrame, cstream);


		RegionNum_CurBatch += RegionNum_CurImage;
	}


	int ROIAddrOffset = TotalROINumber*ROIWholeSize;

	cudaMemcpyAsync(&h_ROIMem[ROIAddrOffset], d_ROIMem, RegionNum_CurBatch * ROIWholeSize * sizeof(short), cudaMemcpyDeviceToHost, cstream);
	cudaStreamSynchronize(cstream); // wait task of this stream finish

#if(WLE_ENABLE == 1)

	// estimate WLE parameter
	int WLEParaAddrOffset = TotalROINumber*WLE_ParaNumber;

	WLEParameterEstimator->WLEParameterEstimate(&h_ROIMem[ROIAddrOffset], ROISize, RegionNum_CurBatch, cstream);

	cudaMemcpyAsync(&WLEParameterEstimator->h_WLEPara[WLEParaAddrOffset], WLEParameterEstimator->d_WLEPara, RegionNum_CurBatch * WLE_ParaNumber * sizeof(float), cudaMemcpyDeviceToHost, cstream);
	cudaStreamSynchronize(cstream); // wait task of this stream finish

#endif // WLE_ENABLE


	TotalROINumber += RegionNum_CurBatch;
}

void LDROIExtractData_TypeDef::FilterInit()
{
	// sigma 2
	h_LineFilterH_Bkg[0] = 0.15247f;
	h_LineFilterH_Bkg[1] = 0.22184f;
	h_LineFilterH_Bkg[2] = 0.25138f;
	h_LineFilterH_Bkg[3] = 0.22184f;
	h_LineFilterH_Bkg[4] = 0.15247f;


	// sigma 0.7
	h_LineFilterH_Signal[0] = 0.00962f;
	h_LineFilterH_Signal[1] = 0.20542f;
	h_LineFilterH_Signal[2] = 0.56991f;
	h_LineFilterH_Signal[3] = 0.20542f;
	h_LineFilterH_Signal[4] = 0.00962f;


	cudaStream_t loc_stream1;
	cudaStreamCreate(&loc_stream1);

	cudaMemcpyAsync(d_LineFilterH_Bkg, h_LineFilterH_Bkg, ImageFilter_KernelSize*sizeof(float), cudaMemcpyHostToDevice, loc_stream1);
	cudaMemcpyAsync(d_LineFilterH_Signal, h_LineFilterH_Signal, ImageFilter_KernelSize * sizeof(float), cudaMemcpyHostToDevice, loc_stream1);

	cudaStreamDestroy(loc_stream1);
}

void LDROIExtractData_TypeDef::ImageVarianceCalc(unsigned short *d_iRawImg, int ImageWidth, int ImageHigh, cudaStream_t cstream)
{
	int ImageCenterX = ImageWidth / 2;
	int ImageCenterY = ImageHigh / 2;

	int StdImageWidth = ImageWidth / 4;
	int StdImageHigh = ImageHigh / 4;

	StdImageWidth = Min(StdImageWidth, 64);
	StdImageHigh = Min(StdImageHigh, 64);

	int StdImageSize = Min(StdImageWidth, StdImageHigh);

	cudaMemsetAsync(d_MeanDataX, 0, sizeof(float), cstream);
	cudaMemsetAsync(d_MeanDataX2, 0, sizeof(float), cstream);

	ImageVarianceCalcOnSelectRegion(d_iRawImg, d_MeanDataX, d_MeanDataX2, ImageCenterX, ImageCenterY, StdImageSize, ImageWidth, ImageHigh, cstream);


	cudaMemcpyAsync(h_MeanDataX, d_MeanDataX, sizeof(float), cudaMemcpyDeviceToHost, cstream);
	cudaMemcpyAsync(h_MeanDataX2, d_MeanDataX2, sizeof(float), cudaMemcpyDeviceToHost, cstream);

	cudaStreamSynchronize(cstream);


	ImageVariance = *h_MeanDataX2 - (*h_MeanDataX)*(*h_MeanDataX);

	if (ImageVariance < 100.0f)ImageVariance = 100.0f;

}



int LDROIExtractData_TypeDef::GetAccumulatedROINum()
{
	return TotalROINumber;
}


void LDROIExtractData_TypeDef::ResetROINum()
{
	TotalROINumber = 0;

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
	int MaxBatchImgNum = (2048 * 2048 / ImageWidth / ImageHigh) + 2;

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

	// background calculation
//	err = cudaMallocHost((void **)&h_BackgroundImage, MaxBatchedImageSize * sizeof(short));
//	HandleErr(err, "cudaMalloc h_BackgroundImage");


	err = cudaMalloc((void **)&d_BackgroundImage, MaxBatchedImageSize * sizeof(short));
	HandleErr(err, "cudaMalloc d_BackgroundImage");

	err = cudaMalloc((void **)&d_LineFilterImage_t, MaxBatchedImageSize * sizeof(short));

	// image std use

	err = cudaMallocHost((void **)&h_MeanDataX, sizeof(float));
	err = cudaMallocHost((void **)&h_MeanDataX2, sizeof(float));

	err = cudaMalloc((void **)&d_MeanDataX, sizeof(float));
	err = cudaMalloc((void **)&d_MeanDataX2, sizeof(float));


	// image filter

	err = cudaMallocHost((void **)&h_LineFilterH_Bkg, ImageFilter_KernelSize * sizeof(float));
	err = cudaMalloc((void **)&d_LineFilterH_Bkg, ImageFilter_KernelSize * sizeof(float));

	err = cudaMallocHost((void **)&h_LineFilterH_Signal, ImageFilter_KernelSize * sizeof(float));
	err = cudaMalloc((void **)&d_LineFilterH_Signal, ImageFilter_KernelSize * sizeof(float));


	// extracted molecular ROI data
	const int ROIWholeSize = LocPara.ROISize*(LocPara.ROISize + 1);

	err = cudaMallocHost((void **)&h_ROIMem, MaxPointNum * ROIWholeSize * sizeof(unsigned short));
	err = cudaMalloc((void **)&d_ROIMem, MaxPointNum * ROIWholeSize * sizeof(unsigned short));
	
	err = cudaMallocHost((void **)&h_ROINumPerImage, MaxBatchImgNum * sizeof(int));
	err = cudaMalloc((void **)&d_ROINumPerImage, MaxBatchImgNum * sizeof(int));

	err = cudaMalloc((void **)&d_ROIPosArray, MaxPointNum * RegionPosInfNum * sizeof(int));

	// consecutive fitting, to merge ROI by assign each consecutive ROI the average
	err = cudaMalloc((void **)&d_ForwardLinkID, MaxPointNum * sizeof(int));
	err = cudaMalloc((void **)&d_BackwardLinkID, MaxPointNum * sizeof(int));
	err = cudaMalloc((void **)&d_ConsecutiveNum, MaxPointNum * sizeof(int));


	//
	FilterInit();

	ResetROINum();

}


void LDROIExtractData_TypeDef::Deinit()
{
	cudaError_t err;

	WLEParameterEstimator->Deinit();

	delete WLEParameterEstimator;


	// raw image
	err = cudaFreeHost(h_RawImg);
	HandleErr(err, "cudaFreeHost h_RawImg");

	err = cudaFree(d_RawImg);
	HandleErr(err, "cudaFree d_RawImg");

	err = cudaFree(d_RawImg_Smoothed);

	// background calculation
//	err = cudaFreeHost(h_BackgroundImage);
//	HandleErr(err, "cudaFreeHost h_RawImg");

	err = cudaFree(d_BackgroundImage);
	HandleErr(err, "cudaFree d_BackgroundImage");

	err = cudaFree(d_LineFilterImage_t);

	// image std use

	err = cudaFreeHost(h_MeanDataX);
	err = cudaFreeHost(h_MeanDataX2);

	err = cudaFree(d_MeanDataX);
	err = cudaFree(d_MeanDataX2);


	// image filter
	err = cudaFreeHost(h_LineFilterH_Bkg);
	err = cudaFree(d_LineFilterH_Bkg);

	err = cudaFreeHost(h_LineFilterH_Signal);
	err = cudaFree(d_LineFilterH_Signal);


	// extracted molecular ROI data
	err = cudaFreeHost(h_ROIMem);
	err = cudaFree(d_ROIMem);

	err = cudaFreeHost(h_ROINumPerImage);
	err = cudaFree(d_ROINumPerImage);

	err = cudaFree(d_ROIPosArray);


	// consecutive fitting, to merge ROI by assign each consecutive ROI the average
	err = cudaFree(d_ForwardLinkID);
	err = cudaFree(d_BackwardLinkID);
	err = cudaFree(d_ConsecutiveNum);

}


/*
void LDROIExtractData_TypeDef::ROIMergeForConsecutiveFitting(int ROISize, int FluoNum, cudaStream_t cstream)
{
const int ROIWholeSize = ROISize*(ROISize + 1);

cudaMemcpyAsync(d_ROIMem, h_ROIMem, FluoNum * ROIWholeSize * sizeof(short), cudaMemcpyHostToDevice, cstream);
cudaMemcpyAsync(WLEParameterEstimator->d_WLEPara, WLEParameterEstimator->h_WLEPara, FluoNum * WLE_ParaNumber * sizeof(float), cudaMemcpyHostToDevice, cstream);


MergeConsecutiveROI(d_ROIMem, WLEParameterEstimator->d_WLEPara, ROISize, d_ForwardLinkID, d_BackwardLinkID, d_ConsecutiveNum, FluoNum, cstream);


cudaMemcpyAsync(h_ROIMem, d_ROIMem, FluoNum * ROIWholeSize * sizeof(short), cudaMemcpyDeviceToHost, cstream);
cudaMemcpyAsync(WLEParameterEstimator->h_WLEPara, WLEParameterEstimator->d_WLEPara, FluoNum * WLE_ParaNumber * sizeof(float), cudaMemcpyDeviceToHost, cstream);


cudaStreamSynchronize(cstream); // wait task of this stream finish

}
*/

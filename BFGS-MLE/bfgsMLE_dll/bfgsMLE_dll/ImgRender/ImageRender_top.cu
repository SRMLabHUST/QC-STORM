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

#include "ImageRender_2d.h"
#include "ImageRender_3d.h"
#include "ImageRender_top.h"
#include <malloc.h>


/*
RGBImageEncodeMode 0 : peak photon as weight for each molecule, rendered by localization precision calculated by CRLB
RGBImageEncodeMode 1 : 1 as weight for each molecule, rendered by localization precision calculated by CRLB
RGBImageEncodeMode 2 : 1 as weight for each molecule, rendered by fixed localization precision
*/


void ImageRenderData_TypeDef::FluoRenderTop(float *h_LocArry, LocalizationPara & LocPara, int RenderingMode, float FixedlocPrec, int FluoNum, cudaStream_t cstream)
{

	// for get display max intensity
	tRendFluoNum += FluoNum;
	if (tRendFluoNum > 20 * MaxPointNum)
	{
		tRendFluoNum = 0;
		ResetFillMaxVal(1);
	}

	cudaMemcpyAsync(d_LocArry, h_LocArry, FluoNum*OutParaNumGS2D*sizeof(float), cudaMemcpyHostToDevice, cstream);

	switch (LocPara.LocType)
	{
	case LocType_GS2D:
		// for 2d localization
		ImgRender_FluoRenderWithPrec_2D(h_LocArry, this, LocPara, RenderingMode, FixedlocPrec, FluoNum, cstream);

		break;

	case LocType_AS3D:
		// for 3d localization
	case LocType_DH3D:
		// for 3d localization
		ImgRender_FluoRender_3D(h_LocArry, this, LocPara, RenderingMode, FixedlocPrec, FluoNum, cstream);

		break;


	default:
		break;
	}


}


void ImageRenderData_TypeDef::ResetFillMaxVal(int Mode)
{
	cudaStream_t cstream;
	cudaStreamCreate(&cstream);

	if (Mode == 0)
	{
		cudaMemsetAsync(d_MaxImageVal, 0, sizeof(int), cstream);
		cudaMemsetAsync(d_HistMaxDat, 0, ImgRend_MaxDatHistLen*sizeof(int), cstream);
		cudaStreamSynchronize(cstream);

	}
	else
	{
		cudaMemcpyAsync(h_HistMaxDat, d_HistMaxDat, ImgRend_MaxDatHistLen * sizeof(int), cudaMemcpyDeviceToHost, cstream);
		cudaStreamSynchronize(cstream);

		int NData = h_HistMaxDat[0];
		if (NData <= 1)NData = 1;

		for (int cnt = 0; cnt < ImgRend_MaxDatHistLen; cnt++)
		{
			h_HistMaxDat[cnt] = h_HistMaxDat[cnt] * 200 / NData;
		}

		cudaMemcpyAsync(d_HistMaxDat, h_HistMaxDat, ImgRend_MaxDatHistLen * sizeof(int), cudaMemcpyHostToDevice, cstream);
		cudaStreamSynchronize(cstream);

	}

	cudaStreamDestroy(cstream);

}


void ImageRenderData_TypeDef::GetDispImgTop(LocalizationPara & LocPara, float BrightRatio, int oImgWidth, int oImgHigh, int cposX, int cposY, float DispZoom, cudaStream_t cstream)
{
	
	if (RendType_Is2D(LocPara.LocType))
	{
		//2d image 
		ImgRender_GetDispImg_2D(this, BrightRatio, oImgWidth, oImgHigh, cposX, cposY, DispZoom, LocPara.PixelZoom, LocPara.SRImageWidth, LocPara.SRImageHigh, cstream);

	}
	else
	{
		// astigmatism 3d image
		ImgRender_GetDispImg_3D(this, LocPara.MinZDepth, LocPara.MaxZDepth, LocPara.ColorMode_3D, BrightRatio, oImgWidth, oImgHigh, cposX, cposY, DispZoom, LocPara.PixelZoom, LocPara.SRImageWidth, LocPara.SRImageHigh, cstream);

	}
	cudaStreamSynchronize(cstream);

}

void ImageRenderData_TypeDef::GetSaveImgTop(LocalizationPara & LocPara, float BrightRatio, int RGBImageEncodeMode, cudaStream_t cstream)
{
	if (RendType_Is2D(LocPara.LocType))
	{
		//2d image 
		ImgRender_GetSaveImg_2D(this, BrightRatio, RGBImageEncodeMode, LocPara.SRImageWidth, LocPara.SRImageHigh, cstream);

	}
	else
	{
		// astigmatism 3d image
		ImgRender_GetSaveImg_3D(this, LocPara.MinZDepth, LocPara.MaxZDepth, LocPara.ColorMode_3D, BrightRatio, RGBImageEncodeMode, LocPara.SRImageWidth, LocPara.SRImageHigh, cstream);

	}
	cudaStreamSynchronize(cstream);

}

int ImageRenderData_TypeDef::GetDispMaxVal()
{
	cudaStream_t cstream;
	cudaStreamCreate(&cstream);

	cudaMemcpyAsync(h_MaxImageVal, d_MaxImageVal, sizeof(int), cudaMemcpyDeviceToHost, cstream);
	cudaStreamSynchronize(cstream);

	int MaxImageVal = *h_MaxImageVal;

	float LightRatio;
	LightRatio = GetLightRatio(h_HistMaxDat, d_HistMaxDat);
	
	MaxImageVal = MaxImageVal*LightRatio;

	if (MaxImageVal > 20000000)
	{
		ResetFillMaxVal(1);
	}
	cudaStreamDestroy(cstream);

	return MaxImageVal;
}

void ImageRenderData_TypeDef::ResetFillImgTop(LocalizationPara & LocPara)
{
	if (RendType_Is2D(LocPara.LocType))
	{
		//2d image 
		ImgRender_ResetFillImg_2D(this, LocPara.SRImageWidth, LocPara.SRImageHigh);

	}
	else
	{
		// color encoded depth 3d image
		ImgRender_ResetFillImg_3D(this, LocPara.SRImageWidth, LocPara.SRImageHigh);

	}
}


void ImageRenderData_TypeDef::GetMaxImgSizeFromLocArry(float *h_LocArry, float *d_LocArry, int *MaxImgWidth, int *MaxImgHigh, int FluoNum, cudaStream_t cstream)
{
	
	int *d_MaxImgWidth;
	int *d_MaxImgHigh;


	cudaMalloc((void **)&d_MaxImgWidth, sizeof(int));
	cudaMalloc((void **)&d_MaxImgHigh, sizeof(int));
	cudaMalloc((void **)&d_LocArry, PointNumTh * 2 * OutParaNumGS2D*sizeof(float));


	cudaMemsetAsync(d_MaxImgWidth, 0, sizeof(int), cstream);
	cudaMemsetAsync(d_MaxImgHigh, 0, sizeof(int), cstream);


	cudaMemcpyAsync(d_LocArry, h_LocArry, FluoNum*OutParaNumGS2D*sizeof(float), cudaMemcpyHostToDevice, cstream);

	int BlockDim = ThreadsPerBlock;
	int BlockNum = ((FluoNum + ThreadsPerBlock - 1) / ThreadsPerBlock);

	FindMaxImgSize << <BlockNum, BlockDim, 0, cstream >> >(d_LocArry, d_MaxImgWidth, d_MaxImgHigh, FluoNum, OutParaNumGS2D);


	cudaMemcpyAsync(MaxImgWidth, d_MaxImgWidth, sizeof(int), cudaMemcpyDeviceToHost, cstream);
	cudaMemcpyAsync(MaxImgHigh, d_MaxImgHigh, sizeof(int), cudaMemcpyDeviceToHost, cstream);
	cudaStreamSynchronize(cstream);


	*MaxImgWidth = (*MaxImgWidth) / 4 * 4;
	*MaxImgHigh = (*MaxImgHigh) / 4 * 4;


	cudaFree(d_MaxImgWidth);
	cudaFree(d_MaxImgHigh);

}

void ImageRenderData_TypeDef::Init(LocalizationPara & LocPara, int MaxDispImgWidth, int MaxDispImgHigh)
{
	// host and gpu
	cudaError_t err;

	
	cudaMallocHost((void **)&h_LocArry, MaxPointNum*OutParaNumGS2D*sizeof(float));
	cudaMalloc((void **)&d_LocArry, MaxPointNum*OutParaNumGS2D*sizeof(float));

	err = cudaMalloc((void **)&d_SRIntensityImg, LocPara.SRImageWidth*LocPara.SRImageHigh*sizeof(float));
	HandleErr(err, "cudaMalloc d_SRIntensityImg");

	if (RendType_Is2D(LocPara.LocType))
	{
		//2d image 
		d_SRColorMapImg = NULL;
	}
	else
	{
		// color encoded depth 3d image
		err = cudaMalloc((void **)&d_SRColorMapImg, LocPara.SRImageWidth*LocPara.SRImageHigh*sizeof(float));
		HandleErr(err, "cudaMalloc d_SRColorMapImg");
	}


	cudaMallocHost((void **)&h_MaxImageVal, sizeof(int));
	cudaMalloc((void **)&d_MaxImageVal, sizeof(int));

	cudaMallocHost((void **)&h_HistMaxDat, ImgRend_MaxDatHistLen * sizeof(int));
	cudaMalloc((void **)&d_HistMaxDat, ImgRend_MaxDatHistLen * sizeof(int));


	err = cudaMalloc((void **)&d_DispRendImg, MaxDispImgWidth*MaxDispImgHigh * 4 * sizeof(char));
	HandleErr(err, "cudaMalloc d_DispRendImg");

	err = cudaMalloc((void **)&d_SaveRendImg, LocPara.SRImageWidth*LocPara.SRImageHigh * 4);
	HandleErr(err, "cudaMalloc d_SaveRendImg");

	/*
	err = cudaMallocHost((void **)&h_DispRendImg, MaxDispImgWidth*MaxDispImgHigh * 4 * sizeof(char));
	HandleErr(err, "cudaMalloc h_DispRendImg");

	err = cudaMallocHost((void **)&h_SaveRendImg, RawImgWidth*RawImgHigh*PixelZoom*PixelZoom * 4);
	HandleErr(err, "cudaMalloc h_SaveRendImg");

	*/

	h_DispRendImg = (char *)malloc(MaxDispImgWidth*MaxDispImgHigh * 4 * sizeof(char));
	if (h_DispRendImg == NULL)
	{
		printf("malloc h_DispRendImg error\n");

	}
	h_SaveRendImg = (char *)malloc(LocPara.SRImageWidth*LocPara.SRImageHigh * 4);
	if (h_SaveRendImg == NULL)
	{
		printf("malloc h_SaveRendImg error\n");

	}


	// initial some parameters
	ResetFillImgTop(LocPara);
	ResetFillMaxVal(0);
	tRendFluoNum = 0;
}

void ImageRenderData_TypeDef::Deinit(LocalizationPara & LocPara)
{
	cudaError_t err;

	cudaFreeHost(h_LocArry);
	cudaFree(d_LocArry);

	err = cudaFree(d_SRIntensityImg);
	HandleErr(err, "cudaFree d_SRIntensityImg");

	if (RendType_Is2D(LocPara.LocType))
	{
		//2d image 

	}
	else
	{
		if (d_SRColorMapImg != NULL)
		{
			// color encoded depth 3d image
			err = cudaFree(d_SRColorMapImg);
			HandleErr(err, "cudaFree d_SRColorMapImg");

		}

	}


	cudaFreeHost(h_MaxImageVal);
	cudaFree(d_MaxImageVal);
	cudaFreeHost(h_HistMaxDat);
	cudaFree(d_HistMaxDat);

	cudaFree(d_DispRendImg);
	err = cudaFree(d_SaveRendImg);

	HandleErr(err, "cudaFree d_SaveRendImg");

//	cudaFreeHost(h_DispRendImg);
//	cudaFreeHost(h_SaveRendImg);

	free(h_DispRendImg);
	free(h_SaveRendImg);

}



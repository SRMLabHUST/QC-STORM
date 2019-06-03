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

// image render



#define RendROIWidth			7
#define RendROIWidth_half		((int)RendROIWidth/2)

__device__ static void CreateRenderPSF(float oPSFArry[RendROIWidth*RendROIWidth], float RendXPos, float RendYPos, float RendSigma)
{
	float(*RenderPSF)[RendROIWidth] = (float(*)[RendROIWidth])oPSFArry;
	int xcnt, ycnt;

	RendSigma = -0.5f / (RendSigma*RendSigma);

	RendXPos = RendXPos - (int)RendXPos; // shift to 0-0.99
	RendYPos = RendYPos - (int)RendYPos; // shift to 0-0.99
	RendXPos = RendXPos + RendROIWidth_half; // 0 is center
	RendYPos = RendYPos + RendROIWidth_half; // 0 is center

#pragma unroll
	for (ycnt = 0; ycnt < RendROIWidth; ycnt++)
	{
#pragma unroll
		for (xcnt = 0; xcnt < RendROIWidth; xcnt++)
		{
			RenderPSF[ycnt][xcnt] = __expf(((xcnt - RendXPos)*(xcnt - RendXPos) + (ycnt - RendYPos)*(ycnt - RendYPos)) * RendSigma);

		}
	}
}

/*
RGBImageEncodeMode 0 : peak photon as weight for each molecule, rendered by localization precision calculated by CRLB
RGBImageEncodeMode 1 : 1 as weight for each molecule, rendered by localization precision calculated by CRLB
RGBImageEncodeMode 2 : 1 as weight for each molecule, rendered by fixed localization precision
*/

__global__ void FluoRenderWithPrec_2D(float *d_LocArry, float *d_SRIntensityImg, int *d_MaxImageVal, int *d_HistMaxDat, int RenderingMode, float FixedlocPrec, float SNR_th, float PixelSize, float PixelZoom, int SRImageWidth, int SRImageHigh, int FluoNum)
{

	//    int tid = threadIdx.x;
	int gid = threadIdx.x + blockDim.x*blockIdx.x;

	// only for 2d
	float PSFArry[RendROIWidth*RendROIWidth];
	float(*RenderPSF)[RendROIWidth] = (float(*)[RendROIWidth])PSFArry;

	// parameters for 2d and 3d
	float(*pLocArry)[OutParaNumGS2D]; 
	pLocArry = (float(*)[OutParaNumGS2D])d_LocArry;



	int HistPos = 0;
	int Offset;


	if (gid < FluoNum)
	{
		float PeakPhoton = pLocArry[gid][Pos_PPho]; // peak photon intensity as weight 
		float RendXPos = pLocArry[gid][Pos_XPos] * PixelZoom;
		float RendYPos = pLocArry[gid][Pos_YPos] * PixelZoom;

		float CurSNR = pLocArry[gid][Pos_PSNR];
		float LocPrec = pLocArry[gid][Pos_CrbX]; // unit is nm

		CurSNR = isnan(LocPrec) ? 0 : CurSNR;
		PeakPhoton = isnan(LocPrec) ? 0 : PeakPhoton;
		PeakPhoton = (LocPrec <= 0) ? 0 : PeakPhoton;

		// valid fluos
		if ((CurSNR >= SNR_th) && (PeakPhoton > 1.0f) && (RendXPos > RendROIWidth) && (RendXPos < SRImageWidth - 1 - RendROIWidth) && (RendYPos > RendROIWidth) && (RendYPos < SRImageHigh - 1 - RendROIWidth))
		{

			if ((RenderingMode == RenderingMode_1Photon_CalculatedLocPrec) || (RenderingMode == RenderingMode_1Photon_FixedLocPrec))
			{
				// use the same weight but not the peak photon
				PeakPhoton = 1.0f;
			}
			if (RenderingMode == RenderingMode_1Photon_FixedLocPrec)
			{
				LocPrec = FixedlocPrec;
			}


			// localization precision is localization error
			if (LocPrec > 30.0f)
			{
				LocPrec = 30.0f; // set maximum resolution
				PeakPhoton = 0.0f; // do not render it
			}

			// rendered psf don't be too small
			LocPrec = max(LocPrec, PixelSize / PixelZoom / 2);


			// calculate rendering PSF
			float RendSigma = LocPrec / (PixelSize / PixelZoom);

			CreateRenderPSF(PSFArry, RendXPos, RendYPos, RendSigma);

			if (RenderingMode == RenderingMode_FittedPhoton_1Pixel)
			{
				// only render one pixel
				Offset = SRImageWidth*((int)RendYPos) + (int)RendXPos;
				atomicAdd(&d_SRIntensityImg[Offset], PeakPhoton);

			}
			else
			{
#pragma unroll
				for (int rcnt = 0; rcnt < RendROIWidth; rcnt++)
				{
					Offset = SRImageWidth*((int)RendYPos - RendROIWidth_half + rcnt) + (int)RendXPos - RendROIWidth_half;

#pragma unroll
					for (int ccnt = 0; ccnt < RendROIWidth; ccnt++)
					{
						atomicAdd(&d_SRIntensityImg[Offset + ccnt], PeakPhoton*RenderPSF[rcnt][ccnt]);

					}
				}
			}


			// for display automatic intensity adjustment
			Offset = SRImageWidth*(int)RendYPos + (int)RendXPos;
			atomicMax(d_MaxImageVal, (int)d_SRIntensityImg[Offset]);

			HistPos = (int)d_SRIntensityImg[Offset] * ImgRend_MaxDatHistLen / (*d_MaxImageVal);
			if ((HistPos >= 0) && (HistPos < ImgRend_MaxDatHistLen))
			{
				atomicAdd(&d_HistMaxDat[HistPos], 1);
			}

		}
	}
}

// render frags, mean several lines to avoid conflict with GPU localization,
// since rendering takes longer time

__global__ void RenderDisplayImg_2D(float *d_SRIntensityImg, char* d_RendImg, int *d_MaxImageVal, float LightRatio, int oImgWidth, int oImgHigh, int cposX, int cposY, int PixelStep, float PixelZoom, int SRImageWidth, int SRImageHigh)
{
	// 5 10 25 50 80 100
	//    int tid = threadIdx.x;
	int gid = threadIdx.x + blockDim.x*blockIdx.x;


	int StartPosX = cposX*PixelZoom - (oImgWidth / 2)*PixelStep;
	int StartPosY = cposY*PixelZoom - (oImgHigh / 2)*PixelStep;

	int OImgOffset = 0;
	int OImgOffset1 = 0;

	int curx;
	int cury;
	int xcnt, ycnt;

	int PixelStep1 = 1;

	PixelStep1 = PixelStep / 4;
	if (PixelStep1 < 1)PixelStep1 = 1;


	float CurIntensity = 0;
	float CurZDepth = 0;
	float tIntensity = 0;

	int cnt;
	int addcount = 0;
	int datr, datg, datb;

	for (cnt = 0; cnt < oImgHigh; cnt++)
	{
		curx = StartPosX + PixelStep*gid;
		cury = StartPosY + PixelStep*cnt;

		OImgOffset = (cnt*oImgWidth + gid) * 3;
		CurIntensity = 0;

		if ((curx >= 0) && (curx < SRImageWidth - 1 - PixelStep) && (cury >= 0) && (cury < SRImageHigh - 1 - PixelStep))
		{
			addcount = 0;
			for (ycnt = 0; ycnt < PixelStep; ycnt += PixelStep1)
			{
				OImgOffset1 = SRImageWidth*(cury + ycnt);
				for (xcnt = 0; xcnt < PixelStep; xcnt += PixelStep1)
				{
					tIntensity = d_SRIntensityImg[OImgOffset1 + (curx + xcnt)];
					CurIntensity += tIntensity;

				}
				addcount++;
			}
			//			if(addcount>1)addcount=addcount/2;

			// color temparature method, a little higher zoom ratio
			CurIntensity = CurIntensity / addcount / sqrtf(addcount);

			CurIntensity = CurIntensity*900.0f / ((*d_MaxImageVal)*LightRatio);// CurIntensity*pi/512
			if (CurIntensity > 765.0f)CurIntensity = 765.0f;

			if (CurIntensity > 510)
			{
				datr = 255;
				datg = 255;
				datb = CurIntensity - 510;
			}
			else if (CurIntensity > 255)
			{
				datr = 255;
				datg = CurIntensity - 255;
				datb = 0;
			}
			else
			{
				datr = CurIntensity;
				datg = 0;
				datb = 0;
			}
		}
		else
		{
			datr = 0;
			datg = 0;
			datb = 0;

		}
		d_RendImg[OImgOffset + 0] = datb; //0;			//b for 24 color bitmap
		d_RendImg[OImgOffset + 1] = datg; //CurIntensity;	//g for 24 color bitmap
		d_RendImg[OImgOffset + 2] = datr; //CurIntensity;	//r for 24 color bitmap
	}
}

// note the bitmap image of RGB storage for R, G, B is not 2,1,0

__global__ void RenderSaveImg_2D(float *d_SRIntensityImg, char* d_RendImg, int *d_MaxImageVal, float LightRatio, int RGBImageEncodeMode, int SRImageWidth, int SRImageHigh)
{
	// 5 10 25 50 80 100
	//    int tid = threadIdx.x;
	int gid = threadIdx.x + blockDim.x*blockIdx.x;


	int datr, datg, datb;


	int curx = gid;

	if (curx < SRImageWidth)
	{
		for (int cury = 0; cury < SRImageHigh; cury++)
		{

			// filt img Offset pos
			int FillImgOffset = SRImageWidth*cury + curx;
			float CurIntensity = d_SRIntensityImg[FillImgOffset];


			int RendImgOffset = 0;

			// result img Offset pos
			if (RGBImageEncodeMode == RGBImage_EncodeMode_3Bytes)
			{
				RendImgOffset = (cury * SRImageWidth + curx)* 3;
			}
			else
			{
				RendImgOffset = (cury * SRImageWidth + curx) * 4;
			}

			if (CurIntensity == 0)
			{
				if (RGBImageEncodeMode == RGBImage_EncodeMode_3Bytes)
				{
					d_RendImg[RendImgOffset + 0] = 0; //0;			//b for 24 color bitmap
					d_RendImg[RendImgOffset + 1] = 0; //CurIntensity;	//g for 24 color bitmap
					d_RendImg[RendImgOffset + 2] = 0; //CurIntensity;	//r for 24 color bitmap
				}
				else
				{
					d_RendImg[RendImgOffset + 0] = 0; //0;			//b for 24 color bitmap
					d_RendImg[RendImgOffset + 1] = 0; //CurIntensity;	//g for 24 color bitmap
					d_RendImg[RendImgOffset + 2] = 0; //CurIntensity;	//r for 24 color bitmap
					d_RendImg[RendImgOffset + 3] = 0xff; // alpha
				}
				continue;
			}


			// color temparature method, a little higher zoom ratio
			CurIntensity = CurIntensity*900.0f / ((*d_MaxImageVal)*LightRatio);// CurIntensity*pi/512
			if (CurIntensity>765.0f)CurIntensity = 765.0f;

			if (CurIntensity > 510)
			{
				datr = 255;
				datg = 255;
				datb = CurIntensity - 510;
			}
			else if (CurIntensity > 255)
			{
				datr = 255;
				datg = CurIntensity - 255;
				datb = 0;
			}
			else
			{
				datr = CurIntensity;
				datg = 0;
				datb = 0;
			}
			if (RGBImageEncodeMode == RGBImage_EncodeMode_3Bytes)
			{
				d_RendImg[RendImgOffset + 0] = datb; //0;			//b for 24 color bitmap
				d_RendImg[RendImgOffset + 1] = datg; //CurIntensity;	//g for 24 color bitmap
				d_RendImg[RendImgOffset + 2] = datr; //CurIntensity;	//r for 24 color bitmap
			}
			else
			{
				d_RendImg[RendImgOffset + 0] = datb; //0;			//b for 24 color bitmap
				d_RendImg[RendImgOffset + 1] = datg; //CurIntensity;	//g for 24 color bitmap
				d_RendImg[RendImgOffset + 2] = datr; //CurIntensity;	//r for 24 color bitmap
				d_RendImg[RendImgOffset + 3] = 0xff; //CurIntensity;	//r for 24 color bitmap
			}
		}

	}

}


//////////////////////////////////////////////////////////////////////////////////////////////

// code wrapper
void ImgRender_FluoRenderWithPrec_2D(float *h_LocArry, ImageRenderData_TypeDef *h_RendData, LocalizationPara & LocPara, int RenderingMode, float FixedlocPrec, int FluoNum, cudaStream_t cstream)
{

	int BlockDim = ThreadsPerBlock;
	int BlockNum = ((FluoNum + ThreadsPerBlock - 1) / ThreadsPerBlock);

	FluoRenderWithPrec_2D << <BlockNum, BlockDim, 0, cstream >> >(h_RendData->d_LocArry, h_RendData->d_SRIntensityImg, h_RendData->d_MaxImageVal, h_RendData->d_HistMaxDat, RenderingMode, FixedlocPrec, LocPara.SNR_th, LocPara.PixelSize, LocPara.PixelZoom, LocPara.SRImageWidth, LocPara.SRImageHigh, FluoNum);
	

	cudaStreamSynchronize(cstream);
}



void ImgRender_GetDispImg_2D(ImageRenderData_TypeDef *h_RendData, float BrightRatio, int oImgWidth, int oImgHigh, int cposX, int cposY, float DispZoom, float PixelZoom, int SRImageWidth, int SRImageHigh, cudaStream_t cstream)
{

	int PixelStep = (int)(1.0f / DispZoom + 0.5f);
	float LightRatio = GetLightRatio(h_RendData->h_HistMaxDat, h_RendData->d_HistMaxDat);
	LightRatio = LightRatio*BrightRatio;


	int BlockDim = ThreadsPerBlock;
	int BlockNum = ((oImgWidth + ThreadsPerBlock - 1) / ThreadsPerBlock);

	RenderDisplayImg_2D << <BlockNum, BlockDim, 0, cstream >> >(h_RendData->d_SRIntensityImg, h_RendData->d_DispRendImg, h_RendData->d_MaxImageVal, LightRatio, oImgWidth, oImgHigh, cposX, cposY, PixelStep, PixelZoom, SRImageWidth, SRImageHigh);

	cudaMemcpyAsync(h_RendData->h_DispRendImg, h_RendData->d_DispRendImg, oImgWidth*oImgHigh * 3 * sizeof(char), cudaMemcpyDeviceToHost, cstream);
	cudaStreamSynchronize(cstream);

}



void ImgRender_GetSaveImg_2D(ImageRenderData_TypeDef *h_RendData, float BrightRatio, int RGBImageEncodeMode, int SRImageWidth, int SRImageHigh, cudaStream_t cstream)
{
#if ImageRender_Debug
	cudaError_t err;
#endif

	float LightRatio;
	LightRatio = GetLightRatio(h_RendData->h_HistMaxDat, h_RendData->d_HistMaxDat);
	LightRatio = LightRatio*BrightRatio;

	int BlockDim = ThreadsPerBlock;
	int BlockNum = ((SRImageWidth + ThreadsPerBlock - 1) / ThreadsPerBlock);

	RenderSaveImg_2D << <BlockNum, BlockDim, 0, cstream >> >(h_RendData->d_SRIntensityImg, h_RendData->d_SaveRendImg, h_RendData->d_MaxImageVal, LightRatio, RGBImageEncodeMode, SRImageWidth, SRImageHigh);



	cudaMemcpyAsync(h_RendData->h_SaveRendImg, h_RendData->d_SaveRendImg, SRImageWidth*SRImageHigh * 4 * sizeof(char), cudaMemcpyDeviceToHost, cstream);
	cudaStreamSynchronize(cstream);


#if ImageRender_Debug
	HandleErr(err, "get result rend img");
#endif
}




void ImgRender_GetFillMaxVal(int *h_MaxImageVal, int *d_MaxImageVal, cudaStream_t cstream)
{
	cudaMemcpyAsync(h_MaxImageVal, d_MaxImageVal, sizeof(int), cudaMemcpyDeviceToHost, cstream);
	cudaStreamSynchronize(cstream);

}

void ImgRender_GetHistDat(int *h_HistMaxDat, int *d_HistMaxDat, cudaStream_t cstream)
{
	cudaMemcpyAsync(h_HistMaxDat, d_HistMaxDat, ImgRend_MaxDatHistLen*sizeof(int), cudaMemcpyDeviceToHost, cstream);
	cudaStreamSynchronize(cstream);
}

float GetLightRatio(int *h_HistMaxDat, int *d_HistMaxDat)
{
	cudaStream_t cstream;
	cudaStreamCreate(&cstream);

	int cnt = 0;
	float TotalDat = 0;
	float CurDat = 0.0f;
	float LightRatio = 0.1f;

	ImgRender_GetHistDat(h_HistMaxDat, d_HistMaxDat, cstream);

	if (h_HistMaxDat[1] < 1) h_HistMaxDat[1] = 1;

	for (cnt = 0; cnt < ImgRend_MaxDatHistLen - 1; cnt++) // - 1
	{
		TotalDat += h_HistMaxDat[cnt];
	}
	if (TotalDat < 1)TotalDat = 1;
	for (cnt = 0; cnt < ImgRend_MaxDatHistLen - 1; cnt++)
	{
		CurDat += h_HistMaxDat[cnt];
		if (CurDat / TotalDat > 0.95f) //98
		{
			break;
		}
	}

	if (cnt < 2)cnt = 2;
	LightRatio = (float)cnt / (float)ImgRend_MaxDatHistLen;
	
	/*
	printf("LightRatio:%d-%f\n", cnt, LightRatio);
	for (cnt = 0; cnt < ImgRend_MaxDatHistLen; cnt++) // - 1
	{
		printf("%d ", h_HistMaxDat[cnt]);
	}
*/
	cudaStreamSynchronize(cstream);
	cudaStreamDestroy(cstream);
	return LightRatio;
}




void ImgRender_ResetFillImg_2D(ImageRenderData_TypeDef *h_RendData, int SRImageWidth, int SRImageHigh)
{
	cudaStream_t cstream;
	cudaStreamCreate(&cstream);


	cudaMemsetAsync(h_RendData->d_SRIntensityImg, 0, SRImageWidth*SRImageHigh*sizeof(float), cstream);
	cudaMemsetAsync(h_RendData->d_MaxImageVal, 0, sizeof(int), cstream);
	cudaMemsetAsync(h_RendData->d_HistMaxDat, 0, ImgRend_MaxDatHistLen*sizeof(int), cstream);

	cudaStreamSynchronize(cstream);
	cudaStreamDestroy(cstream);

}


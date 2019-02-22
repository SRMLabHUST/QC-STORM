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

#include "ImageRender_3d.h"

// image render


#define RendROIWidth			7
#define RendROIWidth_half		((int)RendROIWidth/2)

__device__ void HSVtoRGB(int *r, int *g, int *b, int h, int s, int v);
__device__ void ColorGenerate(float CurZDepth, float CurIntensity, float MinZDepth, float MaxZDepth, float MaxIntensity, int ColorMode_3D, int *r, int *g, int *b);


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

__global__ void FluoRenderWithPrec_3D(float *d_LocArry, float *d_SRIntensityImg, float *d_SRColorMapImg, int *d_MaxImageVal, int *d_HistMaxDat, int RenderingMode, float FixedlocPrec, float SNR_th, float PixelSize, float PixelZoom, int SRImageWidth, int SRImageHigh, int FluoNum)
{
	//    int tid = threadIdx.x;
	int gid = threadIdx.x + blockDim.x*blockIdx.x;

	float PSFArry[RendROIWidth*RendROIWidth];
	float(*RenderPSF)[RendROIWidth] = (float(*)[RendROIWidth])PSFArry;


	// parameters for 2d and 3d
	float(*pLocArry)[OutParaNumGS2D];
	pLocArry = (float(*)[OutParaNumGS2D])d_LocArry;


	int HistPos = 0;
	int Offset;


	if (gid<FluoNum)
	{
		float PeakPhoton = pLocArry[gid][Pos_PPho]; // peak photon intensity as weight 
		float RendXPos = pLocArry[gid][Pos_XPos] * PixelZoom;
		float RendYPos = pLocArry[gid][Pos_YPos] * PixelZoom;
		float RendZPos = pLocArry[gid][Pos_ZPos]; // nm, for color

		float CurSNR = pLocArry[gid][Pos_PSNR];
		float LocPrecX = pLocArry[gid][Pos_CrbX]; // unit is nm
		float LocPrecY = pLocArry[gid][Pos_CrbY]; // unit is nm
		float LocPrec = (LocPrecX + LocPrecY) / 2;

		CurSNR = isnan(LocPrec) ? 0 : CurSNR;

		// valid fluos
		if ((CurSNR >= SNR_th) && (PeakPhoton > 1.0f) && (RendXPos > RendROIWidth) && (RendXPos < SRImageWidth - RendROIWidth) && (RendYPos > RendROIWidth) && (RendYPos < SRImageHigh - RendROIWidth))
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

			// localization precision
			if ((LocPrec > 0.01f) && (LocPrec < 40.0f))
			{				
				// avoid NAN precision
			}
			else
			{
				LocPrec = 40.0f; // set maximum 30nm resolution
				PeakPhoton = 0.0f; // do not render it
			}

			// rendered psf don't be too small
			LocPrec = max(LocPrec, PixelSize / PixelZoom / 2);

			// calculate rendering PSF
			float RendSigma = LocPrec / (PixelSize / PixelZoom);

			CreateRenderPSF(PSFArry, RendXPos, RendYPos, RendSigma);


#pragma unroll
			for (int rcnt = 0; rcnt < RendROIWidth; rcnt++)
			{
				Offset = SRImageWidth*((int)RendYPos - RendROIWidth_half + rcnt) + (int)RendXPos - RendROIWidth_half;

#pragma unroll
				for (int ccnt = 0; ccnt < RendROIWidth; ccnt++)
				{
					// light intensity weighted depth (L1*D1 + L2*D2)/(L1 + L2) for d_SRColorMapImg
					atomicAdd(&d_SRIntensityImg[Offset + ccnt], PeakPhoton*RenderPSF[rcnt][ccnt]);
					atomicAdd(&d_SRColorMapImg[Offset + ccnt], RendZPos*PeakPhoton*RenderPSF[rcnt][ccnt]);
				}
			}

			// for display automatic intensity adjustment
			Offset = SRImageWidth*(int)RendYPos + (int)RendXPos;
			atomicMax(d_MaxImageVal, (int)d_SRIntensityImg[Offset]);

			// get the relative intensity histgram, the histgram region is dynamic from 0 to d_MaxImageVal
			HistPos = (int)d_SRIntensityImg[Offset] * ImgRend_MaxDatHistLen / (*d_MaxImageVal);
			if ((HistPos >= 0) && (HistPos < ImgRend_MaxDatHistLen))
			{
				atomicAdd(&d_HistMaxDat[HistPos], 1);
			}
		}
	}
}


__global__ void RenderDisplayImg_3d(float *d_SRIntensityImg, float *d_SRColorMapImg, char* d_RendImg, int *d_MaxImageVal, float MinZDepth, float MaxZDepth, int ColorMode_3D, float LightRatio, int oImgWidth, int oImgHigh, int cposX, int cposY, int PixelStep, float PixelZoom, int SRImageWidth, int SRImageHigh)
{
	// 5 10 25 50 80 100
	//    int tid = threadIdx.x;
	int gid = threadIdx.x + blockDim.x*blockIdx.x;

	int MaxIntensity = (*d_MaxImageVal)*LightRatio;


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


	float CurIntensity=0;
	float CurZDepth=0;

	float tIntensity=0;
	float tDepth=0;

	int cnt;
	int addcount = 0;

	int datr, datg, datb;

	for (cnt = 0; cnt<oImgHigh; cnt++)
	{
		curx = StartPosX + PixelStep*gid;
		cury = StartPosY + PixelStep*cnt;

		OImgOffset = (cnt*oImgWidth + gid) * 3;

		CurIntensity = 0;
		CurZDepth = 0;

		if ((curx >= 0) && (curx < SRImageWidth - PixelStep) && (cury >= 0) && (cury < SRImageHigh - PixelStep))
		{
			addcount = 0;
			for (ycnt = 0; ycnt < PixelStep; ycnt += PixelStep1)
			{
				OImgOffset1 = SRImageWidth*(cury + ycnt);
				for (xcnt = 0; xcnt < PixelStep; xcnt += PixelStep1)
				{
					tIntensity = d_SRIntensityImg[OImgOffset1 + (curx + xcnt)];
					CurIntensity += tIntensity;

					if (tIntensity > 0)tDepth = d_SRColorMapImg[OImgOffset1 + (curx + xcnt)] / tIntensity;
					else tDepth = 0;

					CurZDepth += tDepth;
					addcount++;

				}
			}

			CurIntensity /= addcount;
			CurZDepth /=  addcount;
		}
		else
		{
			CurIntensity = 0;
			CurZDepth = 0;
		}

		ColorGenerate(CurZDepth, CurIntensity, MinZDepth, MaxZDepth, MaxIntensity, ColorMode_3D, &datr, &datg, &datb);


		d_RendImg[OImgOffset + 0] = datb; //0;			//b for 24 color bitmap
		d_RendImg[OImgOffset + 1] = datg; //CurIntensity;	//g for 24 color bitmap
		d_RendImg[OImgOffset + 2] = datr; //CurIntensity;	//r for 24 color bitmap
	}
}

__global__ void RenderSaveImg_3D(float *d_SRIntensityImg, float *d_SRColorMapImg, char* d_RendImg, int *d_MaxImageVal, float MinZDepth, float MaxZDepth, int ColorMode_3D, float LightRatio, int RGBImageEncodeMode, int SRImageWidth, int SRImageHigh)
{
	//    int tid = threadIdx.x;
	int gid = threadIdx.x + blockDim.x*blockIdx.x;

	int MaxIntensity = (*d_MaxImageVal)*LightRatio;


	int FillImgOffset = 0;
	int RendImgOffset = 0;

	int curx;
	int cury;


	int datr, datg, datb;

	float CurIntensity;
	float CurZDepth;


	int cnt;

	curx = gid;

	for (cnt = 0; cnt< SRImageHigh; cnt++)
	{
		cury = cnt;

		// filt img Offset pos
		FillImgOffset = SRImageWidth*cury + curx;
		CurIntensity = d_SRIntensityImg[FillImgOffset];


		// result img Offset pos
		if (RGBImageEncodeMode == RGBImage_EncodeMode_3Bytes)
		{
			RendImgOffset = cnt * SRImageWidth * 3 + gid * 3;
		}
		else
		{
			RendImgOffset = cnt * SRImageWidth * 4 + gid * 4;
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
				d_RendImg[RendImgOffset + 3] = 0xff; //CurIntensity;	//r for 24 color bitmap
			}
			continue;
		}

		CurZDepth = d_SRColorMapImg[FillImgOffset] / d_SRIntensityImg[FillImgOffset]; // in nm

		ColorGenerate(CurZDepth, CurIntensity, MinZDepth, MaxZDepth, MaxIntensity, ColorMode_3D, &datr, &datg, &datb);


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


//////////////////////////////////////////////////////////////////////////////////////////////

// code wrapper

void ImgRender_FluoRender_3D(float *h_LocArry, ImageRenderData_TypeDef *h_RendData, LocalizationPara & LocPara, int RenderingMode, float FixedlocPrec, int FluoNum, cudaStream_t cstream)
{

	int BlockDim = ThreadsPerBlock;
	int BlockNum = ((FluoNum + ThreadsPerBlock - 1) / ThreadsPerBlock);

	FluoRenderWithPrec_3D << <BlockNum, BlockDim, 0, cstream >> >(h_RendData->d_LocArry, h_RendData->d_SRIntensityImg, h_RendData->d_SRColorMapImg, h_RendData->d_MaxImageVal, h_RendData->d_HistMaxDat, RenderingMode, FixedlocPrec, LocPara.SNR_th, LocPara.PixelSize, LocPara.PixelZoom, LocPara.SRImageWidth, LocPara.SRImageHigh, FluoNum);


	cudaStreamSynchronize(cstream);
}



void ImgRender_GetDispImg_3D(ImageRenderData_TypeDef *h_RendData, float MinZDepth, float MaxZDepth, int ColorMode_3D, float BrightRatio, int oImgWidth, int oImgHigh, int cposX, int cposY, float DispZoom, float PixelZoom, int SRImageWidth, int SRImageHigh, cudaStream_t cstream)
{
#if ImageRender_Debug
	cudaError_t err;
#endif

	float LightRatio = GetLightRatio(h_RendData->h_HistMaxDat, h_RendData->d_HistMaxDat);
	LightRatio = LightRatio*BrightRatio;


	int PixelStep = (int)(1.0f / DispZoom + 0.5f);

	int BlockDim = ThreadsPerBlock;
	int BlockNum = ((oImgWidth + ThreadsPerBlock - 1) / ThreadsPerBlock);

	RenderDisplayImg_3d << <BlockNum, BlockDim, 0, cstream >> >(h_RendData->d_SRIntensityImg, h_RendData->d_SRColorMapImg, h_RendData->d_DispRendImg, h_RendData->d_MaxImageVal, MinZDepth, MaxZDepth, ColorMode_3D, LightRatio, oImgWidth, oImgHigh, cposX, cposY, PixelStep, PixelZoom, SRImageWidth, SRImageHigh);


	cudaMemcpyAsync(h_RendData->h_DispRendImg, h_RendData->d_DispRendImg, oImgWidth*oImgHigh * 3 * sizeof(char), cudaMemcpyDeviceToHost, cstream);
	cudaStreamSynchronize(cstream);

#if ImageRender_Debug
	HandleErr(err, "memcpy h_MaxImageVal");
#endif

}

void ImgRender_GetSaveImg_3D(ImageRenderData_TypeDef *h_RendData, float MinZDepth, float MaxZDepth, int ColorMode_3D, float BrightRatio, int RGBImageEncodeMode, int SRImageWidth, int SRImageHigh, cudaStream_t cstream)
{
#if ImageRender_Debug
	cudaError_t err;
#endif


	float LightRatio;
	LightRatio = GetLightRatio(h_RendData->h_HistMaxDat, h_RendData->d_HistMaxDat);
	LightRatio = LightRatio*BrightRatio;

	int BlockDim = ThreadsPerBlock;
	int BlockNum = (SRImageWidth + ThreadsPerBlock - 1) / ThreadsPerBlock;

	RenderSaveImg_3D << <BlockNum, BlockDim, 0, cstream >> >(h_RendData->d_SRIntensityImg, h_RendData->d_SRColorMapImg, h_RendData->d_SaveRendImg, h_RendData->d_MaxImageVal, MinZDepth, MaxZDepth, ColorMode_3D, LightRatio, RGBImageEncodeMode, SRImageWidth, SRImageHigh);


	cudaMemcpyAsync(h_RendData->h_SaveRendImg, h_RendData->d_SaveRendImg, SRImageWidth*SRImageHigh * 4 * sizeof(char), cudaMemcpyDeviceToHost, cstream);
	cudaStreamSynchronize(cstream);

#if ImageRender_Debug
	HandleErr(err, "get result rend img");
#endif
}


void ImgRender_ResetFillImg_3D(ImageRenderData_TypeDef *h_RendData, int SRImageWidth, int SRImageHigh)
{
	cudaStream_t cstream;
	cudaStreamCreate(&cstream);

	cudaMemsetAsync(h_RendData->d_SRIntensityImg, 0, SRImageWidth*SRImageHigh*sizeof(float), cstream);
	cudaMemsetAsync(h_RendData->d_SRColorMapImg, 0, SRImageWidth*SRImageHigh*sizeof(float), cstream);
	cudaMemsetAsync(h_RendData->d_MaxImageVal, 0, sizeof(int), cstream);
	cudaMemsetAsync(h_RendData->d_HistMaxDat, 0, ImgRend_MaxDatHistLen*sizeof(int), cstream);

	cudaStreamSynchronize(cstream);
	cudaStreamDestroy(cstream);
}

__device__ void ColorGenerate(float CurZDepth, float CurIntensity,float MinZDepth, float MaxZDepth, float MaxIntensity, int ColorMode_3D, int *r, int *g, int *b)
{
	int hval, sval, vval;

#define MaxVVal		100.0f

	sval = 100; // 100 saturation

	float MaxZRange = MaxZDepth - MinZDepth;

	// get the light intensity
	vval = CurIntensity * 120.0f / MaxIntensity;
	if (vval >= MaxVVal)vval = MaxVVal;

	// zdepth color rendering
	if (CurZDepth < MinZDepth)CurZDepth = MinZDepth;
	if (CurZDepth > MaxZDepth)CurZDepth = MaxZDepth;

	// 0-240
	hval = (CurZDepth - MinZDepth) * 240 / MaxZRange; // -130-0,0-130

	if (ColorMode_3D == ImgRend_ColorMode_3D_BlueToRed)
	{
		hval = 240 - hval; // HSV model, v 0 is red, 240 is blue
	}


	if (hval < 0)hval = 0;
	if (hval > 240)hval = 240;

	HSVtoRGB(r, g, b, hval, sval, vval);

}

// pure red:h=0,pure blue:h=120,pure blue:h=240
__device__ void HSVtoRGB(int *r, int *g, int *b, int h, int s, int v)
{
	// convert from HSV/HSB to RGB color
	// R,G,B from 0-255, H from 0-260, S,V from 0-100
	// ref http://colorizer.org/

	// The hue (H) of a color refers to which pure color it resembles
	// The saturation (S) of a color describes how white the color is
	// The value (V) of a color, also called its lightness, describes how dark the color is

	if (h < 0)h = 0;
	if (s < 0)s = 0;
	if (v < 0)v = 0;

	int i;

	float RGB_min, RGB_max;
	RGB_max = v*2.55f;
	RGB_min = RGB_max*(100 - s) / 100.0f;

	i = h / 60;
	int difs = h % 60; // factorial part of h

	// RGB adjustment amount by hue 
	float RGB_Adj = (RGB_max - RGB_min)*difs / 60.0f;

	switch (i) {
	case 0:
		*r = RGB_max + 0.5f; // + 0.5f for round off
		*g = RGB_min + RGB_Adj + 0.5f;
		*b = RGB_min + 0.5f;
		break;
	case 1:
		*r = RGB_max - RGB_Adj + 0.5f;
		*g = RGB_max + 0.5f;
		*b = RGB_min + 0.5f;
		break;
	case 2:
		*r = RGB_min + 0.5f;
		*g = RGB_max + 0.5f;
		*b = RGB_min + RGB_Adj + 0.5f;
		break;
	case 3:
		*r = RGB_min + 0.5f;
		*g = RGB_max - RGB_Adj + 0.5f;
		*b = RGB_max + 0.5f;
		break;
	case 4:
		*r = RGB_min + RGB_Adj + 0.5f;
		*g = RGB_min + 0.5f;
		*b = RGB_max + 0.5f;
		break;
	default:		// case 5:
		*r = RGB_max + 0.5f;
		*g = RGB_min + 0.5f;
		*b = RGB_max - RGB_Adj + 0.5f;
		break;
	}

	*r = max(*r, 0);
	*g = max(*g, 0);
	*b = max(*b, 0);


	*r = min(*r, 255);
	*g = min(*g, 255);
	*b = min(*b, 255);

}


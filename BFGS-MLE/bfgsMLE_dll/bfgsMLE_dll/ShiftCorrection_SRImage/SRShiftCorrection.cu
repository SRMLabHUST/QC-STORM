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

#include "SRShiftCorrection.h"

// image render



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


__global__ void PointFillRegionWithPrec1(float *d_LocArry, float *d_SRIntensityImg, float QE, float SNR_th, float PixelSize, float PixelZoom, int SRImageWidth, int SRImageHigh, int FluoNum)
{
	//    int tid = threadIdx.x;
	int gid = threadIdx.x + blockDim.x*blockIdx.x;

	// only for 2d
	float PSFArry[RendROIWidth*RendROIWidth];
	float(*RenderPSF)[RendROIWidth] = (float(*)[RendROIWidth])PSFArry;

	float(*pLocArry)[OutParaNumGS2D]; // for parameter array

	float RendXPos;
	float RendYPos;
	float PeakPhoton;

	float LocPrec;
	float CurSNR;
	float RendSigma;

	int rcnt = 0;
	int Offset;

	// parameters for 2d and 3d
	pLocArry = (float(*)[OutParaNumGS2D])d_LocArry;

	//PeakPhoton = 1.0f; // same weight 

	if (gid < FluoNum)
	{
		PeakPhoton = pLocArry[gid][Pos_PPho]; // peak photon intensity as weight 
		RendXPos = pLocArry[gid][Pos_XPos] * PixelZoom;
		RendYPos = pLocArry[gid][Pos_YPos] * PixelZoom;

		CurSNR = pLocArry[gid][Pos_PSNR];
		LocPrec = pLocArry[gid][Pos_CrbX]; // unit is nm


		// valid fluos
		if ((CurSNR >= SNR_th) && (PeakPhoton > 1.0f) && (RendXPos > RendROIWidth) && (RendXPos < SRImageWidth - RendROIWidth) && (RendYPos > RendROIWidth) && (RendYPos < SRImageHigh - RendROIWidth))
		{

			// localization precision is localization error
			if ((LocPrec > 0.01f) && (LocPrec < 35.0f))
			{
				// avoid NAN precision
			}
			else
			{
				LocPrec = 35.0f; // set maximum 30nm resolution
				PeakPhoton = 0.0f; // do not render it
			}


			RendSigma = LocPrec / (PixelSize / PixelZoom);

			// calculate rendering PSF
			CreateRenderPSF(PSFArry, RendXPos, RendYPos, RendSigma);

			// use the same weight but not the peak photon
			// seems 1 is better than peak photon
			PeakPhoton = 1.0f;

#pragma unroll
			for (rcnt = 0; rcnt<RendROIWidth; rcnt++)
			{
				Offset = SRImageWidth*((int)RendYPos - RendROIWidth_half + rcnt) + (int)RendXPos - RendROIWidth_half;

				atomicAdd(&d_SRIntensityImg[Offset + 0], PeakPhoton*RenderPSF[rcnt][0]);
				atomicAdd(&d_SRIntensityImg[Offset + 1], PeakPhoton*RenderPSF[rcnt][1]);
				atomicAdd(&d_SRIntensityImg[Offset + 2], PeakPhoton*RenderPSF[rcnt][2]);
				atomicAdd(&d_SRIntensityImg[Offset + 3], PeakPhoton*RenderPSF[rcnt][3]);
				atomicAdd(&d_SRIntensityImg[Offset + 4], PeakPhoton*RenderPSF[rcnt][4]);
				atomicAdd(&d_SRIntensityImg[Offset + 5], PeakPhoton*RenderPSF[rcnt][5]);
				atomicAdd(&d_SRIntensityImg[Offset + 6], PeakPhoton*RenderPSF[rcnt][6]);
			}
		}
	}
}

__global__ void gpuApplyShiftCorr(float *d_LocArry, float *d_XFrameShift, float *d_YFrameShift, float *d_ZFrameShift, int ShiftCorrEnable, int FluoNum)
{
	int gid = threadIdx.x + blockDim.x*blockIdx.x;

	float(*pLocArry)[OutParaNumGS2D]; // for parameter array

	// parameters for 2d and 3d
	pLocArry = (float(*)[OutParaNumGS2D])d_LocArry;

	int CurFrame=0;

	if (gid < FluoNum)
	{
		CurFrame = pLocArry[gid][OutParaNumGS2D - 1];

		if (ShiftCorrEnable & 0x01)
		{
			// x,y is pixel
			pLocArry[gid][Pos_XPos] = pLocArry[gid][Pos_XPos] - d_XFrameShift[CurFrame - 1];
		}
		if (ShiftCorrEnable & 0x02)
		{
			// x,y is pixel
			pLocArry[gid][Pos_YPos] = pLocArry[gid][Pos_YPos] - d_YFrameShift[CurFrame - 1];
		}

		if (ShiftCorrEnable & 0x04)
		{
			// z is nm
//			pLocArry[gid][Pos_ZPos] = pLocArry[gid][Pos_ZPos] - d_ZFrameShift[CurFrame];
		}

	}
}

__global__ void gpuCrossCorr_Mul(float *d_FillImg1, float *d_FillImg2, float *d_MulImg, int ShiftX, int ShiftY, int CorrShiftBiasX, int CorrShiftBiasY, int SRImageWidth, int SRImageHigh)
{
	int gid = threadIdx.x + blockDim.x*blockIdx.x;
	int curx2 = gid % SRImageWidth;
	int cury2 = gid / SRImageWidth;

	int curx1 = curx2 + CorrShiftBiasX + ShiftX;
	int cury1 = cury2 + CorrShiftBiasY + ShiftY;

	if (gid < SRImageWidth*SRImageHigh)
	{
		if ((curx1 >= 0) && (curx1 < SRImageWidth) && (cury1 >= 0) && (cury1 < SRImageHigh))
		{
			d_MulImg[gid] = d_FillImg1[gid] * d_FillImg2[cury1*SRImageWidth + curx1];
		}
	}
}

__global__ void gpuCrossCorr_SumOfMul(float *d_MulImg, float *d_SumLine, int SRImageWidth, int SRImageHigh)
{
	int gid = threadIdx.x + blockDim.x*blockIdx.x;

	int ProcNum = SRImageHigh / ThreadsPerBlock;
	int ResidualNum = SRImageHigh % ThreadsPerBlock;

	int pcnt = 0;
	int ycnt = 0;
	int cury;
	int curx = gid;
	float SumDat = 0;


	if (gid < SRImageWidth)
	{
		for (pcnt = 0; pcnt < ProcNum; pcnt++)
		{
			for (ycnt = 0; ycnt < ThreadsPerBlock; ycnt++)
			{
				cury = pcnt*ThreadsPerBlock + ycnt;
				SumDat += d_MulImg[cury*SRImageWidth + curx];
			}
		}
		for (ycnt = 0; ycnt < ResidualNum; ycnt++)
		{
			cury = pcnt*ThreadsPerBlock + ycnt;
			SumDat += d_MulImg[cury*SRImageWidth + curx];
		}
		d_SumLine[curx] = SumDat;
	}
}



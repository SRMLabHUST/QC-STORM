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

#include "cuda_runtime.h"
#include "device_launch_parameters.h"


#include <stdio.h>

#include "bfgs_CommonPara.h"


#define ThreadsPerBlock						32


#define GaussianCRLB_ParaNum				6

#define CRLB_Para_TPhn						0
#define CRLB_Para_Bkgn						1
#define CRLB_Para_XPos						2
#define CRLB_Para_YPos						3
#define CRLB_Para_SigX						4
#define CRLB_Para_SigY						5

#define GaussianCRLB_ROISize2D				9
#define GaussianCRLB_ROISize3D				15


#define TakeReadNoiseToBackground			1



template<int GaussianCRLB_ROI>
__device__ void GaussianPSF_generate2D(float *pPSFDat, float* ParaArry);

template<int DatLen>
__device__ void PSFData_Log_Calc(float *pPSFDat, float *pPSFDat_log);

template<int DatLen>
__device__ void Derivative2_calc(float *Deriv2, float *dat_0, float*dat_p, float *dat_n, float dh);

template<int DatLen>
__device__ float GaussianCRLB_GetCRLB(float* PSFModelData, float *PSFDeriv2, float *PSFDeriv2_log);

template<int GaussianCRLB_ROI, int DatLen>
__device__ void GaussianCRLB_Derivative2(float *pPSFDeviv2, float *pPSFDeviv2_log, float *ParaArry, int ParaID);


///////////////////////////////////////////////////////////////////////////////////////
// for round PSF, we only calculate the CRLB of x
// x0, y0 is set to 0, to only consider the effect of photon and background, psf width etc.
template<int GaussianCRLB_ROI>
__global__ void GaussianCRLB_Calc_top(float *d_LocArry, float ReadNoise_e, float QE, float PixelSize)
{
	enum {
		GaussianCRLB_HalfROI = (int)(GaussianCRLB_ROI / 2),
		GaussianCRLB_ROIArea = (GaussianCRLB_ROI*GaussianCRLB_ROI),
	};

	float PSFModelData[GaussianCRLB_ROIArea];
	float PSFDeriv2[GaussianCRLB_ROIArea];
	float PSFDeriv2_log[GaussianCRLB_ROIArea];

	int gid = threadIdx.x + blockDim.x*blockIdx.x;
	//	int tid = threadIdx.x;

	float(*pLocArry)[OutParaNumGS2D] = (float(*)[OutParaNumGS2D])d_LocArry; // for parameter array

	float TotalPhoton_e = pLocArry[gid][Pos_TPho] * QE; //
	float Background_e = pLocArry[gid][Pos_Bakg] * QE;  //
	float PSFSigmaX = pLocArry[gid][Pos_SigX]; // pixel
	float PSFSigmaY = pLocArry[gid][Pos_SigY]; // pixel

	float ParaStdX = 100.0f;
	float ParaStdY = 100.0f;

	if ((TotalPhoton_e > 1) && (Background_e > 1) && (PSFSigmaX > 0.1) && (PSFSigmaY > 0.1) && (PSFSigmaX < 10) && (PSFSigmaY < 10))
	{
		float ParaArry[GaussianCRLB_ParaNum];

		// for simplicity, count read noise to background level
		if (TakeReadNoiseToBackground)
		{
			Background_e = (Background_e*QE + ReadNoise_e*ReadNoise_e) / QE;
		}

		ParaArry[CRLB_Para_TPhn] = TotalPhoton_e;
		ParaArry[CRLB_Para_Bkgn] = Background_e;

		ParaArry[CRLB_Para_XPos] = 0; // x0,pixel
		ParaArry[CRLB_Para_YPos] = 0; // y0,pixel
		ParaArry[CRLB_Para_SigX] = PSFSigmaX; // pixel
		ParaArry[CRLB_Para_SigY] = PSFSigmaY; // pixel


											  // Gaussian PSF
		GaussianPSF_generate2D<GaussianCRLB_ROI>(PSFModelData, ParaArry);

		// second derivative
		GaussianCRLB_Derivative2<GaussianCRLB_ROI, GaussianCRLB_ROIArea>(PSFDeriv2, PSFDeriv2_log, ParaArry, CRLB_Para_XPos);
		ParaStdX = GaussianCRLB_GetCRLB<GaussianCRLB_ROIArea>(PSFModelData, PSFDeriv2, PSFDeriv2_log);

		// second derivative
		GaussianCRLB_Derivative2<GaussianCRLB_ROI, GaussianCRLB_ROIArea>(PSFDeriv2, PSFDeriv2_log, ParaArry, CRLB_Para_YPos);
		ParaStdY = GaussianCRLB_GetCRLB<GaussianCRLB_ROIArea>(PSFModelData, PSFDeriv2, PSFDeriv2_log);

		// write to loc array
		ParaStdX = ParaStdX*PixelSize; // nm
		ParaStdY = ParaStdY*PixelSize; // nm

	}

	pLocArry[gid][Pos_CrbX] = ParaStdX;
	pLocArry[gid][Pos_CrbY] = ParaStdY;

}

template<int GaussianCRLB_ROI>
__device__ void GaussianPSF_generate2D(float *pPSFDat, float* ParaArry)
{
	enum {
		GaussianCRLB_HalfROI = (int)(GaussianCRLB_ROI / 2),
	};

	float(*pPSFROI)[GaussianCRLB_ROI] = (float(*)[GaussianCRLB_ROI])pPSFDat;


	float TotalPhoton_e = ParaArry[CRLB_Para_TPhn];
	float Background_e = ParaArry[CRLB_Para_Bkgn];

	float x0 = ParaArry[CRLB_Para_XPos];
	float y0 = ParaArry[CRLB_Para_YPos];

	float PSFSigmaX = ParaArry[CRLB_Para_SigX];
	float PSFSigmaY = ParaArry[CRLB_Para_SigY];

	float pi = 3.141592654f;
	float PeakPn = 1 / (2 * pi * PSFSigmaX * PSFSigmaY) * TotalPhoton_e;

	float PSFSigma_ix = 1 / (2 * PSFSigmaX * PSFSigmaX);
	float PSFSigma_iy = 1 / (2 * PSFSigmaY * PSFSigmaY);


	float CurX, CurY;
	for (int ycnt = 0; ycnt < GaussianCRLB_ROI; ycnt++)
	{
		CurY = ycnt - GaussianCRLB_HalfROI;
		for (int xcnt = 0; xcnt < GaussianCRLB_ROI; xcnt++)
		{
			CurX = xcnt - GaussianCRLB_HalfROI;
			pPSFROI[ycnt][xcnt] = PeakPn *__expf(-((CurX - x0)*(CurX - x0)*PSFSigma_ix + (CurY - y0)*(CurY - y0)*PSFSigma_iy)) + Background_e;
		}
	}
}

template<int DatLen>
__device__ void PSFData_Log_Calc(float *pPSFDat, float *pPSFDat_log)
{
	int cnt;
#pragma unroll
	for (cnt = 0; cnt < DatLen; cnt++)
	{
		// can't use __logf
		pPSFDat_log[cnt] = logf(pPSFDat[cnt]);
	}
}

template<int DatLen>
__device__ void Derivative2_calc(float *Deriv2, float *dat_0, float*dat_p, float *dat_n, float dh)
{
	int cnt;
#pragma unroll
	for (cnt = 0; cnt < DatLen; cnt++)
	{
		Deriv2[cnt] = (dat_p[cnt] - 2 * dat_0[cnt] + dat_n[cnt]) / dh / dh;
	}
}

template<int GaussianCRLB_ROI, int DatLen>
__device__ void GaussianCRLB_Derivative2(float *pPSFDeviv2, float *pPSFDeviv2_log, float *ParaArry, int ParaID)
{
	float tPSFData_0[DatLen];
	float tPSFData_log_0[DatLen];

	float tPSFData_p[DatLen];
	float tPSFData_log_p[DatLen];

	float tPSFData_n[DatLen];
	float tPSFData_log_n[DatLen];

	const float dh = 0.02f;

	// second deriative by numerical  method
	GaussianPSF_generate2D<GaussianCRLB_ROI>(tPSFData_0, ParaArry);
	PSFData_Log_Calc<DatLen>(tPSFData_0, tPSFData_log_0);


	ParaArry[ParaID] += dh;
	GaussianPSF_generate2D<GaussianCRLB_ROI>(tPSFData_p, ParaArry);
	PSFData_Log_Calc<DatLen>(tPSFData_p, tPSFData_log_p);


	ParaArry[ParaID] -= 2 * dh;
	GaussianPSF_generate2D<GaussianCRLB_ROI>(tPSFData_n, ParaArry);
	PSFData_Log_Calc<DatLen>(tPSFData_n, tPSFData_log_n);

	ParaArry[ParaID] += dh;

	//
	Derivative2_calc<DatLen>(pPSFDeviv2, tPSFData_0, tPSFData_p, tPSFData_n, dh);
	Derivative2_calc<DatLen>(pPSFDeviv2_log, tPSFData_log_0, tPSFData_log_p, tPSFData_log_n, dh);

}

template<int DatLen>
__device__ float GaussianCRLB_GetCRLB(float* PSFModelData, float *PSFDeriv2, float *PSFDeriv2_log)
{
	int cnt;
	float FishirInf = 0;
	for (cnt = 0; cnt < DatLen; cnt++)
	{
		FishirInf += (PSFDeriv2[cnt] - PSFModelData[cnt] * PSFDeriv2_log[cnt]);
	}

	if (FishirInf == 0.0f)FishirInf = 1;

	// standard deviation of parameter estimation
	float ParaStd = sqrtf(1 / FishirInf);

	return ParaStd;
}

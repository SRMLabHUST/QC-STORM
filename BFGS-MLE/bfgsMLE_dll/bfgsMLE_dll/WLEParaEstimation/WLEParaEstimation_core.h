#pragma once

#include "cuda_runtime.h"
#include "device_launch_parameters.h"

#include <stdio.h>


#include "WLEParaEstimation_Parameters.h"



#include "SplineInterpolation.h"


__device__ int FindMaxPos(float *iData, int DataLen, int LSel, int RSel);

__device__ void SigmaWidthCalc_TwoSides(float *SigmaL, float *SigmaR, float *MeanData_Intp, int DataLen, int MaxPos);


// find the index i where xn >= x_i[i]
__host__ __device__ int FindPositionID_search(float *x_i, float xn, int InputDataLen, float UserPara);

__host__ __device__ int FindPositionID_0(float *x_i, float xn, int InputDataLen, float InterplotGap);



template <int ROISize>
__global__ void gpu_CalculatePSFWidth(unsigned short * d_ROIMem, float *d_WLEPara, int FluoNum)
{
	// margin, valid data, margin
	enum {
		MarginLen = (ROISize == 7) ? 1 : ((ROISize < 7) ? 0 : 2),
		AvgDataLen = ROISize - 2 * MarginLen,

		InterpolatedDatalen = (AvgDataLen - 1)*InterpolationRatio + 1,
	};


	int gid = threadIdx.x + blockDim.x*blockIdx.x;
	int tid = threadIdx.x;

	const int ROIWholeSize = ROISize*(ROISize + 1);

	const int CurROIAddr = ROIWholeSize*gid;

	unsigned short(*pROI)[ROISize] = (unsigned short(*)[ROISize])&d_ROIMem[CurROIAddr];


	float(*pWLEPara)[WLE_ParaNumber] = (float(*)[WLE_ParaNumber])d_WLEPara;


	// mean value along x,y direction
	float MeanX[AvgDataLen];
	float MeanY[AvgDataLen];

	float x0[AvgDataLen];

	float xn[InterpolatedDatalen];

	float yn_MeanX[InterpolatedDatalen];
	float yn_MeanY[InterpolatedDatalen];


	// PSF width along four direction
	float SigmaL = 0;
	float SigmaR = 0;
	float SigmaU = 0;
	float SigmaD = 0;


	if (gid < FluoNum)
	{
		// get mean value along x,y direction
		for (int cnt = 0; cnt < AvgDataLen; cnt++)
		{
			MeanX[cnt] = 0;
			MeanY[cnt] = 0;
		}

		for (int r = 0; r < AvgDataLen; r++)
		{
			for (int c = 0; c < AvgDataLen; c++)
			{
				MeanX[c] = MeanX[c] + pROI[r + MarginLen][c + MarginLen];
				MeanY[r] = MeanY[r] + pROI[r + MarginLen][c + MarginLen];
			}
		}


		float MinX = MeanX[0];
		float MinY = MeanY[0];

		for (int cnt = 0; cnt < AvgDataLen; cnt++)
		{
			MinX = min(MinX, MeanX[cnt]);
			MinY = min(MinY, MeanY[cnt]);
		}

		for (int cnt = 0; cnt < AvgDataLen; cnt++)
		{
			MeanX[cnt] = (MeanX[cnt] - MinX) / AvgDataLen;
			MeanY[cnt] = (MeanY[cnt] - MinY) / AvgDataLen;
		}

		// spline interpolation
		for (int cnt = 0; cnt < AvgDataLen; cnt++)
		{
			x0[cnt] = cnt;
		}

		for (int cnt = 0; cnt < InterpolatedDatalen; cnt++)
		{
			xn[cnt] = cnt * InterpolationGap;
		}

		SplineInterpolation<AvgDataLen, InterpolatedDatalen>(x0, MeanX, xn, yn_MeanX, FindPositionID_0, 0);
		SplineInterpolation<AvgDataLen, InterpolatedDatalen>(x0, MeanY, xn, yn_MeanY, FindPositionID_0, 0);

		int FindRegionlen = 1.5f*InterpolationRatio;
		int InterpCenter = InterpolatedDatalen / 2;

		int MaxPos_X = FindMaxPos(yn_MeanX, InterpolatedDatalen, InterpCenter - FindRegionlen, InterpCenter + FindRegionlen);
		int MaxPos_Y = FindMaxPos(yn_MeanY, InterpolatedDatalen, InterpCenter - FindRegionlen, InterpCenter + FindRegionlen);


		SigmaWidthCalc_TwoSides(&SigmaL, &SigmaR, yn_MeanX, InterpolatedDatalen, MaxPos_X);
		SigmaWidthCalc_TwoSides(&SigmaU, &SigmaD, yn_MeanY, InterpolatedDatalen, MaxPos_Y);

		// use small ROI size, thus PSF width is estimated smaller
		pWLEPara[gid][WLE_Para_SigmaL] = SigmaL*1.25f;
		pWLEPara[gid][WLE_Para_SigmaR] = SigmaR*1.25f;
		pWLEPara[gid][WLE_Para_SigmaU] = SigmaU*1.25f;
		pWLEPara[gid][WLE_Para_SigmaD] = SigmaD*1.25f;

		// the high discrepency comes from low SNR fluo or high density
//		float diff1 = abs(SigmaL - SigmaR) / min(SigmaL, SigmaR);
//		float diff2 = abs(SigmaU - SigmaD) / min(SigmaU, SigmaD);

	}
}

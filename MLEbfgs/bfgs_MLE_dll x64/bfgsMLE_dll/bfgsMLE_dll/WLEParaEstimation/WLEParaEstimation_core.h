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

#pragma once

#include "cuda_runtime.h"
#include "device_launch_parameters.h"

#include <stdio.h>


#include "WLEParaEstimation_Parameters.h"

#include "LDROIExtraction_Param.h"


#include "SplineInterpolation.h"


__device__ int FindMaxPos(float *iData, int DataLen, int LSel, int RSel);

__device__ void SigmaWidthCalc_TwoSides(float *SigmaL, float *SigmaR, float *MeanData_Intp, int DataLen, int MaxPos);


// find the index i where xn >= x_i[i]
__host__ __device__ int FindPositionID_search(float *x_i, float xn, int InputDataLen, float UserPara);

__host__ __device__ int FindPositionID_0(float *x_i, float xn, int InputDataLen, float InterplotGap);



template <int ROISize>
__global__ void gpu_CalculatePSFWidth(unsigned short * d_ImageROI, float *d_WLEPara, int FluoNum)
{
	// margin, valid data, margin
	enum {
//		MarginLen = (ROISize < 7) ? 0 : 1,
		MarginLen = 0,

		AvgDataLen = ROISize - 2 * MarginLen,

		InterpolatedDatalen = (AvgDataLen - 1)*InterpolationRatio + 1,
	};


	int gid = threadIdx.x + blockDim.x*blockIdx.x;
	int tid = threadIdx.x;


	const int ROIWholeSize = ROISize*(ROISize + 1);

	const int CurROIAddr = ROIWholeSize*gid;

	unsigned short(*pROI)[ROISize] = (unsigned short(*)[ROISize])&d_ImageROI[CurROIAddr];


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

		int FindRegionlen = InterpolationRatio * PeakFindRangeRatio;
		int InterpCenter = InterpolatedDatalen / 2;

		int MaxPos_X = FindMaxPos(yn_MeanX, InterpolatedDatalen, InterpCenter - FindRegionlen, InterpCenter + FindRegionlen);
		int MaxPos_Y = FindMaxPos(yn_MeanY, InterpolatedDatalen, InterpCenter - FindRegionlen, InterpCenter + FindRegionlen);


		SigmaWidthCalc_TwoSides(&SigmaL, &SigmaR, yn_MeanX, InterpolatedDatalen, MaxPos_X);
		SigmaWidthCalc_TwoSides(&SigmaU, &SigmaD, yn_MeanY, InterpolatedDatalen, MaxPos_Y);


		pWLEPara[gid][WLE_Para_SigmaL] = SigmaL;
		pWLEPara[gid][WLE_Para_SigmaR] = SigmaR;
		pWLEPara[gid][WLE_Para_SigmaU] = SigmaU;
		pWLEPara[gid][WLE_Para_SigmaD] = SigmaD;

	}
}

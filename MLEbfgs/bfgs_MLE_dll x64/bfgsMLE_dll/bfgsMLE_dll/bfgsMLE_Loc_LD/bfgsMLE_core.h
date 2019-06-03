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



#include "bfgs_CommonPara.h"

#include "bfgs_base.h"


#include "WLEParaEstimation_Parameters.h"



#define ThreadsPerBlock			32 //Threads Per Block



// enable weighted MLE
#define WLE_ENABLE				1




// enable for fair comparasion bwteen WLE and MLE

#define SPEED_TEST				0


#define Math_PI					3.14159265358979f



// scaling factor of MLE fitting
#define AScalFactor				(128.0f)
#define BScalFactor				(64.0f)


#define rAScalFactor			(1.0f/AScalFactor)
#define rBScalFactor			(1.0f/BScalFactor)

#define SScalFactor				(0.1f)
#define rSScalFactor			(1.0f/SScalFactor)




template <int FitParaNum>
__device__ void VectorAddMul1(float oVector[], float Ininf[][ThreadsPerBlock], float d0[][ThreadsPerBlock], float coeff, int tid)
{
	int cnt = 0;
#pragma unroll
	for (cnt = 0; cnt < FitParaNum; cnt++)
	{
		oVector[cnt] = d0[cnt][tid] * coeff + Ininf[cnt][tid];
	}
}

template <int FitParaNum>
__device__ void IninfConstrain(float Ininf[][ThreadsPerBlock], int tid)
{
	int cnt = 0;
#pragma unroll
	for (cnt = 0; cnt < FitParaNum; cnt++)
	{
		if (Ininf[cnt][tid]<0.01f) Ininf[cnt][tid] = 0.01f; // A
	}
}

template <int FitParaNum>
__device__ void ConstructD0(float D0[], float sk[FitParaNum], float yk[FitParaNum], int tid)
{
	float divdat;
	float skyk[FitParaNum*FitParaNum]; //  I-(sk*yk')/(yk'*sk)
	//	float yksk[25]; // the same with skyk but transposition
	float sksk[FitParaNum*FitParaNum];
	float tD0[FitParaNum*FitParaNum];

	float(*pskyk)[FitParaNum] = (float(*)[FitParaNum])skyk;
	float(*psksk)[FitParaNum] = (float(*)[FitParaNum])sksk;
	float(*ptD0)[FitParaNum] = (float(*)[FitParaNum])tD0;
	float(*pD0)[FitParaNum] = (float(*)[FitParaNum])&D0[0];

	int row = 0;
	int col = 0;
	int cnt = 0;

	// = sk.*yk
	divdat = 0;
#pragma unroll
	for (cnt = 0; cnt < FitParaNum; cnt++)
	{
		divdat += sk[cnt] * yk[cnt];
	}

	divdat = 1.0f / divdat;
	// divdat = __fdividef(1.0f , divdat); // 

	float tdat10[FitParaNum];
	float tdat20[FitParaNum];

#pragma unroll
	for (cnt = 0; cnt < FitParaNum; cnt++)
	{
		tdat10[cnt] = yk[cnt] * divdat;
		tdat20[cnt] = sk[cnt] * divdat;
	}

#pragma unroll
	for (row = 0; row<FitParaNum; row++)
	{
#pragma unroll
		for (col = 0; col<FitParaNum; col++)
		{
			// I-(sk*yk')/(yk'*sk)
			if (row == col) pskyk[row][col] = 1.0f - sk[row] * tdat10[col];
			else         pskyk[row][col] = 0.0f - sk[row] * tdat10[col];

			// (sk*sk')/(yk'*sk)
			psksk[row][col] = sk[row] * tdat20[col];
		}
	}
	// 
#pragma unroll
	for (row = 0; row<FitParaNum; row++)
	{
#pragma unroll
		for (col = 0; col<FitParaNum; col++)
		{
			// tD0 = skyk*D0
			ptD0[row][col] = 0;
#pragma unroll
			for (cnt = 0; cnt < FitParaNum; cnt++)
			{
				ptD0[row][col] += pskyk[row][cnt] * pD0[cnt][col];
			}
		}
	}

#pragma unroll
	for (row = 0; row<FitParaNum; row++)
	{
#pragma unroll
		for (col = 0; col<FitParaNum; col++)
		{
			// D0 = D0*yksk
			pD0[row][col] = 0;
#pragma unroll
			for (cnt = 0; cnt < FitParaNum; cnt++)
			{
				pD0[row][col] += ptD0[row][cnt] * pskyk[col][cnt];
			}
		}
	}
	// D0=D0+sksk;
#pragma unroll
	for (row = 0; row<FitParaNum; row++)
	{
#pragma unroll
		for (cnt = 0; cnt < FitParaNum; cnt++)
		{
			pD0[row][cnt] = pD0[row][cnt] + psksk[row][cnt];
		}
	}
}

template <int FitParaNum>
__device__ void MatMulVector(float D0[], float grad[][ThreadsPerBlock], float d0[][ThreadsPerBlock], int tid)
{
	int row;
	int cnt = 0;

	float(*pD0)[FitParaNum] = (float(*)[FitParaNum])&D0[0];

	// d0=-D0*grad;    %search direction
#pragma unroll
	for (row = 0; row<FitParaNum; row++)
	{
		d0[row][tid] = 0;
#pragma unroll
		for (cnt = 0; cnt < FitParaNum; cnt++)
		{
			d0[row][tid] -= pD0[row][cnt] * grad[cnt][tid];
		}
	}

}

template <int ROISize>
__device__ void WLEWeightsCalc(float* WLE_Weight, float WLE_SigmaX, float WLE_SigmaY, int MoleculeType)
{
	float(*pWLE_Weight)[ROISize] = (float(*)[ROISize])WLE_Weight;

	// use pre-estimated center pos don't help improve the localization precision
	float ROICenter = ROISize / 2.0f;


	float WLE_SigmaX1 = 1.0f / (2.0f * WLE_SigmaX * WLE_SigmaX);
	float WLE_SigmaY1 = 1.0f / (2.0f * WLE_SigmaY * WLE_SigmaY);


	float rowpos, colpos;

	for (int row = 0; row < ROISize; row++)
	{
		for (int col = 0; col < ROISize; col++)
		{
			rowpos = row + 0.5f; // pixel center position
			colpos = col + 0.5f;

			if (MoleculeType == MoleculeType_MLEFit)
			{
				pWLE_Weight[row][col] = 1.0f;
			}
			else
			{
				pWLE_Weight[row][col] = __expf(-((colpos - ROICenter)*(colpos - ROICenter)*WLE_SigmaX1 + (rowpos - ROICenter)*(rowpos - ROICenter)*WLE_SigmaY1));
			}
		}
	}
}

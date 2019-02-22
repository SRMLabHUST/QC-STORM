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


#define ThreadsPerBlock			32 //Threads Per Block


// enable weighted MLE
#define WLE_ENABLE				1



#define Math_PI					3.14159265358979f

// for 2d loc


#define FitParaNum_2D			5

#define IterateNum_2D			8  // total iteration number, fixed iteration number
#define IterateNum_bs_2D		11 // bisection iteration to find best walk length

// fitting parameters gaussian 2d
#define Fit2D_Peak				0
#define Fit2D_XPos				1
#define Fit2D_YPos				2
#define Fit2D_Sigm				3
#define Fit2D_Bakg				4

// fitting parameters for astigmatism 3d
#define FitParaNum_AS3D			6

#define IterateNum_AS3D			11  // total iteration number, fixed iteration number
#define IterateNum_bs_AS3D		11 // bisection iteration to find best walk length



#define FitAS3D_Peak				0
#define FitAS3D_XPos				1
#define FitAS3D_YPos				2
#define FitAS3D_SigX				3
#define FitAS3D_SigY				4
#define FitAS3D_Bakg				5


// fitting parameters for double helix 3d
#define FitParaNum_DH3D			7

#define IterateNum_DH3D			11  // total iteration number, fixed iteration number
#define IterateNum_bs_DH3D		11 // bisection iteration to find best walk length


#define DHSystemPSF				2.0f
#define DHParaSigma				(1.0f/(2.0f*DHSystemPSF*DHSystemPSF))

#define FitDH3D_Peak				0
#define FitDH3D_XPos				1
#define FitDH3D_YPos				2
#define FitDH3D_SigX				3
#define FitDH3D_SigY				4
#define FitDH3D_Bakg				5
#define FitDH3D_Angl				6


//
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

/////////////////////////////////////////////////////////////////////
// for weighted MLE

#define PSFInt_Num				10
#define PSFInt_Gap				(1.0f/PSFInt_Num)



template <int RegionLen>
__device__ void RemoveMeanDatBkg(float *MeanDat)
{
	float MinDat = MeanDat[0];

#pragma unroll
	for (int cnt = 0; cnt < RegionLen; cnt++)
	{
		MinDat = min(MinDat, MeanDat[cnt]);
	}

#pragma unroll
	for (int cnt = 0; cnt < RegionLen; cnt++)
	{
		MeanDat[cnt] = MeanDat[cnt] - MinDat;
	}
}

template <int RegionLen>
__device__ float GetMeanDatCenterPos_Rough(float *MeanDat)
{
	float wsum = 0;
	float sum = 0;

	int ROICenter = RegionLen / 2;
	int DatSelL = ROICenter - 2;
	int DatSelR = ROICenter + 2;

	DatSelL = max(DatSelL, 0);
	DatSelR = min(DatSelR, RegionLen - 1);

#pragma unroll
	for (int cnt = DatSelL; cnt <= DatSelR; cnt++)
	{
		wsum += (cnt + 0.5f)*MeanDat[cnt];
		sum += MeanDat[cnt];
	}

	sum = max(sum, 1.0f);

	return wsum / sum;

}

template <int RegionLen>
__device__ float GetMeanDatCenterPos(float *MeanDat)
{
	float wsum = 0;
	float sum = 0;

#pragma unroll
	for (int cnt = 0; cnt < RegionLen; cnt++)
	{
		wsum += (cnt + 0.5f)*MeanDat[cnt];
		sum += MeanDat[cnt];
	}

	sum = max(sum, 1.0f);

	return wsum / sum;

}

template <int RegionLen>
__device__ void PSFInterpolation(float *MeanDat, float *PSFIntData, int IntArrayLen)
{

	int y0 = 0;
	int y1 = 0;
	float d0 = 0;
	float d1 = 0;
	for (int cnt = 0; cnt < IntArrayLen; cnt++)
	{
		y0 = cnt / PSFInt_Num;
		y1 = y0+1;
		y1 = min(y1, RegionLen - 1);

		d0 = (cnt % 10)*PSFInt_Gap;
		d1 = 1.0f - d0;

		PSFIntData[cnt] = MeanDat[y0] * d1 + MeanDat[y1] * d0;
	}

}



template <int RegionLen>
__device__ float GetPSFSigmaEstimation(float *MeanDat, float MaxPos)
{
	// for proper storage
	enum {
		IntArrayLen = (RegionLen - 1)*PSFInt_Num + 2,
	};

	float PSFIntData[IntArrayLen];

	int IntArrayLen1 = IntArrayLen - 1;

	PSFInterpolation<RegionLen>(MeanDat, PSFIntData, IntArrayLen1);

	// MaxPos is center position by default

	float SigmaL;
	float SigmaR;

	float WDat = 0;
	float SDat = 0;

	float CurPos = 0;
	float CurSel = 0;

	// left side
	WDat = 0;
	SDat = 0;


	for (int cnt = 0; cnt < IntArrayLen1; cnt++)
	{
		CurPos = 0.5f + cnt*PSFInt_Gap;
		CurSel = CurPos < (MaxPos - PSFInt_Gap);

		WDat = WDat + PSFIntData[cnt] * (MaxPos - CurPos) * CurSel;
		SDat = SDat + PSFIntData[cnt] * CurSel;
	}
	SDat = max(SDat, 1.0f);

	SigmaL = 1.414f*WDat / SDat;

	// right side
	WDat = 0;
	SDat = 0;
	for (int cnt = 0; cnt < IntArrayLen1; cnt++)
	{
		CurPos = 0.5f + cnt*PSFInt_Gap;
		CurSel = CurPos > (MaxPos + PSFInt_Gap);

		WDat = WDat + PSFIntData[cnt] * (CurPos - MaxPos)* CurSel;
		SDat = SDat + PSFIntData[cnt] * CurSel;
	}
	SDat = max(SDat, 1.0f);

	SigmaR = 1.414f* WDat / SDat;

	float SigmaDiff = abs(SigmaL - SigmaR) / min(SigmaL, SigmaR);


//	printf("esti:%.4f %.4f %.4f %.4f\n", MaxPos, SigmaL, SigmaR, SigmaDiff/ min(SigmaL, SigmaR));

	// get the minimul one to avoid signal contamination
	SigmaL = min(SigmaL, SigmaR);


#define PSFWidth_Torlation		0.20f

	float Ratio = SigmaDiff - PSFWidth_Torlation;
	Ratio = max(Ratio, 0.0f);

	// not that necessary to apply WLE
	if (SigmaDiff > PSFWidth_Torlation)
	{
		SigmaL = SigmaL / (1.414f + Ratio / 2.0f); // 
	}

	//
	if (SigmaL <= 0)SigmaL = RegionLen / 2.0f / 2.355f;
	if (SigmaL <= 1.0f)SigmaL = 1.0f;

	return SigmaL;
}

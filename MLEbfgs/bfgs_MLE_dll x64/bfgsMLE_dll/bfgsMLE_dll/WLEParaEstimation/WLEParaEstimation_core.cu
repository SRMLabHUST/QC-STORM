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

#include "WLEParaEstimation.h"



// find the index i where xn >= x_i[i]
__host__ __device__ int FindPositionID_0(float *x_i, float xn, int InputDataLen, float InterplotGap)
{
#define X0_i_FirstPixel		0.0f

	int i = xn - X0_i_FirstPixel;

	return i;
}


// find the index i where xn >= x_i[i]
__host__ __device__ int FindPositionID_search(float *x_i, float xn, int InputDataLen, float UserPara)
{
	int i;
	for (i = 0; i < InputDataLen; i++)
	{
		if (xn >= x_i[i])
		{
			break;
		}
	}

	return i;
}


void CalculatePSFWidth(unsigned short * d_ImageROI, float *d_WLEPara, int FluoNum, int ROISize, cudaStream_t cstream)
{

	int BlockDim = ThreadsPerBlock;
	int BlockNum = (FluoNum + ThreadsPerBlock - 1) / ThreadsPerBlock;

	switch (ROISize)
	{
	case 5:
		gpu_CalculatePSFWidth<5> << < BlockNum, BlockDim, 0, cstream >> >(d_ImageROI, d_WLEPara, FluoNum);

		break;
	case 7:
		gpu_CalculatePSFWidth<7> << < BlockNum, BlockDim, 0, cstream >> >(d_ImageROI, d_WLEPara, FluoNum);

		break;
	case 9:
		gpu_CalculatePSFWidth<9> << < BlockNum, BlockDim, 0, cstream >> >(d_ImageROI, d_WLEPara, FluoNum);

		break;
	case 11:
		gpu_CalculatePSFWidth<11> << < BlockNum, BlockDim, 0, cstream >> >(d_ImageROI, d_WLEPara, FluoNum);

		break;
	case 13:
		gpu_CalculatePSFWidth<13> << < BlockNum, BlockDim, 0, cstream >> >(d_ImageROI, d_WLEPara, FluoNum);

		break;
	case 15:
		gpu_CalculatePSFWidth<15> << < BlockNum, BlockDim, 0, cstream >> >(d_ImageROI, d_WLEPara, FluoNum);

		break;
	case 17:
		gpu_CalculatePSFWidth<17> << < BlockNum, BlockDim, 0, cstream >> >(d_ImageROI, d_WLEPara, FluoNum);

		break;
	case 19:
		gpu_CalculatePSFWidth<19> << < BlockNum, BlockDim, 0, cstream >> >(d_ImageROI, d_WLEPara, FluoNum);

		break;
	}
}


__device__ void SigmaWidthCalc_TwoSides(float *SigmaL, float *SigmaR, float *MeanData_Intp, int DataLen, int MaxPos)
{
	float MaxValue = MeanData_Intp[MaxPos];

	int LLen = MaxPos;
	int RLen = DataLen - MaxPos - 1;

	int CmpLen = min(LLen, RLen);

	// left side
	float SDat_L = 0;

	for (int i = MaxPos - CmpLen; i <= MaxPos; i++)
	{
		SDat_L = SDat_L + MeanData_Intp[i];
	}

	// right side
	float SDat_R = 0;

	for (int i = MaxPos; i <= MaxPos + CmpLen; i++)
	{
		SDat_R = SDat_R + MeanData_Intp[i];
	}

	*SigmaL = SDat_L * InterpolationGap / MaxValue / 1.2533f; // sqrt(pi/2)
	*SigmaR = SDat_R * InterpolationGap / MaxValue / 1.2533f; // sqrt(pi/2)
}

__device__ int FindMaxPos(float *iData, int DataLen, int LSel, int RSel)
{
	LSel = max(LSel, 0);
	RSel = min(RSel, DataLen - 1);

	int MaxIdx = DataLen / 2;

	for (int cnt = LSel; cnt <= RSel; cnt++)
	{
		if (iData[cnt] > iData[MaxIdx])
		{
			MaxIdx = cnt;
		}
	}

	return MaxIdx;
}



__global__ void gpu_CalculateNearestNeighborDistance(unsigned short * d_ImageROI, int ROISize, float *d_WLEPara, int FluoNum)
{
	int gid = threadIdx.x + blockDim.x*blockIdx.x;

	float(*pWLEPara)[WLE_ParaNumber] = (float(*)[WLE_ParaNumber])d_WLEPara;

	const int ROIWholeSize = ROISize*(ROISize + 1);

	const int AddrOffset = ROIWholeSize*gid + ROISize*ROISize;

	if (gid < FluoNum)
	{
		int CurX = d_ImageROI[AddrOffset + 0];
		int CurY = d_ImageROI[AddrOffset + 1];
		int CurF = d_ImageROI[AddrOffset + 2] + d_ImageROI[AddrOffset + 3] * 65536;

		int NextX = 0;
		int NextY = 0;
		int NextF = 0;


		float NearestNeighborDistance = 100000000.0f;

		for (int cnt = 0; cnt < FluoNum; cnt++)
		{
			int NextAddrOffset = ROIWholeSize*cnt + ROISize*ROISize;

			NextX = d_ImageROI[NextAddrOffset + 0];
			NextY = d_ImageROI[NextAddrOffset + 1];
			NextF = d_ImageROI[NextAddrOffset + 2] + d_ImageROI[NextAddrOffset + 3] * 65536;

			int DistanceBias = (NextF != CurF)*100000000.0f + (cnt == gid)*100000000.0f;

			float curDistance = (CurX - NextX)* (CurX - NextX) + (CurY - NextY)* (CurY - NextY) + DistanceBias;

			NearestNeighborDistance = min(NearestNeighborDistance, curDistance);
		}

		NearestNeighborDistance = sqrtf(NearestNeighborDistance);

		pWLEPara[gid][WLE_Para_NearDistance] = NearestNeighborDistance;

	}
}

void CalculateNearestNeighborDistance(unsigned short * d_ImageROI, int ROISize, float *d_WLEPara, int FluoNum, cudaStream_t cstream)
{

	int BlockDim = ThreadsPerBlock;
	int BlockNum = (FluoNum + ThreadsPerBlock - 1) / ThreadsPerBlock;

	gpu_CalculateNearestNeighborDistance << < BlockNum, BlockDim, 0, cstream >> > (d_ImageROI, ROISize, d_WLEPara, FluoNum);

}


// add imaging type
__global__ void gpu_MoleculeTypeClasify(unsigned short * d_ImageROI, int LocType, int MultiEmitterFitEn, int ROISize, float *d_WLEPara, int FluoNum)
{
	int gid = threadIdx.x + blockDim.x*blockIdx.x;

	float(*pWLEPara)[WLE_ParaNumber] = (float(*)[WLE_ParaNumber])d_WLEPara;

	const int ROIWholeSize = ROISize*(ROISize + 1);

	const int AddrOffset = ROIWholeSize*gid + ROISize*ROISize;

	int ROISize_Half = ROISize / 2;

	float MaxMeanSigmaTh = 0;
	

//	MaxMeanSigmaTh = ROISize / 2.0f / 2.35f;


	if (LocType == LocType_GS2D)
	{
		MaxMeanSigmaTh = ROISize / 2.1f / 2.35f;
	}
	else
	{
		MaxMeanSigmaTh = ROISize / 1.9f / 2.35f;
	}


	if (gid < FluoNum)
	{
		int ROIType_Detection = d_ImageROI[AddrOffset + 4];


		float NearestNeighborDistance = pWLEPara[gid][WLE_Para_NearDistance];
		float SigmaL = pWLEPara[gid][WLE_Para_SigmaL];
		float SigmaR = pWLEPara[gid][WLE_Para_SigmaR];
		float SigmaU = pWLEPara[gid][WLE_Para_SigmaU];
		float SigmaD = pWLEPara[gid][WLE_Para_SigmaD];

	
		float diff1 = abs(SigmaL - SigmaR) / min(SigmaL, SigmaR);
		float diff2 = abs(SigmaU - SigmaD) / min(SigmaU, SigmaD);

		int IsContamination = (diff1 >= MaxTolerablePSFWidthDiff) || (diff2 >= MaxTolerablePSFWidthDiff);

		float MeanSigma = (SigmaL + SigmaR + SigmaU + SigmaD) / 4;
		//		float MeanSigma = (min(SigmaL, SigmaR) + min(SigmaU, SigmaD)) / 2;


		int MoleculeType = MoleculeType_MLEFit; // single molecule, 1: multiXYPos

		if ((ROIType_Detection == ROIType_Fit_Single) && (IsContamination == 0) && (NearestNeighborDistance >= ROISize) && (MeanSigma <= MaxMeanSigmaTh))
		{
			MoleculeType = MoleculeType_MLEFit; // well isolated molecule
		}
		else if ((ROIType_Detection == ROIType_Fit_Multi) || (MeanSigma > MaxMeanSigmaTh) || (NearestNeighborDistance < (ROISize / 2.0f + 1.0f)))
		{
			MoleculeType = MoleculeType_MultiFit;
		}
		else
		{
			MoleculeType = MoleculeType_WLEFit;
		}

		pWLEPara[gid][WLE_Para_FluoType] = MoleculeType;

	}
}


void MoleculeTypeClasify(unsigned short * d_ImageROI, int LocType, int MultiEmitterFitEn, int ROISize, float *d_WLEPara, int FluoNum, cudaStream_t cstream)
{

	int BlockDim = ThreadsPerBlock;
	int BlockNum = (FluoNum + ThreadsPerBlock - 1) / ThreadsPerBlock;

	gpu_MoleculeTypeClasify << < BlockNum, BlockDim, 0, cstream >> > (d_ImageROI, LocType, MultiEmitterFitEn, ROISize, d_WLEPara, FluoNum);
}


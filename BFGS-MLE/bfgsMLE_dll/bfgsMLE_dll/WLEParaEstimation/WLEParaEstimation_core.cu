#include "WLEParaEstimation.h"


__global__ void gpu_CalculateNearestNeighborDistance(unsigned short * d_ROIMem, int ROISize, float *d_WLEPara, int FluoNum);

__global__ void gpu_MoleculeTypeClasify(int ROISize, float *d_WLEPara, int FluoNum);


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


void CalculatePSFWidth(unsigned short * d_ROIMem, float *d_WLEPara, int FluoNum, int ROISize, cudaStream_t cstream)
{

	int BlockDim = ThreadsPerBlock;
	int BlockNum = (FluoNum + ThreadsPerBlock - 1) / ThreadsPerBlock;

	switch (ROISize)
	{
	case 5:
		gpu_CalculatePSFWidth<5> << < BlockNum, BlockDim, 0, cstream >> >(d_ROIMem, d_WLEPara, FluoNum);

		break;
	case 7:
		gpu_CalculatePSFWidth<7> << < BlockNum, BlockDim, 0, cstream >> >(d_ROIMem, d_WLEPara, FluoNum);

		break;
	case 9:
		gpu_CalculatePSFWidth<9> << < BlockNum, BlockDim, 0, cstream >> >(d_ROIMem, d_WLEPara, FluoNum);

		break;
	case 11:
		gpu_CalculatePSFWidth<11> << < BlockNum, BlockDim, 0, cstream >> >(d_ROIMem, d_WLEPara, FluoNum);

		break;
	case 13:
		gpu_CalculatePSFWidth<13> << < BlockNum, BlockDim, 0, cstream >> >(d_ROIMem, d_WLEPara, FluoNum);

		break;
	case 15:
		gpu_CalculatePSFWidth<15> << < BlockNum, BlockDim, 0, cstream >> >(d_ROIMem, d_WLEPara, FluoNum);

		break;
	case 17:
		gpu_CalculatePSFWidth<17> << < BlockNum, BlockDim, 0, cstream >> >(d_ROIMem, d_WLEPara, FluoNum);

		break;
	case 19:
		gpu_CalculatePSFWidth<19> << < BlockNum, BlockDim, 0, cstream >> >(d_ROIMem, d_WLEPara, FluoNum);

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

void CalculateNearestNeighborDistance(unsigned short * d_ROIMem, int ROISize, float *d_WLEPara, int FluoNum, cudaStream_t cstream)
{

	int BlockDim = ThreadsPerBlock;
	int BlockNum = (FluoNum + ThreadsPerBlock - 1) / ThreadsPerBlock;
	
	gpu_CalculateNearestNeighborDistance << < BlockNum, BlockDim, 0, cstream >> > (d_ROIMem, ROISize, d_WLEPara, FluoNum);

}


__global__ void gpu_CalculateNearestNeighborDistance(unsigned short * d_ROIMem, int ROISize, float *d_WLEPara, int FluoNum)
{
	int gid = threadIdx.x + blockDim.x*blockIdx.x;

	float(*pWLEPara)[WLE_ParaNumber] = (float(*)[WLE_ParaNumber])d_WLEPara;

	const int ROIWholeSize = ROISize*(ROISize + 1);

	const int AddrOffset = ROIWholeSize*gid + ROISize*ROISize;

	if (gid < FluoNum)
	{
		int CurX = d_ROIMem[AddrOffset + 0];
		int CurY = d_ROIMem[AddrOffset + 1];
		int CurF = d_ROIMem[AddrOffset + 2] + d_ROIMem[AddrOffset + 3] * 65536;

		int NextX = 0;
		int NextY = 0;
		int NextF = 0;


		float NearestNeighborDistance = 100000000.0f;

		for (int cnt = 0; (cnt < FluoNum) && (NextF <= CurF); cnt++)
		{
			int NextAddrOffset = ROIWholeSize*cnt + ROISize*ROISize;

			NextX = d_ROIMem[NextAddrOffset + 0];
			NextY = d_ROIMem[NextAddrOffset + 1];
			NextF = d_ROIMem[NextAddrOffset + 2] + d_ROIMem[NextAddrOffset + 3] * 65536;

			int DistanceBias = (NextF != CurF)*100000000.0f + (cnt == gid)*100000000.0f;

			float curDistance = (CurX - NextX)* (CurX - NextX) + (CurY - NextY)* (CurY - NextY) + DistanceBias;

			NearestNeighborDistance = min(NearestNeighborDistance, curDistance);
		}

		NearestNeighborDistance = sqrtf(NearestNeighborDistance);

		pWLEPara[gid][WLE_Para_NearDistance] = NearestNeighborDistance;

	}
}

void MoleculeTypeClasify(int ROISize, float *d_WLEPara, int FluoNum, cudaStream_t cstream)
{

	int BlockDim = ThreadsPerBlock;
	int BlockNum = (FluoNum + ThreadsPerBlock - 1) / ThreadsPerBlock;

	gpu_MoleculeTypeClasify << < BlockNum, BlockDim, 0, cstream >> > (ROISize, d_WLEPara, FluoNum);
}


__global__ void gpu_MoleculeTypeClasify(int ROISize,float *d_WLEPara, int FluoNum)
{
	int gid = threadIdx.x + blockDim.x*blockIdx.x;

	float(*pWLEPara)[WLE_ParaNumber] = (float(*)[WLE_ParaNumber])d_WLEPara;


	float MultiEmitterFit_Distance_Th = (ROISize - 1) / 2.0f;
	float MaxSigmaTh = 7 / 1.9f / 2.35f;

	if (gid < FluoNum)
	{
		float NearestNeighborDistance = pWLEPara[gid][WLE_Para_NearDistance];
		float SigmaL = pWLEPara[gid][WLE_Para_SigmaL];
		float SigmaR = pWLEPara[gid][WLE_Para_SigmaR];
		float SigmaU = pWLEPara[gid][WLE_Para_SigmaU];
		float SigmaD = pWLEPara[gid][WLE_Para_SigmaD];


		float diff1 = abs(SigmaL - SigmaR) / min(SigmaL, SigmaR);
		float diff2 = abs(SigmaU - SigmaD) / min(SigmaU, SigmaD);

		int IsContamination = (diff1 >= MaxTolerablePSFWidthDiff) || (diff2 >= MaxTolerablePSFWidthDiff);


		SigmaL = max(SigmaL, SigmaR);
		SigmaU = max(SigmaU, SigmaD);


		int MoleculeType = MoleculeType_MLEFit; // single molecule, 1: multi


		if ((IsContamination == 0) && (NearestNeighborDistance >= ROISize))
		{
			MoleculeType = MoleculeType_MLEFit; // well isolated molecule

		}
		else if ((NearestNeighborDistance <= MultiEmitterFit_Distance_Th))
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

#include "LDROIExtraction.h"

// merge ROI at the same position in adjecent frames in a batch frame, by average ROI
__global__ void gpuFindConsecutiveROI(unsigned short * d_ROIMem, int ROISize, int *d_ForwardLinkID, int *d_BackwardLinkID, int FluoNum);
__global__ void gpuFindConsecutiveNum(int *d_ForwardLinkID, int *d_BackwardLinkID, int *d_ConsecutiveNum, int FluoNum);
__global__ void gpuMergeROIData(unsigned short * d_ROIMem, int ROISize, int *d_ForwardLinkID, int *d_BackwardLinkID, int *d_ConsecutiveNum, int FluoNum);
__global__ void gpuMergeWLEPara(float *d_WLEPara, int *d_ForwardLinkID, int *d_BackwardLinkID, int *d_ConsecutiveNum, int FluoNum);



void MergeConsecutiveROI(unsigned short * d_ROIMem, float *d_WLEPara, int ROISize, int *d_ForwardLinkID, int *d_BackwardLinkID, int *d_ConsecutiveNum, int FluoNum, cudaStream_t cstream)
{

	int BlockDim = ThreadsPerBlock;
	int BlockNum = (FluoNum + ThreadsPerBlock - 1) / ThreadsPerBlock;

	cudaMemsetAsync(d_ForwardLinkID, 0, FluoNum * sizeof(int), cstream);
	cudaMemsetAsync(d_BackwardLinkID, 0, FluoNum * sizeof(int), cstream);
	cudaMemsetAsync(d_ConsecutiveNum, 0, FluoNum * sizeof(int), cstream);

	gpuFindConsecutiveROI << < BlockNum, BlockDim, 0, cstream >> > (d_ROIMem, ROISize, d_ForwardLinkID, d_BackwardLinkID, FluoNum);
	gpuFindConsecutiveNum << < BlockNum, BlockDim, 0, cstream >> > (d_ForwardLinkID, d_BackwardLinkID, d_ConsecutiveNum, FluoNum);
	gpuMergeROIData << < BlockNum, BlockDim, 0, cstream >> > (d_ROIMem, ROISize, d_ForwardLinkID, d_BackwardLinkID, d_ConsecutiveNum, FluoNum);
	gpuMergeWLEPara << < BlockNum, BlockDim, 0, cstream >> > (d_WLEPara, d_ForwardLinkID, d_BackwardLinkID, d_ConsecutiveNum, FluoNum);

}


// merge ROI at the same position in adjecent frames in a batch frame, by average ROI
__global__ void gpuFindConsecutiveROI(unsigned short * d_ROIMem, int ROISize, int *d_ForwardLinkID, int *d_BackwardLinkID, int FluoNum)
{
	//    int tid = threadIdx.x;
	int gid = threadIdx.x + blockDim.x*blockIdx.x;

	int gid_0 = blockDim.x*blockIdx.x;

	const int ROIDataLen = ROISize*(ROISize + 1);

	const int AddrOffset = ROIDataLen*gid + ROISize*ROISize;


	if (gid < FluoNum)
	{

		int CurX = d_ROIMem[AddrOffset + 0];
		int CurY = d_ROIMem[AddrOffset + 1];
		int CurF = d_ROIMem[AddrOffset + 2] + d_ROIMem[AddrOffset + 3] * 65536;

		int NextX = 0;
		int NextY = 0;
		int NextF = 0;

		int PairID = -2;
		int SimilarFluoNum = 0;

		// molecule's position in next frame is larger
		for (int cnt = gid_0; (cnt < FluoNum) && (NextF <= CurF + 1); cnt++)
		{
			int NextAddrOffset = ROIDataLen*cnt + ROISize*ROISize;

			NextX = d_ROIMem[NextAddrOffset + 0];
			NextY = d_ROIMem[NextAddrOffset + 1];
			NextF = d_ROIMem[NextAddrOffset + 2] + d_ROIMem[NextAddrOffset + 3] * 65536;

			// find the closest in next frame
			if ((CurX == NextX) && (CurY == NextY) && (NextF - CurF == 1))
			{
				PairID = cnt;
				SimilarFluoNum++;

//				printf("consec: %d %d - %d %d %d, %d %d %d\n", gid, PairID, CurX, CurY, CurF, NextX, NextY, NextF);
			}
		}
		// there is consecutive molecule
		if (SimilarFluoNum > 0)
		{
			d_BackwardLinkID[gid] = PairID;
			d_ForwardLinkID[PairID] = 1; // has foward
		}
		else
		{
			// invalid molecules
			d_BackwardLinkID[gid] = -2;
//			d_ForwardLinkID[gid] = -2;
		}
	}
}

__global__ void gpuFindConsecutiveNum(int *d_ForwardLinkID, int *d_BackwardLinkID, int *d_ConsecutiveNum, int FluoNum)
{
	//    int tid = threadIdx.x;
	int gid = threadIdx.x + blockDim.x*blockIdx.x;

	if (gid < FluoNum)
	{
		int ConsecutiveROINum = 0;

		// the first fluo is valid fluo and don't have a forward
		if (d_ForwardLinkID[gid] == 0)
		{
			int CurId = gid;
			int NextID = d_BackwardLinkID[CurId];

			while (NextID > 0)
			{
				NextID = d_BackwardLinkID[NextID];

				ConsecutiveROINum++;
			}
		}

		d_ConsecutiveNum[gid] = ConsecutiveROINum;
	}
}


__global__ void gpuMergeROIData(unsigned short * d_ROIMem, int ROISize, int *d_ForwardLinkID, int *d_BackwardLinkID, int *d_ConsecutiveNum, int FluoNum)
{
	//    int tid = threadIdx.x;
	int gid = threadIdx.x + blockDim.x*blockIdx.x;


	const int ROIDataLen = ROISize*(ROISize + 1);
	const int ROIDataSize = ROISize*ROISize;


	if (gid < FluoNum)
	{
		int ConsecutiveROINum = d_ConsecutiveNum[gid];

		float Weights = 1.0f / (ConsecutiveROINum + 1);

		int CurId = 0;
		int NextValid = 0;


		// the first fluo is valid fluo and don't have a forward
		if (d_ForwardLinkID[gid] == 0)
		{
			if (ConsecutiveROINum > 0)
			{
				for (int i = 0; i < ROIDataSize; i++)
				{
					// loop of roi data, get merged data
					float SumDdata = 0;

					CurId = gid;
					NextValid = 1;

					while (NextValid)
					{
						int AddrOffset = ROIDataLen*CurId;

						// ROI data
						SumDdata += Weights*d_ROIMem[AddrOffset + i];


						CurId = d_BackwardLinkID[CurId];
						NextValid = CurId > 0;
					}

					// set merged data
					CurId = gid;
					NextValid = 1;

					while (NextValid)
					{
						int AddrOffset = ROIDataLen*CurId;

						// ROI data
						d_ROIMem[AddrOffset + i] = SumDdata + 0.5f;// round off


						CurId = d_BackwardLinkID[CurId];
						NextValid = CurId > 0;
					}
				}
			}
		}
	}
}



__global__ void gpuMergeWLEPara(float *d_WLEPara, int *d_ForwardLinkID, int *d_BackwardLinkID, int *d_ConsecutiveNum,  int FluoNum)
{
	//    int tid = threadIdx.x;
	int gid = threadIdx.x + blockDim.x*blockIdx.x;

	float(*pWLEPara)[WLE_ParaNumber] = (float(*)[WLE_ParaNumber])d_WLEPara;


	if (gid < FluoNum)
	{
		int ConsecutiveROINum = d_ConsecutiveNum[gid];
		float Weights = 1.0f / (ConsecutiveROINum + 1);

		int CurId = 0;
		int NextValid = 0;


		// the first fluo is valid fluo and don't have a forward
		if (d_ForwardLinkID[gid] == 0)
		{
			if (ConsecutiveROINum > 0)
			{
				float NearestNeighborDistance = 1000;
				float SigmaL = 0;
				float SigmaR = 0;
				float SigmaU = 0;
				float SigmaD = 0;
				float MoleculeType = 0;

				CurId = gid;
				NextValid = 1;

				while (NextValid)
				{
					// WLE para
					NearestNeighborDistance = min(NearestNeighborDistance, pWLEPara[CurId][WLE_Para_NearDistance]);

					SigmaL += Weights * pWLEPara[CurId][WLE_Para_SigmaL];
					SigmaR += Weights * pWLEPara[CurId][WLE_Para_SigmaR];
					SigmaU += Weights * pWLEPara[CurId][WLE_Para_SigmaU];
					SigmaD += Weights * pWLEPara[CurId][WLE_Para_SigmaD];

					MoleculeType = max(MoleculeType, pWLEPara[CurId][WLE_Para_FluoType]);

					CurId = d_BackwardLinkID[CurId];
					NextValid = CurId > 0;
				}

				// set merged data
				CurId = gid;
				NextValid = 1;

				while (NextValid)
				{
					// WLE para
					pWLEPara[CurId][WLE_Para_NearDistance] = NearestNeighborDistance;
					pWLEPara[CurId][WLE_Para_SigmaL] = SigmaL;
					pWLEPara[CurId][WLE_Para_SigmaR] = SigmaR;
					pWLEPara[CurId][WLE_Para_SigmaU] = SigmaU;
					pWLEPara[CurId][WLE_Para_SigmaD] = SigmaD;
					pWLEPara[CurId][WLE_Para_FluoType] = MoleculeType;

					CurId = d_BackwardLinkID[CurId];
					NextValid = CurId > 0;
				}
			}
		}
	}
}


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

#include "ConsecutiveFilter.h"

__global__ void gpuConsecutiveFilterPair(float *d_LocArry, int *d_ForwardLinkID, int *d_BackwardLinkID, float Distance_th_pixel, int FluoNum);
__global__ void gpuFindMaxSNRID(float *d_LocArry, int *d_ForwardLinkID, int *d_BackwardLinkID, int *d_MaxSNRID, int FluoNum);
__global__ void gpuRemoveConsecutiveFluo_KeepMaxSNR(float *d_LocArry, int *d_ForwardLinkID, int *d_BackwardLinkID, int *d_MaxSNRID, int FluoNum);
__global__ void gpuRemoveConsecutiveFluo_KeepFirst(float *d_LocArry, int *d_ForwardLinkID, int *d_BackwardLinkID, int FluoNum);


// d_iLocArry come from localization data 

void NyqConsecutiveFilter_TypeDef::FilterConsecutiveFluo(float * d_LocArry, int FluoNum, float Distance_th_pixel, cudaStream_t cstream)
{

	int BlockDim = ThreadsPerBlock;

	int BlockNum = (FluoNum + ThreadsPerBlock - 1) / ThreadsPerBlock;


	cudaMemsetAsync(d_ForwardLinkID, 0, FluoNum * sizeof(int), cstream);
	cudaMemsetAsync(d_BackwardLinkID, 0, FluoNum * sizeof(int), cstream);
	cudaMemsetAsync(d_MaxSNRID, 0, FluoNum * sizeof(int), cstream);


	gpuConsecutiveFilterPair << <BlockNum, BlockDim, 0, cstream >> >(d_LocArry, d_ForwardLinkID, d_BackwardLinkID, Distance_th_pixel, FluoNum);

	//	gpuRemoveConsecutiveFluo_KeepFirst << <BlockNum, BlockDim, 0, cstream >> >(d_LocArry, d_ForwardLinkID, d_BackwardLinkID, FluoNum);

	gpuFindMaxSNRID << <BlockNum, BlockDim, 0, cstream >> >(d_LocArry, d_ForwardLinkID, d_BackwardLinkID, d_MaxSNRID, FluoNum);
	gpuRemoveConsecutiveFluo_KeepMaxSNR << <BlockNum, BlockDim, 0, cstream >> >(d_LocArry, d_ForwardLinkID, d_BackwardLinkID, d_MaxSNRID, FluoNum);


	cudaStreamSynchronize(cstream);

}



void NyqConsecutiveFilter_TypeDef::Init(int TotalFluoNum)
{
	//	cudaMallocHost((void **)&h_LocArry, MaxFluoNum * OutParaNumGS2D*sizeof(float));
	//	cudaMalloc((void **)&d_LocArry, TotalFluoNum * OutParaNumGS2D*sizeof(float));


	cudaMalloc((void **)&d_ForwardLinkID, TotalFluoNum  * sizeof(int));
	cudaMalloc((void **)&d_BackwardLinkID, TotalFluoNum  * sizeof(int));
	cudaMalloc((void **)&d_MaxSNRID, TotalFluoNum  * sizeof(int));

}


void NyqConsecutiveFilter_TypeDef::DeInit()
{
	//	cudaFreeHost(h_LocArry);
	//	cudaFree(d_LocArry);


	cudaFree(d_ForwardLinkID);
	cudaFree(d_BackwardLinkID);
	cudaFree(d_MaxSNRID);

}

////////////////////////////////////////
// cuda low level implements //


__global__ void gpuConsecutiveFilterPair(float *d_LocArry, int *d_ForwardLinkID, int *d_BackwardLinkID, float Distance_th_pixel, int FluoNum)
{
	//    int tid = threadIdx.x;

	int gid = blockDim.x*blockIdx.x + threadIdx.x;

	// all thread in a warp read from the same memory to improve efficiency
	int gid_0 = gid / 32 * 32;


	float(*pLocArry)[OutParaNumGS2D]; //
	pLocArry = (float(*)[OutParaNumGS2D])d_LocArry;


	float XPos1;
	float XPos2;

	float YPos1;
	float YPos2;

	int Frame1 = 0;
	int Frame2 = 0;

	int cnt = 0;

	float DistanceTh2; // square of distance

	float CurDistance2;

	int SimilarFluoNum = 0;

	int PairID = 0;

	if (gid < FluoNum)
	{

		XPos1 = pLocArry[gid][Pos_XPos];
		YPos1 = pLocArry[gid][Pos_YPos];
		Frame1 = pLocArry[gid][Pos_Frme];


		// square of distance, to avoid sqrt
		DistanceTh2 = Distance_th_pixel*Distance_th_pixel;

		if ((XPos1 > 5) && (YPos1 > 5))
		{
			// molecule's position in next frame is larger
			for (cnt = gid_0; (Frame2 <= Frame1 + 1) && (cnt < FluoNum); cnt++)
			{
				XPos2 = pLocArry[cnt][Pos_XPos];
				YPos2 = pLocArry[cnt][Pos_YPos];
				Frame2 = pLocArry[cnt][Pos_Frme];


				CurDistance2 = (XPos1 - XPos2)*(XPos1 - XPos2) + (YPos1 - YPos2)*(YPos1 - YPos2);

				// find the closest in next frame
				if ((CurDistance2 < DistanceTh2) && (Frame2 - Frame1 == 1))
				{
					PairID = cnt;
					DistanceTh2 = CurDistance2; // update distance th for find the closest

					SimilarFluoNum++;
				}
			}

			if (SimilarFluoNum > 0)
			{
				d_BackwardLinkID[gid] = PairID;
				d_ForwardLinkID[PairID] = 1; // has foward
			}
			else
			{
				// if don't have a backword consecutive fluo
				d_BackwardLinkID[gid] = -1;
			}
		}
		else
		{
			// invalid molecules
			d_ForwardLinkID[gid] = -2;
			d_BackwardLinkID[gid] = -2;
		}
	}
}

__global__ void gpuFindMaxSNRID(float *d_LocArry, int *d_ForwardLinkID, int *d_BackwardLinkID, int *d_MaxSNRID, int FluoNum)
{
	//    int tid = threadIdx.x;
	int gid = blockDim.x*blockIdx.x + threadIdx.x;


	float(*pLocArry)[OutParaNumGS2D]; // for 9 parameter array
	pLocArry = (float(*)[OutParaNumGS2D])d_LocArry;


	int cnt;

	int NextID;

	int MaxSNRID = 0;
	float CurSNR = 0;
	float NexSNR = 0;

	if (gid < FluoNum)
	{
		// the first fluo is valid fluo and don't have a forward, may have a backword
		if (d_ForwardLinkID[gid] == 0)
		{

			MaxSNRID = gid;
			CurSNR = pLocArry[gid][Pos_PSNR];

			NextID = d_BackwardLinkID[gid];

			while (NextID > 0)
			{
				NexSNR = pLocArry[NextID][Pos_PSNR];

				if (NexSNR > CurSNR)
				{
					CurSNR = NexSNR;
					MaxSNRID = NextID;
				}

				NextID = d_BackwardLinkID[NextID];
			}

			// only write it for the leading fluo in a sequence
			d_MaxSNRID[gid] = MaxSNRID;
		}
	}
}

// keep the max SNR but not the first one in a sequence
__global__ void gpuRemoveConsecutiveFluo_KeepMaxSNR(float *d_LocArry, int *d_ForwardLinkID, int *d_BackwardLinkID, int *d_MaxSNRID, int FluoNum)
{
	//    int tid = threadIdx.x;
	int gid = blockDim.x*blockIdx.x + threadIdx.x;


	float(*pLocArry)[OutParaNumGS2D]; // for parameter array
	pLocArry = (float(*)[OutParaNumGS2D])d_LocArry;


	int cnt;

	int NextID;
	int MaxID;

	if (gid < FluoNum)
	{
		// the first fluo is valid fluo and don't have a forward, may have a backword
		if (d_ForwardLinkID[gid] == 0)
		{

			NextID = gid;
			MaxID = d_MaxSNRID[gid];

			// NextID could be 0 
			while (NextID >= 0)
			{
				// remove fluo in a sequency if it not the max snr id
				if (NextID != MaxID)
				{
#pragma unroll
					for (cnt = 0; cnt < OutParaNumGS2D; cnt++)
					{
						pLocArry[NextID][cnt] = 0;
					}
				}

				NextID = d_BackwardLinkID[NextID];
			}
		}
	}
}

__global__ void gpuRemoveConsecutiveFluo_KeepFirst(float *d_LocArry, int *d_ForwardLinkID, int *d_BackwardLinkID, int FluoNum)
{
	//    int tid = threadIdx.x;
	int gid = blockDim.x*blockIdx.x + threadIdx.x;


	float(*pLocArry)[OutParaNumGS2D]; // for 9 parameter array
	pLocArry = (float(*)[OutParaNumGS2D])d_LocArry;


	int cnt;

	int NextID;

	if (gid < FluoNum)
	{
		// the first fluo is valid fluo and don't have a forward, may have a backword
		if (d_ForwardLinkID[gid] == 0)
		{
			NextID = d_BackwardLinkID[gid];

			// remove fluo in a sequency if it not the first id
			while (NextID > 0)
			{
#pragma unroll
				for (cnt = 0; cnt < OutParaNumGS2D; cnt++)
				{
					pLocArry[NextID][cnt] = 0;
				}

				NextID = d_BackwardLinkID[NextID];
			}
		}
	}
}

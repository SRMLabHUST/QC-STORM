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

#include "bfgs_top.h"
#include "bfgs_OntimeCalc.h"


// code wrapper
void LDLocData_TypeDef::CopyDataToGPU(float * h_iLocArry, int FluoNum, cudaStream_t cstream)
{
	cudaMemcpyAsync(d_LocArry, h_iLocArry, FluoNum*OutParaNumGS2D*sizeof(float), cudaMemcpyHostToDevice, cstream);
	cudaStreamSynchronize(cstream);

}

void LDLocData_TypeDef::OntimeCalc(LocalizationPara & LocPara, int FluoNum, cudaStream_t cstream)
{
	// ontime calculate
	memset(h_OntimeRatio, 0, MaxOnTimeConsecutiveNum*sizeof(float));

	// calculate on time after localizaiton
	int StartFrame = GetFirstFrame(h_LocArry, FluoNum);
	int EndFrame = GetLastFrame(h_LocArry, FluoNum);


	int TotalFrame = EndFrame - StartFrame + 1;
	int FilterStartFrame = StartFrame + TotalFrame * 3 / 10;
	int FilterEndFrame = StartFrame + TotalFrame * 7 / 10;

	cudaMemsetAsync(d_ForwardLinkID, 0, FluoNum * sizeof(int), cstream);
	cudaMemsetAsync(d_BackwardLinkID, 0, FluoNum * sizeof(int), cstream);
	cudaMemsetAsync(d_ConsecutiveNum, 0, FluoNum * sizeof(int), cstream);


	//
	cudaMemsetAsync(d_OntimeDistrib, 0, MaxOnTimeConsecutiveNum*sizeof(int), cstream);
	cudaMemsetAsync(d_ValidFluoNum, 0, sizeof(int), cstream);

	int BlockDim = ThreadsPerBlock;
	int BlockNum = ((FluoNum + ThreadsPerBlock - 1) / ThreadsPerBlock);


	// most time comsuming, because compare each molecule to others
	gpuOntimeCalcFluoPair << <BlockNum, BlockDim, 0, cstream >> >(d_LocArry, d_ForwardLinkID, d_BackwardLinkID, LocPara.PixelSize, LocPara.QE, FluoNum);

	gpuOntimeCalc << <BlockNum, BlockDim, 0, cstream >> >(d_LocArry, d_ForwardLinkID, d_BackwardLinkID, d_ConsecutiveNum, LocPara.PixelSize, LocPara.QE, FluoNum);

	gpuGetOntimeDistribution << <BlockNum, BlockDim, 0, cstream >> >(d_LocArry, d_ConsecutiveNum, d_OntimeDistrib, d_ValidFluoNum, FilterStartFrame, FilterEndFrame, FluoNum);


	// calculate on time (emitting in consecutive frames) distribution
	cudaMemcpyAsync(h_OntimeDistrib, d_OntimeDistrib, MaxOnTimeConsecutiveNum*sizeof(int), cudaMemcpyDeviceToHost, cstream);
	cudaMemcpyAsync(h_ValidFluoNum, d_ValidFluoNum, sizeof(int), cudaMemcpyDeviceToHost, cstream);

	cudaStreamSynchronize(cstream);

	// ontime distribution ratio
	int cnt;
	if (*h_ValidFluoNum <= 0)*h_ValidFluoNum = 1;
	for (cnt = 0; cnt < MaxOnTimeConsecutiveNum; cnt++)
	{
		h_OntimeRatio[cnt] = h_OntimeDistrib[cnt] * 1.0f / (*(h_ValidFluoNum));
	}

}


//////////////////////////////////////


__global__ void gpuOntimeCalcFluoPair(float *d_LocArry, int *d_ForwardLinkID, int *d_BackwardLinkID, float PixelSize, float QE, int FluoNum)
{
	//    int tid = threadIdx.x;
	int gid = threadIdx.x + blockDim.x*blockIdx.x;
	int gid_0 = blockDim.x*blockIdx.x;

	float(*pLocArry)[OutParaNumGS2D]; // for 9 parameter array

	pLocArry = (float(*)[OutParaNumGS2D])d_LocArry;


	int PairID;
	float XPos1;
	float XPos2;

	float YPos1;
	float YPos2;


	int Frame1 = 0;
	int Frame2 = 0;

	int cnt = 0;

	float DistanceTh2 = 0; // square of distance

	float CurDistance2 = 0;

	int SimilarFluoNum = 0;


	if (gid < FluoNum)
	{
		XPos1 = pLocArry[gid][Pos_XPos];
		YPos1 = pLocArry[gid][Pos_YPos];
		Frame1 = pLocArry[gid][Pos_Frme];

		// can be revised to be loc precision
//		DistanceTh2 = 30.0f / PixelSize; // 30nm
		DistanceTh2 = 0.5f; // 0.5 pixel
		DistanceTh2 = DistanceTh2*DistanceTh2;

		if ((XPos1 > 5) && (YPos1 > 5))
		{
			// molecule's position in next frame is larger
			for (cnt = gid_0; (cnt < FluoNum) && (Frame2 <= Frame1 + 1); cnt++)
			{
				XPos2 = pLocArry[cnt][Pos_XPos];
				YPos2 = pLocArry[cnt][Pos_YPos];
				Frame2 = pLocArry[cnt][Pos_Frme];


				CurDistance2 = (XPos1 - XPos2)*(XPos1 - XPos2) + (YPos1 - YPos2)*(YPos1 - YPos2);

				// find the closest in next frame
				if ((CurDistance2 < DistanceTh2) && (Frame2 - Frame1 == 1))
				{
					PairID = cnt;
					SimilarFluoNum++;
					DistanceTh2 = CurDistance2; // update distance th for find the closest
				}

			}

			if (SimilarFluoNum > 0)
			{
				// if no molecule is paired, it's 0
				d_BackwardLinkID[gid] = PairID;
				d_ForwardLinkID[PairID] = 1; // has foward
			}

		}
		else
		{
			// invalid molecules
			d_ForwardLinkID[gid] = -2;
		}
	}
}

__global__ void gpuOntimeCalc(float *d_LocArry, int *d_ForwardLinkID, int *d_BackwardLinkID, int *d_ConsecutiveNum, float PixelSize, float QE, int FluoNum)
{
	//    int tid = threadIdx.x;
	int gid = threadIdx.x + blockDim.x*blockIdx.x;


	int CurId = 0;
	int NextValid = 0;

	int OverlapNum = 0;


	if (gid < FluoNum)
	{
		// the first fluo is valid fluo and don't have a forward
		if (d_ForwardLinkID[gid] == 0)
		{
			CurId = gid;
			NextValid = 1;

			while (NextValid)
			{
				OverlapNum++;

				CurId = d_BackwardLinkID[CurId];
				NextValid = CurId > 0;
			}

			atomicAdd(&d_ConsecutiveNum[gid], OverlapNum);
		}
	}
}

__global__ void gpuGetOntimeDistribution(float *d_LocArry, int *d_ConsecutiveNum, int *d_OntimeDistrib, int *d_ValidFluoNum, int StartFrame, int EndFrame, int FluoNum)
{
	//    int tid = threadIdx.x;
	int gid = threadIdx.x + blockDim.x*blockIdx.x;
	float(*pLocArry)[OutParaNumGS2D]; // for 9 parameter array

	pLocArry = (float(*)[OutParaNumGS2D])d_LocArry;

	int CurFrame = 0;
	int ConsecutiveNum = 0;

	if (gid < FluoNum)
	{
		CurFrame = pLocArry[gid][Pos_Frme];

		if ((CurFrame >= StartFrame) && (CurFrame <= EndFrame))
		{
			ConsecutiveNum = d_ConsecutiveNum[gid];

			// valid molecule
			if (ConsecutiveNum > 0)
			{
				ConsecutiveNum--; // 1-0,2-1
				if (ConsecutiveNum < MaxOnTimeConsecutiveNum)
				{
					atomicAdd(&d_OntimeDistrib[ConsecutiveNum], 1);
				}
				atomicAdd(d_ValidFluoNum, 1);

			}
		}
	}
}



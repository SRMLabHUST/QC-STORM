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

#include "BackgroundNoiseFilter.h"




__global__ void gpuNeighborNumberCalc_NoiseFilter(float *d_LocArry, int* d_NeighborNum_Th1, int* d_NeighborNum_Th2, int* d_NeighborNum_Th3, float Distance_Th1, float Distance_Th2, float Distance_Th3, int FluoNum)
{
	//    int tid = threadIdx.x;
	int gid = blockDim.x*blockIdx.x + threadIdx.x;


	float(*pLocArry)[OutParaNumGS2D]; // for N parameter array
	pLocArry = (float(*)[OutParaNumGS2D])d_LocArry;

	Distance_Th1 = Distance_Th1*Distance_Th1;
	Distance_Th2 = Distance_Th2*Distance_Th2;
	Distance_Th3 = Distance_Th3*Distance_Th3;


	int CurId = gid;
	int NewID = 0;


	if (gid < FluoNum) // FluoNum
	{
		float CurDistance2;

		// unit of x, y is pixel in storage, while z is nm
		float CXPos = pLocArry[CurId][Pos_XPos];
		float CYPos = pLocArry[CurId][Pos_YPos];

		// find the nearest neighboring distance of current molecule
		if ((CXPos > 1) && (CYPos > 1))
		{
			for (int cnt = 0; cnt < FluoNum; cnt++)
			{
				NewID = cnt;

				float NXPos = pLocArry[NewID][Pos_XPos];
				float NYPos = pLocArry[NewID][Pos_YPos];

				// avoid the min distance of comparing to itself
				CurDistance2 = (CXPos - NXPos)*(CXPos - NXPos) + (CYPos - NYPos)*(CYPos - NYPos); // + (CZPos - NZPos)*(CZPos - NZPos) 

				atomicAdd(&d_NeighborNum_Th1[gid], (int)(CurDistance2 <= Distance_Th1));
				atomicAdd(&d_NeighborNum_Th2[gid], (int)(CurDistance2 <= Distance_Th2));
				atomicAdd(&d_NeighborNum_Th3[gid], (int)(CurDistance2 <= Distance_Th3));

			}
		}
	}
}

__global__ void gpuNoiseIdenfity_NoiseFilter(int* d_IsNoise, int* d_NeighborNum_Th1, int* d_NeighborNum_Th2, int* d_NeighborNum_Th3, int FluoNum)
{
	int gid = blockDim.x*blockIdx.x + threadIdx.x;

	if (gid < FluoNum)
	{
		float FluoNum_th1 = d_NeighborNum_Th1[gid];
		float FluoNum_th2 = d_NeighborNum_Th2[gid];
		float FluoNum_th3 = d_NeighborNum_Th3[gid];

		if (FluoNum_th1 < 1)FluoNum_th1 = 1;
		if (FluoNum_th2 < 1)FluoNum_th2 = 1;

		float Ratio1 = FluoNum_th2 / FluoNum_th1;
		float Ratio2 = FluoNum_th3 / FluoNum_th1;
		float Ratio3 = FluoNum_th3 / FluoNum_th2;

		int pos1 = Ratio1 <= 1;
		int pos2 = Ratio2 <= 1.3;
		int pos3 = Ratio3 <= 1;
		int pos4 = FluoNum_th2 <= 2;

		int IsNoise = pos1 || pos2 || pos3 || pos4;

		d_IsNoise[gid] = (int)(IsNoise > 0);
	}
}


__global__ void gpuRemoveNoiseFluo_NoiseFilter(float *d_LocArry, int* d_IsNoise, int FluoNum)
{
	int gid = blockDim.x*blockIdx.x + threadIdx.x;

	float(*pLocArry)[OutParaNumGS2D]; //
	pLocArry = (float(*)[OutParaNumGS2D])d_LocArry;

	if (gid < FluoNum)
	{
		if (d_IsNoise[gid] > 0)
		{
#pragma unroll
			for (int cnt = 0; cnt < OutParaNumGS2D; cnt++)
			{
				pLocArry[gid][cnt] = 0;
			}
		}
	}
}


__global__ void gpuMinDistanceCalc_NoiseFilter(float *d_LocArry, float *d_MinDistance, int SelFluoNum, int FluoNum)
{
	//    int tid = threadIdx.x;
	int gid = blockDim.x*blockIdx.x + threadIdx.x;


	float(*pLocArry)[OutParaNumGS2D]; // for N parameter array
	pLocArry = (float(*)[OutParaNumGS2D])d_LocArry;


	if (gid < SelFluoNum) // FluoNum
	{
		float MinDistanceTh2 = 1000000.0f;
		float CurDistance2;

		// unit of x, y is pixel in storage, while z is nm
		float CXPos = pLocArry[gid][Pos_XPos];
		float CYPos = pLocArry[gid][Pos_YPos];


		// find the nearest neighboring distance of current molecule
		if ((CXPos > 5) && (CYPos > 5))
		{
			for (int cnt = 0; cnt < FluoNum; cnt++)
			{

				float NXPos = pLocArry[cnt][Pos_XPos];
				float NYPos = pLocArry[cnt][Pos_YPos];

				// avoid the min distance of comparing to itself
				CurDistance2 = (CXPos - NXPos)*(CXPos - NXPos) + (CYPos - NYPos)*(CYPos - NYPos); // + (CZPos - NZPos)*(CZPos - NZPos) 
				CurDistance2 = CurDistance2 + (cnt == gid)*1000000.0f;

				MinDistanceTh2 = fminf(CurDistance2, MinDistanceTh2);
			}

			MinDistanceTh2 = __fsqrt_rn(MinDistanceTh2);

			d_MinDistance[gid] = MinDistanceTh2;
		}
	}
}

__global__ void gpuMeanDistanceCalc_NoiseFilter(float *d_MinDistance, float *d_ValidNum, float *d_TotalValue, float DistanceTh_pixel, int FluoNum)
{
	//    int tid = threadIdx.x;
	int gid = blockDim.x*blockIdx.x + threadIdx.x;


	float CurDistance;

	if (gid < FluoNum)
	{
		CurDistance = d_MinDistance[gid];

		if ((CurDistance > 0) && (CurDistance < DistanceTh_pixel))
		{
			atomicAdd(d_ValidNum, 1.0f);
			atomicAdd(d_TotalValue, CurDistance);

		}
	}
}

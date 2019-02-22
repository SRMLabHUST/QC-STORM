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

#include "cuda_runtime.h"
#include "device_launch_parameters.h"

#include "DimensionLocDensityCalc_Para.h"

__global__ void gpuMinDistanceCalc_2D(float *d_LocArry, float *d_MinDistance, float PixelSize_nm, int CalcStartPos, int SelDatLen, int FluoNum)
{
	//    int tid = threadIdx.x;
	int gid = blockDim.x*blockIdx.x + threadIdx.x;


	float(*pLocArry)[OutParaNumGS2D]; // for N parameter array
	pLocArry = (float(*)[OutParaNumGS2D])d_LocArry;


	int CurId = gid + CalcStartPos;
	int NewID = 0;


	if (gid < SelDatLen) // FluoNum
	{
		float MinDistanceTh2 = 1000000.0f;
		float CurDistance2;

		// unit of x, y is pixel in storage, while z is nm
		float CXPos = pLocArry[CurId][Pos_XPos];
		float CYPos = pLocArry[CurId][Pos_YPos];

		// find the nearest neighboring distance of current molecule
		if ((CXPos > 5) && (CYPos > 5))
		{
			for (int cnt = 0; cnt < FluoNum; cnt++)
			{
				NewID = cnt + CalcStartPos;

				float NXPos = pLocArry[NewID][Pos_XPos];
				float NYPos = pLocArry[NewID][Pos_YPos];

				// avoid the min distance of comparing to itself
				CurDistance2 = (CXPos - NXPos)*(CXPos - NXPos) + (CYPos - NYPos)*(CYPos - NYPos) + (NewID == CurId)*1000000.0f; // + (CZPos - NZPos)*(CZPos - NZPos) 


				MinDistanceTh2 = fminf(CurDistance2, MinDistanceTh2);
			}

			MinDistanceTh2 = __fsqrt_rn(MinDistanceTh2);

			d_MinDistance[gid] = MinDistanceTh2;
		}
	}
}

__global__ void gpuMinDistanceCalc_3D(float *d_LocArry, float *d_MinDistance, float PixelSize_nm, int CalcStartPos, int SelDatLen, int FluoNum)
{
	//    int tid = threadIdx.x;
	int gid = blockDim.x*blockIdx.x + threadIdx.x;


	float(*pLocArry)[OutParaNumGS2D]; // for N parameter array
	pLocArry = (float(*)[OutParaNumGS2D])d_LocArry;


	int CurId = gid + CalcStartPos;
	int NewID = 0;


	if (gid < SelDatLen) // FluoNum
	{
		float MinDistanceTh2 = 1000000.0f;
		float CurDistance2;

		// unit of x, y is pixel in storage, while z is nm
		float CXPos = pLocArry[CurId][Pos_XPos];
		float CYPos = pLocArry[CurId][Pos_YPos];
		float CZPos = pLocArry[CurId][Pos_ZPos] / PixelSize_nm; // convert the unit of z as the same with x,y


																// find the nearest neighboring distance of current molecule
		if ((CXPos > 5) && (CYPos > 5))
		{
			for (int cnt = 0; cnt < FluoNum; cnt++)
			{
				NewID = cnt + CalcStartPos;

				float NXPos = pLocArry[NewID][Pos_XPos];
				float NYPos = pLocArry[NewID][Pos_YPos];
				float NZPos = pLocArry[NewID][Pos_ZPos] / PixelSize_nm; // convert the unit of z as the same with x,y

																		// avoid the min distance of comparing to itself
				CurDistance2 = (CXPos - NXPos)*(CXPos - NXPos) + (CYPos - NYPos)*(CYPos - NYPos) + (CZPos - NZPos)*(CZPos - NZPos) + (NewID == CurId)*1000000.0f; // 


				MinDistanceTh2 = fminf(CurDistance2, MinDistanceTh2);
			}

			MinDistanceTh2 = __fsqrt_rn(MinDistanceTh2);

			d_MinDistance[gid] = MinDistanceTh2;
		}
	}
}


__global__ void gpuMeanDistanceCalc(float *d_MinDistance, float *d_ValidNum, float *d_TotalValue, float DistanceTh_pixel, int FluoNum)
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



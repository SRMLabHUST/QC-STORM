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

#include "statisticsInfo.h"

__global__ void FillerMoleculesInCurFrame(float *d_LocArry, float*d_FilteredLocArry, int *d_FilteredFluoNum, int FilterFrame, int LocType, int ImageWidth, int ImageHigh, int FluoNum)
{
	//    int tid = threadIdx.x;
	int gid = threadIdx.x + blockDim.x*blockIdx.x;


	float(*pLocArry)[OutParaNumGS2D]; // for parameter array

	float(*poLocArry10)[OutParaNumGS2D]; // for parameter array

	int CurFrame;
	int xpos;
	int ypos;

	pLocArry = (float(*)[OutParaNumGS2D])d_LocArry;

	poLocArry10 = (float(*)[OutParaNumGS2D])d_FilteredLocArry;

	int CurPos = 0;
	int cnt = 0;

	if (gid < FluoNum)
	{
		xpos = (int)(pLocArry[gid][Pos_XPos]); // pixel
		ypos = (int)(pLocArry[gid][Pos_YPos]); // pixel
		CurFrame = pLocArry[gid][Pos_Frme]; // 3d loc data


		if ((xpos >= 5) && (xpos < ImageWidth - 5) && ((ypos >= 5) && (ypos < ImageHigh - 5)))
		{
			if (CurFrame == FilterFrame)
			{
				CurPos = atomicAdd(d_FilteredFluoNum, 1);

#pragma unroll
				for (cnt = 0; cnt < OutParaNumGS2D; cnt++)
				{
					poLocArry10[CurPos][cnt] = pLocArry[gid][cnt];
				}

			}
		}
	}
}


__global__ void CalcOverlapFluoNum(float* d_LocArry, int *d_FluoNum_n0, int *d_FluoNum_n1, float RadiusTh_pixel, int FluoNum)
{
	int gid = threadIdx.x + blockDim.x*blockIdx.x;

	float(*pLocArry)[OutParaNumGS2D]; // for parameter array
	pLocArry = (float(*)[OutParaNumGS2D])d_LocArry;

	float CurXPos = pLocArry[gid][Pos_XPos];
	float CurYPos = pLocArry[gid][Pos_YPos];


	float DistanceTh2 = RadiusTh_pixel*RadiusTh_pixel;

	
	int cnt;
	float ComXPos;
	float ComYPos;

	float CurDistance2 = 0;
	int NeighborFluoNum = 0;

	if (gid < FluoNum)
	{
		for (cnt = 0; cnt < FluoNum; cnt++)
		{
			ComXPos = pLocArry[cnt][Pos_XPos];
			ComYPos = pLocArry[cnt][Pos_YPos];

			CurDistance2 = (ComXPos - CurXPos)*(ComXPos - CurXPos) + (ComYPos - CurYPos)*(ComYPos - CurYPos) + (cnt == gid)*1000000.0f;

			NeighborFluoNum += (CurDistance2 <= DistanceTh2);

		}

		// one is matched with itself
		// note that two or higher are matched mutually
		if (NeighborFluoNum == 0)atomicAdd(d_FluoNum_n0, 1);
		if (NeighborFluoNum == 1)atomicAdd(d_FluoNum_n1, 1);

	}

}



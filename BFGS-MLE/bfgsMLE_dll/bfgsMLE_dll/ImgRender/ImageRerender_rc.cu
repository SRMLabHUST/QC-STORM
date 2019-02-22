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

#include "ImageRender_top.h"


__global__ void FindMaxImgSize(float *d_LocArry, int *d_MaxImageWidth, int *d_MaxImageHigh, int FluoNum, int FillParaNum)
{
	//    int tid = threadIdx.x;
	int gid = threadIdx.x + blockDim.x*blockIdx.x;


	float(*pLocArry)[OutParaNumGS2D]; // for parameter array

	int xpos = 0;
	int ypos = 0;


	int HistPos = 0;
	int rcnt = 0;
	int Offset;

	pLocArry = (float(*)[OutParaNumGS2D])d_LocArry;


	if (gid<FluoNum)
	{
		xpos = (int)(pLocArry[gid][Pos_XPos] + 16); //
		ypos = (int)(pLocArry[gid][Pos_YPos] + 16); //


		atomicMax(d_MaxImageWidth, xpos);
		atomicMax(d_MaxImageHigh, ypos);
	}

}

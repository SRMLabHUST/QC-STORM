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

#include "MarginDataFilter.h"

__global__ void gpuFilterMarginData(float *d_LocArry, int FluoNum, float MarginPercentage, int ImageWidth, int ImageHigh);


void MarginDataFilter_TypeDef::FilterMarginData(float* ih_LocArry, int FluoNum, float MarginPercentage, int ImageWidth, int ImageHigh, cudaStream_t cstream)
{
	ValidFluoNum = FluoNum;

	int BlockSize = (FluoNum + ThreadsPerBlock - 1) / ThreadsPerBlock;

	cudaMemcpyAsync(d_LocArry, ih_LocArry, FluoNum*OutParaNumGS2D*sizeof(float), cudaMemcpyHostToDevice, cstream);

	gpuFilterMarginData << <BlockSize, ThreadsPerBlock, 0, cstream >> >(d_LocArry, FluoNum, MarginPercentage, ImageWidth, ImageHigh);

	cudaMemcpyAsync(h_LocArry, d_LocArry, FluoNum*OutParaNumGS2D*sizeof(float), cudaMemcpyDeviceToHost, cstream);

	cudaStreamSynchronize(cstream);

}


void MarginDataFilter_TypeDef::Init()
{
	ValidFluoNum = 0;

	cudaMallocHost((void **)&h_LocArry, MaxPointNum*OutParaNumGS2D*sizeof(float));

	cudaMalloc((void **)&d_LocArry, MaxPointNum*OutParaNumGS2D*sizeof(float));

}

void MarginDataFilter_TypeDef::DeInit()
{
	cudaFreeHost(h_LocArry);
	cudaFree(d_LocArry);

}


__global__ void gpuFilterMarginData(float *d_LocArry, int FluoNum, float MarginPercentage, int ImageWidth, int ImageHigh)
{
	//    int tid = threadIdx.x;
	int gid = threadIdx.x + blockDim.x*blockIdx.x;


	float(*pLocArry)[OutParaNumGS2D]; // for N parameter array
	pLocArry = (float(*)[OutParaNumGS2D])d_LocArry;


	float XL = ImageWidth*MarginPercentage;
	float XU = ImageWidth*(1 - MarginPercentage);

	float YL = ImageHigh*MarginPercentage;
	float YU = ImageHigh*(1 - MarginPercentage);

	float CurX = 0;
	float CurY = 0;

	if (gid < FluoNum)
	{
		CurX = pLocArry[gid][Pos_XPos];
		CurY = pLocArry[gid][Pos_YPos];

		if (!((CurX >= XL) && (CurX <= XU) && (CurY >= YL) && (CurY <= YU)))
		{
#pragma unroll
			for (int cnt = 0; cnt < OutParaNumGS2D; cnt++)
			{
				pLocArry[gid][cnt] = 0.0f;
			}
		}

	}

}


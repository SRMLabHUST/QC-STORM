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

#include "ZeroLocalizationsRemove.h"




void ZeroLocalizationsRemovel_TypeDef::RemoveZeroLocalizations(float *ih_LocArry, int iFluoNum, int SortFrameEn, int FirstFrame, int EndFrame, cudaStream_t cstream)
{
	cudaMemcpyAsync(d_LocArry_Raw, ih_LocArry, iFluoNum * OutParaNumGS2D * sizeof(float), cudaMemcpyHostToDevice, cstream);

	if (SortFrameEn == 0)
	{
		RemoveZeroWithoutSortFrame(ih_LocArry, iFluoNum, cstream);
	}
	else
	{
		RemoveZero_SortFrame(ih_LocArry, iFluoNum, FirstFrame, EndFrame, cstream);
	}

	cudaMemcpyAsync(h_LocArry, d_LocArry_Valid, ValidFluoNum * OutParaNumGS2D * sizeof(float), cudaMemcpyDeviceToHost, cstream);

	cudaStreamSynchronize(cstream);


}


void ZeroLocalizationsRemovel_TypeDef::RemoveZero_SortFrame(float *ih_LocArry, int iFluoNum, int FirstFrame, int EndFrame, cudaStream_t cstream)
{

	cudaMemsetAsync(d_FluoNum_Valid, 0, sizeof(int), cstream);

	int BlockDim = ThreadsPerBlock;
	int BlockNum = ((iFluoNum + ThreadsPerBlock - 1) / ThreadsPerBlock);

	for (int CurFrame = FirstFrame; CurFrame <= EndFrame; CurFrame++)
	{
		
		gpuRemoveZeorLoc_SortFrame << <BlockNum, BlockDim, 0, cstream >> > (d_LocArry_Valid, d_LocArry_Raw, d_FluoNum_Valid, CurFrame, iFluoNum);
	}

	cudaMemcpyAsync(h_FluoNum_Valid, d_FluoNum_Valid, sizeof(int), cudaMemcpyDeviceToHost, cstream);
	cudaStreamSynchronize(cstream);

	ValidFluoNum = *h_FluoNum_Valid;

}


void ZeroLocalizationsRemovel_TypeDef::RemoveZeroWithoutSortFrame(float *ih_LocArry, int iFluoNum, cudaStream_t cstream)
{
	// also find ValidFluoNum
	FindCopyID(ih_LocArry, iFluoNum);


	cudaMemcpyAsync(d_FluoID_Valid, h_FluoID_Valid, ValidFluoNum * sizeof(int), cudaMemcpyHostToDevice, cstream);

	int BlockDim = ThreadsPerBlock;
	int BlockNum = ((ValidFluoNum + ThreadsPerBlock - 1) / ThreadsPerBlock);

	gpuRemoveZeorLocalizations << <BlockNum, BlockDim, 0, cstream >> > (d_LocArry_Valid, d_LocArry_Raw, d_FluoID_Valid, ValidFluoNum);


}


void ZeroLocalizationsRemovel_TypeDef::FindCopyID(float *ih_LocArry, int iFluoNum)
{
	float(*pLocArry)[OutParaNumGS2D] = (float(*)[OutParaNumGS2D])ih_LocArry;

	int ValidId = 0;
	float PeakPhoton = 0;

	for (int CurID = 0; CurID < iFluoNum; CurID++)
	{
		PeakPhoton = pLocArry[CurID][Pos_PPho];

		if (PeakPhoton > 1.0f)
		{
			h_FluoID_Valid[ValidId] = CurID;
			ValidId++;
		}
	}

	ValidFluoNum = ValidId;

//	printf("ValidFluoNum:%d \n", ValidFluoNum);
}

void ZeroLocalizationsRemovel_TypeDef::Init()
{
	ValidFluoNum = 0;


	// host and gpu
	cudaMallocHost((void **)&h_LocArry, MaxPointNum*OutParaNumGS2D * sizeof(float));

	cudaMalloc((void **)&d_LocArry_Raw, MaxPointNum*OutParaNumGS2D * sizeof(float));
	cudaMalloc((void **)&d_LocArry_Valid, MaxPointNum*OutParaNumGS2D * sizeof(float));

	cudaMallocHost((void **)&h_FluoID_Valid, MaxPointNum * sizeof(int));

	cudaMalloc((void **)&d_FluoID_Valid, MaxPointNum * sizeof(int));

	cudaMallocHost((void **)&h_FluoNum_Valid, sizeof(int));
	cudaMalloc((void **)&d_FluoNum_Valid, sizeof(int));

}

void ZeroLocalizationsRemovel_TypeDef::Deinit()
{
	cudaFreeHost(h_LocArry);

	cudaFree(d_LocArry_Raw);
	cudaFree(d_LocArry_Valid);

	cudaFreeHost(h_FluoID_Valid);

	cudaFree(d_FluoID_Valid);

	cudaFreeHost(h_FluoNum_Valid);
	cudaFree(d_FluoNum_Valid);

}



__global__ void gpuRemoveZeorLocalizations(float * d_LocArry_Valid, float * d_LocArry_Raw, int *d_FluoID_Valid, int ValidFluoNum)
{
	int gid = blockDim.x*blockIdx.x + threadIdx.x;

	float(*pLocArry_Raw)[OutParaNumGS2D] = (float(*)[OutParaNumGS2D])d_LocArry_Raw;
	float(*pLocArry_Valid)[OutParaNumGS2D] = (float(*)[OutParaNumGS2D])d_LocArry_Valid;

	if (gid < ValidFluoNum)
	{
		int RawDataID = d_FluoID_Valid[gid];

#pragma unroll
		for (int cnt = 0; cnt < OutParaNumGS2D; cnt++)
		{
			pLocArry_Valid[gid][cnt] = pLocArry_Raw[RawDataID][cnt];
		}
	}
}

__global__ void gpuRemoveZeorLoc_SortFrame(float * d_LocArry_Valid, float * d_LocArry_Raw, int *d_FluoNum_Valid, int CurFrame, int FluoNum)
{
	int gid = blockDim.x*blockIdx.x + threadIdx.x;

	float(*pLocArry_Raw)[OutParaNumGS2D] = (float(*)[OutParaNumGS2D])d_LocArry_Raw;
	float(*pLocArry_Valid)[OutParaNumGS2D] = (float(*)[OutParaNumGS2D])d_LocArry_Valid;

	if (gid < FluoNum)
	{

		int PeakPhoton = pLocArry_Raw[gid][Pos_PPho];
		int Frame = pLocArry_Raw[gid][Pos_Frme];

		if ((PeakPhoton > 1.0f) && (Frame == CurFrame))
		{

			int NewStoreID = atomicAdd(d_FluoNum_Valid, 1);

#pragma unroll
			for (int cnt = 0; cnt < OutParaNumGS2D; cnt++)
			{
				pLocArry_Valid[NewStoreID][cnt] = pLocArry_Raw[gid][cnt];
			}
		}
	}
}


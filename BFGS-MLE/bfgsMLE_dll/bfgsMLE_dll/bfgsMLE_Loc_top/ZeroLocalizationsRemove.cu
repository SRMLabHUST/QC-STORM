#include "ZeroLocalizationsRemove.h"




void ZeroLocalizationsRemovel_TypeDef::RemoveZeroLocalizations(float *ih_LocArry, int iFluoNum, cudaStream_t cstream)
{
	FindCopyID(ih_LocArry, iFluoNum);


	cudaMemcpyAsync(d_LocArry_Raw, ih_LocArry, iFluoNum * OutParaNumGS2D * sizeof(float), cudaMemcpyHostToDevice, cstream);
	cudaMemcpyAsync(d_FluoID_Valid, h_FluoID_Valid, ValidFluoNum * sizeof(int), cudaMemcpyHostToDevice, cstream);

	int BlockDim = ThreadsPerBlock;
	int BlockNum = ((ValidFluoNum + ThreadsPerBlock - 1) / ThreadsPerBlock);

	gpuRemoveZeorLocalizations << <BlockNum, BlockDim, 0, cstream >> > (d_LocArry_Valid, d_LocArry_Raw, d_FluoID_Valid, ValidFluoNum);

	cudaMemcpyAsync(h_LocArry, d_LocArry_Valid, ValidFluoNum * OutParaNumGS2D * sizeof(float), cudaMemcpyDeviceToHost, cstream);

	cudaStreamSynchronize(cstream);

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

}

void ZeroLocalizationsRemovel_TypeDef::Deinit()
{
	cudaFreeHost(h_LocArry);

	cudaFree(d_LocArry_Raw);
	cudaFree(d_LocArry_Valid);

	cudaFreeHost(h_FluoID_Valid);

	cudaFree(d_FluoID_Valid);

}



__global__ void gpuRemoveZeorLocalizations(float * d_LocArry_Valid, float * d_LocArry_Raw, int *d_FluoID_Valid, int ValidFluoNum)
{
	int gid = blockDim.x*blockIdx.x + threadIdx.x;

	float(*pLocArry_Valid)[OutParaNumGS2D] = (float(*)[OutParaNumGS2D])d_LocArry_Valid;
	float(*pLocArry_Raw)[OutParaNumGS2D] = (float(*)[OutParaNumGS2D])d_LocArry_Raw;

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


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


void BackgroundNoiseFilter_TypeDef::FilterBackgroundNoise(float * h_iLocArry, int FluoNum, int DataSource, cudaStream_t cstream)
{

	int BlockDim = ThreadsPerBlock;
	int BlockNum = (FluoNum + ThreadsPerBlock - 1) / ThreadsPerBlock;

	if (DataSource == ImageSource_CPU_Pinned)
	{
		cudaMemcpyAsync(d_LocArry, h_iLocArry, FluoNum*OutParaNumGS2D * sizeof(float), cudaMemcpyHostToDevice, cstream);
	}
	else if(DataSource == ImageSource_CPU_Normal)
	{
		cudaMemcpy(d_LocArry, h_iLocArry, FluoNum*OutParaNumGS2D * sizeof(float), cudaMemcpyHostToDevice);
	}
	else if (DataSource == ImageSource_GPU)
	{
		cudaMemcpyAsync(d_LocArry, h_iLocArry, FluoNum*OutParaNumGS2D * sizeof(float), cudaMemcpyDeviceToDevice, cstream);
	}
	else
	{
		return;
	}


	BackgroundNoiseRemove(FluoNum, cstream);

	if (DataSource == ImageSource_CPU_Pinned)
	{
		cudaMemcpyAsync(h_LocArry, d_LocArry, FluoNum*OutParaNumGS2D * sizeof(float), cudaMemcpyDeviceToHost, cstream);
	}
	else if (DataSource == ImageSource_CPU_Normal)
	{
		cudaMemcpy(h_LocArry, d_LocArry, FluoNum*OutParaNumGS2D * sizeof(float), cudaMemcpyDeviceToHost);
	}


	cudaStreamSynchronize(cstream);

}

void BackgroundNoiseFilter_TypeDef::BackgroundNoiseRemove(int FluoNum, cudaStream_t cstream)
{
	int BlockDim = ThreadsPerBlock;
	int BlockNum = (FluoNum + ThreadsPerBlock - 1) / ThreadsPerBlock;

	NeighborNumCalc(FluoNum, cstream);

	gpuNoiseIdenfity_NoiseFilter << <BlockNum, BlockDim, 0, cstream >> > (d_IsNoise, d_NeighborNum_Th1, d_NeighborNum_Th2, d_NeighborNum_Th3, FluoNum);
	gpuRemoveNoiseFluo_NoiseFilter << <BlockNum, BlockDim, 0, cstream >> > (d_LocArry, d_IsNoise, FluoNum);

	cudaStreamSynchronize(cstream);

}

void BackgroundNoiseFilter_TypeDef::NeighborNumCalc(int FluoNum, cudaStream_t cstream)
{
	int BlockDim = ThreadsPerBlock;
	int BlockNum = (FluoNum + ThreadsPerBlock - 1) / ThreadsPerBlock;

	float MeanDistance = GetMeanMinNearestNeighborDistance(FluoNum, cstream);

	float Distance_Th1 = MeanDistance * 4.0f;
	float Distance_Th2 = Distance_Th1 * 2.5;
	float Distance_Th3 = Distance_Th1 * 5;

	cudaMemsetAsync(d_NeighborNum_Th1, 0, FluoNum * sizeof(int), cstream);
	cudaMemsetAsync(d_NeighborNum_Th2, 0, FluoNum * sizeof(int), cstream);
	cudaMemsetAsync(d_NeighborNum_Th3, 0, FluoNum * sizeof(int), cstream);

	gpuNeighborNumberCalc_NoiseFilter << <BlockNum, BlockDim, 0, cstream >> > (d_LocArry, d_NeighborNum_Th1, d_NeighborNum_Th2, d_NeighborNum_Th3, Distance_Th1, Distance_Th2, Distance_Th3, FluoNum);

	cudaStreamSynchronize(cstream);

}


float BackgroundNoiseFilter_TypeDef::GetMeanMinNearestNeighborDistance(int FluoNum, cudaStream_t cstream)
{
	// get min nearest neighboring distance
	int SelFluoNum = min(25000, FluoNum); // 

	GetMinNearestNeighborDistance(FluoNum, SelFluoNum, cstream);

	int BlockDim = ThreadsPerBlock;
	int BlockNum = (SelFluoNum + ThreadsPerBlock - 1) / ThreadsPerBlock;

	// get mean min nearest neighboring distance

	// calculate mean distance again without unnormal molecules filtered
	cudaMemsetAsync(d_ValidNum, 0, sizeof(float), cstream);
	cudaMemsetAsync(d_TotalValue, 0, sizeof(int), cstream);

	float Distance_Th_Pixel = 400;
	gpuMeanDistanceCalc_NoiseFilter << <BlockNum, BlockDim, 0, cstream >> > (d_MinDistance, d_ValidNum, d_TotalValue, Distance_Th_Pixel, SelFluoNum);

	//
	cudaMemcpyAsync(h_ValidNum, d_ValidNum,  sizeof(float), cudaMemcpyDeviceToHost, cstream);
	cudaMemcpyAsync(h_TotalValue, d_TotalValue, sizeof(float), cudaMemcpyDeviceToHost, cstream);
	cudaStreamSynchronize(cstream);

	if (*h_ValidNum < 1)*h_ValidNum = 1;
	float MeanDistance = *h_TotalValue / *h_ValidNum;

//	printf("mean distance:%f %f %f\n", MeanDistance, *h_TotalValue, *h_ValidNum);

	// calculate mean distance again with unnormal molecules filtered
	cudaMemsetAsync(d_ValidNum, 0, sizeof(float), cstream);
	cudaMemsetAsync(d_TotalValue, 0, sizeof(int), cstream);

	Distance_Th_Pixel = MeanDistance * 4;
	gpuMeanDistanceCalc_NoiseFilter << <BlockNum, BlockDim, 0, cstream >> > (d_MinDistance, d_ValidNum, d_TotalValue, Distance_Th_Pixel, SelFluoNum);

	//
	cudaMemcpyAsync(h_ValidNum, d_ValidNum, sizeof(float), cudaMemcpyDeviceToHost, cstream);
	cudaMemcpyAsync(h_TotalValue, d_TotalValue, sizeof(float), cudaMemcpyDeviceToHost, cstream);
	cudaStreamSynchronize(cstream);

	if (*h_ValidNum < 1)*h_ValidNum = 1;
	MeanDistance = *h_TotalValue / *h_ValidNum;

//	printf("mean distance:%f %f %f\n", MeanDistance, *h_TotalValue, *h_ValidNum);

	return MeanDistance;

}


void BackgroundNoiseFilter_TypeDef::GetMinNearestNeighborDistance(int FluoNum, int SelFluoNum, cudaStream_t cstream)
{
	SelFluoNum = min(22000, FluoNum);

	int BlockDim = ThreadsPerBlock;
	int BlockNum = (SelFluoNum + ThreadsPerBlock - 1) / ThreadsPerBlock;

	cudaMemsetAsync(d_MinDistance, 0, SelFluoNum * sizeof(float), cstream);

	gpuMinDistanceCalc_NoiseFilter << <BlockNum, BlockDim, 0, cstream >> > (d_LocArry, d_MinDistance, SelFluoNum, FluoNum);
	cudaStreamSynchronize(cstream);

}



void BackgroundNoiseFilter_TypeDef::Init(int TotalFluoNum)
{
	cudaMallocHost((void **)&h_LocArry, TotalFluoNum * OutParaNumGS2D*sizeof(float));
	cudaMalloc((void **)&d_LocArry, TotalFluoNum * OutParaNumGS2D*sizeof(float));

	//
	cudaMalloc((void **)&d_MinDistance, TotalFluoNum * sizeof(float));

	cudaMallocHost((void **)&h_ValidNum, sizeof(float));
	cudaMallocHost((void **)&h_TotalValue, sizeof(float));

	cudaMalloc((void **)&d_ValidNum,   sizeof(float));
	cudaMalloc((void **)&d_TotalValue,   sizeof(float));

	//
	cudaMalloc((void **)&d_NeighborNum_Th1, TotalFluoNum * sizeof(int));
	cudaMalloc((void **)&d_NeighborNum_Th2, TotalFluoNum * sizeof(int));
	cudaMalloc((void **)&d_NeighborNum_Th3, TotalFluoNum * sizeof(int));

	cudaMalloc((void **)&d_IsNoise, TotalFluoNum * sizeof(int));

}


void BackgroundNoiseFilter_TypeDef::DeInit()
{
	cudaFreeHost(h_LocArry);
	cudaFree(d_LocArry);
	//
	cudaFree(d_MinDistance);

	cudaFreeHost(h_ValidNum);
	cudaFreeHost(h_TotalValue);

	cudaFree(d_ValidNum);
	cudaFree(d_TotalValue);

	//
	cudaFree(d_NeighborNum_Th1);
	cudaFree(d_NeighborNum_Th2);
	cudaFree(d_NeighborNum_Th3);

	cudaFree(d_IsNoise);

}



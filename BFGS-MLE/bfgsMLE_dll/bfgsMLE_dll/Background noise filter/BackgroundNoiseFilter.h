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

#pragma once


#include "cuda_runtime.h"
#include "device_launch_parameters.h"


#include "bfgs_CommonPara.h"


#define ThreadsPerBlock				32	//Threads Per Block



// filter background noised by automatic threshold
class BackgroundNoiseFilter_TypeDef
{
public:
	float * h_LocArry;
	float * d_LocArry;

private:
	// find distance threshold
	float *d_MinDistance;

	float *h_ValidNum;
	float *h_TotalValue;

	float *d_ValidNum;
	float *d_TotalValue;

	// filter noise based on neighbor numbe in threshold
	int* d_NeighborNum_Th1;
	int* d_NeighborNum_Th2;
	int* d_NeighborNum_Th3;

	int* d_IsNoise;

public:
	void Init(int TotalFluoNum); // create CPU&GPU memory
	void DeInit(); // release CPU&GPU memory

	// if the data source is cpu, data is copied out, otherwise, it's not
	void FilterBackgroundNoise(float * h_iLocArry, int FluoNum, int DataSource, cudaStream_t cstream);

private:
	void BackgroundNoiseRemove(int FluoNum, cudaStream_t cstream);
	void GetMinNearestNeighborDistance(int FluoNum, int SelFluoNum, cudaStream_t cstream);

	float GetMeanMinNearestNeighborDistance(int FluoNum, cudaStream_t cstream);
	void NeighborNumCalc(int FluoNum, cudaStream_t cstream);

};



__global__ void gpuNeighborNumberCalc_NoiseFilter(float *d_LocArry, int* d_NeighborNum_Th1, int* d_NeighborNum_Th2, int* d_NeighborNum_Th3, float Distance_Th1, float Distance_Th2, float Distance_Th3, int FluoNum);
__global__ void gpuNoiseIdenfity_NoiseFilter(int* d_IsNoise, int* d_NeighborNum_Th1, int* d_NeighborNum_Th2, int* d_NeighborNum_Th3, int FluoNum);
__global__ void gpuRemoveNoiseFluo_NoiseFilter(float *d_LocArry, int* d_IsNoise, int FluoNum);
__global__ void gpuMinDistanceCalc_NoiseFilter(float *d_LocArry, float *d_MinDistance, int SelFluoNum, int FluoNum);
__global__ void gpuMeanDistanceCalc_NoiseFilter(float *d_MinDistance, float *d_ValidNum, float *d_TotalValue, float DistanceTh_pixel, int FluoNum);


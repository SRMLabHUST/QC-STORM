#pragma once


#include "cuda_runtime.h"
#include "device_launch_parameters.h"


#include "bfgs_MLE_dll.h"


#define ThreadsPerBlock				32	//Threads Per Block



// filter molecules in consecutive frames by a radius threshold, only keep the molecule in the first frame
class ConsecutiveFluoMerger_TypeDef
{
public:
	float * h_LocArry;
	float * d_LocArry;


	//
	int *d_ForwardLinkID;
	int *d_BackwardLinkID;

public:
	void Init(unsigned int TotalFluoNum); // create CPU&GPU memory
	void DeInit(); // release CPU&GPU memory

	// Mode 0: keep first, 1: keep Max SNR
	void MergeConsecutiveFluo(float * h_LocArry, int FluoNum, LocalizationPara & LocPara, int FilterMode, float Distance_th_pixel, cudaStream_t cstream);


};



__global__ void gpuFindConsecutiveFilterPair(float *d_LocArry, int *d_ForwardLinkID, int *d_BackwardLinkID, float Distance_th_pixel, int FluoNum);

__global__ void gpuConsecutiveFit(float *d_LocArry, int *d_ForwardLinkID, int *d_BackwardLinkID, float QE, int FluoNum);

__global__ void gpuRemoveConsecutiveFluo_KeepFirst(float *d_LocArry, int *d_ForwardLinkID, int *d_BackwardLinkID, int FluoNum);


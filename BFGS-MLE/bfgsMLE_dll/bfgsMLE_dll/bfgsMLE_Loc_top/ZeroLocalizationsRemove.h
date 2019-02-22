#pragma once

#include <iostream>
using namespace std;


#include "cuda_runtime.h"
#include "device_launch_parameters.h"


#include "bfgs_CommonPara.h"

#define ThreadsPerBlock					32 //Threads Per Block


// remove invalid localizations with zero value
class ZeroLocalizationsRemovel_TypeDef {
public:
	int ValidFluoNum;

	float *h_LocArry;

public:
	void RemoveZeroLocalizations(float *ih_LocArry,int iFluoNum, cudaStream_t cstream);

	void Init();

	void Deinit();


private:
	void FindCopyID(float *ih_LocArry, int iFluoNum);

private:
	float *d_LocArry_Raw;
	float *d_LocArry_Valid;


	int *h_FluoID_Valid;
	int *d_FluoID_Valid;

};



__global__ void gpuRemoveZeorLocalizations(float * d_LocArry_Valid, float * d_LocArry_Raw, int *d_FluoID_Valid, int ValidFluoNum);



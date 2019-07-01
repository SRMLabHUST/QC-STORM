#pragma once
#include "cuda_runtime.h"
#include "device_launch_parameters.h"

#include <stdio.h>

#include "bfgs_base.h"


#define ThreadsPerBlock				32	//Threads Per Block



#define PAIR_ID_LEN			2

#define PAIR_ID_1st			0
#define PAIR_ID_2nd			1



// based on localized or even consecutive fitted position from bfgs_top
class DH3D_MoleculePair_TypeDef
{

public:
	float *h_LocArry;
	float *d_LocArry;

	float *h_oLocArry;
	float *d_oLocArry;

	int oValidFluoNum;


private:
	// paired molecule id for the same molecule
	int *h_PairID;
	int *d_PairID;

	int * h_ValidoFluoNum;
	int * d_ValidoFluoNum;

public:

	void Init();

	void Deinit();

	void MoleculePair(float *h_iLocArry, int FluoNum, LocalizationPara & LocPara, cudaStream_t cstream);

};




__global__ void gpu_MoleculePair(float *d_LocArry, int FluoNum, int *d_PairID, int * d_ValidoFluoNum, float MeanDistance, float DistanceVaryTh);

__global__ void gpu_MoleculeMerge(float *d_LocArry, float *d_oLocArry, int *d_PairID, int ValidFluoNum, int RotateMode, float p4, float p3, float p2, float p1, float p0);

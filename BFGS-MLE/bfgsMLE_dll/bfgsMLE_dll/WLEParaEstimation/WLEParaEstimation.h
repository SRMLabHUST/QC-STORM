#pragma once


#include "cuda_runtime.h"
#include "device_launch_parameters.h"

#include <stdio.h>


#include "bfgs_base.h"
#include "cudaWrapper.h"


#include "WLEParaEstimation_core.h"


#define ThreadsPerBlock				32 //Threads Per Block



// estimage WLE para, contained in the ROI extraction
class WLEParameterEstimation_TypeDef
{
public:
	unsigned short * d_ROIMem;

	float *h_WLEPara;
	float *d_WLEPara;

public:

	void Init(LocalizationPara & LocPara);
	void Deinit();

	void WLEParameterEstimate(unsigned short * h_ROIMem, int ROISize, int FluoNum, cudaStream_t cstream);
};




void CalculatePSFWidth(unsigned short * d_ROIMem, float *d_WLEPara, int FluoNum, int ROISize, cudaStream_t cstream);

void CalculateNearestNeighborDistance(unsigned short * d_ROIMem, int ROISize, float *d_WLEPara, int FluoNum, cudaStream_t cstream);

void MoleculeTypeClasify(int ROISize, float *d_WLEPara, int FluoNum, cudaStream_t cstream);

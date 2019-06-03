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

#include <stdio.h>


#include "bfgs_base.h"
#include "cudaWrapper.h"


#include "WLEParaEstimation_core.h"


#define ThreadsPerBlock				32 //Threads Per Block



// estimage WLE para, contained in the ROI extraction
class WLEParameterEstimation_TypeDef
{
public:

	float *h_WLEPara;
	float *d_WLEPara;

public:

	void Init(LocalizationPara & LocPara);
	void Deinit();

	void WLEParameterEstimate(unsigned short * d_ImageROI, int LocType, int MultiEmitterFitEn, int ROISize, int FluoNum, cudaStream_t cstream);
};




void CalculatePSFWidth(unsigned short * d_ImageROI, float *d_WLEPara, int FluoNum, int ROISize, cudaStream_t cstream);

void CalculateNearestNeighborDistance(unsigned short * d_ImageROI, int ROISize, float *d_WLEPara, int FluoNum, cudaStream_t cstream);

void MoleculeTypeClasify(unsigned short * d_ImageROI, int LocType, int MultiEmitterFitEn, int ROISize, float *d_WLEPara, int FluoNum, cudaStream_t cstream);

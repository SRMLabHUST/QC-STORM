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
#include "bfgs_base.h"

#include "cudaWrapper.h"


#include "bfgsMLE_Single_Emitter_Fit.h"
#include "bfgsMLE_Multi_Emitter_Fit.h"



// both 2d and 3d localization data structure
class LDLocData_TypeDef
{
public:
	unsigned short * h_ImageROI;
	unsigned short * d_ImageROI;

	float * h_LocArry;
	float * d_LocArry;

	float *d_WLEPara;


	// on time calculation
	float *h_OntimeRatio; // for display and time vatiation


	// valid number after localization, still include filtered molecule number
	int oValidFluoNum;

public:

	// multi emitter fitting
	int * h_MultiFitFluoNum;
	int * d_MultiFitFluoNum;
	int * d_MultiFitFluoPos; // position id

	float MultiFitRatio;
	
private:
	// on time calculation, to find Consecutive molecules in adjecent frames
	int *d_ForwardLinkID;
	int *d_BackwardLinkID;
	int *d_ConsecutiveNum;


	int *h_OntimeDistrib;
	int *d_OntimeDistrib;

	int *h_ValidFluoNum; // for ontime distribution 
	int *d_ValidFluoNum; // for ontime distribution 


	// for loc filter
	float *h_SNRSumUp;
	int *h_ValidNum;

	float *d_SNRSumUp;
	int *d_ValidNum;


public:
	void Init(LocalizationPara & LocPara); // create CPU&GPU memory
	void Deinit(LocalizationPara & LocPara); // release CPU&GPU memory

	void BFGS_MLELocalization(unsigned short * h_ImageROI, float *h_WLEPara, LocalizationPara & LocPara, int FluoNum, cudaStream_t cstream);

	
	void OntimeCalc(LocalizationPara & LocPara, int FluoNum, cudaStream_t cstream);

	void CopyDataToGPU(float * h_LocArry, int FluoNum, cudaStream_t cstream);

public:

	static int GetFirstFrame(float * h_LocArry, int FluoNum);
	static int GetLastFrame(float * h_LocArry, int FluoNum);

	static int GetFirstFrameFromROI(unsigned short * h_ImageROI, int ROISize, int FluoNum);
	static int GetLastFrameFromROI(unsigned short * h_ImageROI, int ROISize, int FluoNum);


	// two optional localization precision method, only suitable for 2d localization with symmetric Gaussian PSF
	static void LocPrecCalc_GaussianCRLB(float* d_LocArry, LocalizationPara & LocPara, int FluoNum, cudaStream_t cstream);

private:
	void FilterBadFit(LocalizationPara & LocPara, int FluoNum, cudaStream_t cstream);

};



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




// both 2d and 3d localization data structure
class LDLocData_TypeDef
{
public:
	unsigned short * h_SubRegion;
	unsigned short * d_SubRegion;

	float * h_LocArry;
	float * d_LocArry;

	float *d_WLEPara;


	// Consecutive finding from adjecent frames
	// same with consecutive fitting method
	int *d_ForwardLinkID;
	int *d_BackwardLinkID;
	int *d_ConsecutiveNum;

	int *h_OntimeDistrib; 
	int *d_OntimeDistrib; 

	int *h_ValidFluoNum; // for ontime distribution 
	int *d_ValidFluoNum; // for ontime distribution 

	float *h_OntimeRatio; // for display and time vatiation

	// valid number after localization, still include filtered molecule number
	int oValidFluoNum;

	
private:
	// for loc filter
	float *h_SNRSumUp;
	int *h_ValidNum;

	float *d_SNRSumUp;
	int *d_ValidNum;


public:
	void Init(LocalizationPara & LocPara); // create CPU&GPU memory
	void Deinit(LocalizationPara & LocPara); // release CPU&GPU memory

	void BFGS_MLELocalization(unsigned short * h_SubRegion, float *h_WLEPara, LocalizationPara & LocPara, int FluoNum, cudaStream_t cstream);

	
	void OntimeCalc(LocalizationPara & LocPara, int FluoNum, cudaStream_t cstream);

	void CopyDataToGPU(float * h_LocArry, int FluoNum, cudaStream_t cstream);

public:

	static int GetFirstFrame(float * h_LocArry, int FluoNum);
	static int GetLastFrame(float * h_LocArry, int FluoNum);

	static int GetFirstFrameFromROI(unsigned short * h_SubRegion, int ROISize, int FluoNum);
	static int GetLastFrameFromROI(unsigned short * h_SubRegion, int ROISize, int FluoNum);


	// two optional localization precision method, only suitable for 2d localization with symmetric Gaussian PSF
	static void LocPrecCalc_GaussianCRLB(float* d_LocArry, LocalizationPara & LocPara, int FluoNum, cudaStream_t cstream);

private:
	void FilterBadFit(LocalizationPara & LocPara, int FluoNum, cudaStream_t cstream);

};



// bfgs 2D and as3D loc function
void LDLoc_BFGS_MLELocalizationGS2D(float * d_LocArry, unsigned short * d_SubRegion, float *d_WLEPara, LocalizationPara& LocPara, int FluoNum, cudaStream_t cstream);

void LDLoc_BFGS_MLELocalizationAS3D(float * d_LocArry, unsigned short * d_SubRegion, float *d_WLEPara, LocalizationPara& LocPara, int FluoNum, cudaStream_t cstream);

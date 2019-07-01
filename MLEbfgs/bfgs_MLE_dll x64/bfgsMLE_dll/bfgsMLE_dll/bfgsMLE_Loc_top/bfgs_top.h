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

	float *d_WLEPara; // copy data from roi extraction


	float * h_LocArry;
	float * d_LocArry;


	// valid number after localization
	int oValidFluoNum;


public:

	// ratio of ROI to all detected ROIs, include molecules fitted by multi modality
	float FitRatio_1E; // single molecule fitting ratio
	float FitRatio_2E; // two emitter fitting ratio
	float FitRatio_3E; // three emitter fitting ratio

	// ratio of ROI by last fitting modality to all detected ROIs
	float FitRatio_Final_1E; // single molecule fitting ratio
	float FitRatio_Final_2E; // two emitter fitting ratio
	float FitRatio_Final_3E; // three emitter fitting ratio
	float FitRatio_Final_4E; // four or more emitter fitting ratio

private:

	// for convinient parameter transimission
	CoreFittingPara* h_FitPara;
	CoreFittingPara* d_FitPara;

	FitPosInf_TypeDef* h_FitPosInf;
	FitPosInf_TypeDef* d_FitPosInf;



public:
	void Init(LocalizationPara & LocPara); // create CPU&GPU memory
	void Deinit(LocalizationPara & LocPara); // release CPU&GPU memory

	void BFGS_MLELocalization(unsigned short * h_ImageROI, float *h_WLEPara, LocalizationPara & LocPara, int FluoNum, cudaStream_t cstream);

private:

	void CopyFittingPara(LocalizationPara & LocPara, cudaStream_t cstream);


	void MoleculePreFitClasify(int ROISize, int MultiEmitterFitEn, int FluoNum, cudaStream_t cstream);

	void ResetNumbers(cudaStream_t cstream);

public:
	// note the frame must be sorted by ZeroLocalizationsRemove class
	static int GetFirstFrame(float * h_LocArry, int FluoNum);
	static int GetLastFrame(float * h_LocArry, int FluoNum);


	// two optional localization precision method, only suitable for 2d localization with symmetric Gaussian PSF
	static void LocPrecCalc_GaussianCRLB(float* d_LocArry, LocalizationPara & LocPara, int FluoNum, cudaStream_t cstream);

};



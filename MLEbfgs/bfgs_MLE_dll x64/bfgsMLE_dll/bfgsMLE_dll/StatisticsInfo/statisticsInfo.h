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

#include <iostream>
#include <vector>

using namespace std;



#include <stdio.h>
#include <string.h>
#include <math.h>

#include "bfgs_CommonPara.h"


#include "bfgs_top.h"



// calculate the statistic information of localized data

#define ThreadsPerBlock					32 //Threads Per Block



#define StatInf_Hist_DatLen				400
#define StatInf_Hist_DatLenMax			(StatInf_Hist_DatLen*4)


#define Hist_MaxTotalPhoton				10000.0f
#define Hist_PhotonGap					(Hist_MaxTotalPhoton/StatInf_Hist_DatLen)


#define Hist_MaxLocPrecision			60.0f
#define Hist_LocPrecisionGap			(Hist_MaxLocPrecision/StatInf_Hist_DatLen)


#define Hist_MaxPeakSNR					50.0f
#define Hist_PeakSNRGap					(Hist_MaxPeakSNR/StatInf_Hist_DatLen)


#define Hist_MaxPSFWidth				5.0f
#define Hist_PSFWidthGap				(Hist_MaxPSFWidth/StatInf_Hist_DatLen)


#define MaxOnTimeConsecutiveNum			10


#define LocDensityCalc_RadiusTh_pixel_2D		10
#define LocDensityCalc_RadiusTh_pixel_3D		20


// stastical information class for both 2d and 3d
class FluoStatisticData_TypeDef
{
public:
	int FirstFrame;
	int EndFrame;

	// Distribution of current data
	int *h_Hist_TotalPhoton; 
	int *h_Hist_LocPrecisionXY;	
	int *h_Hist_PeakSNR;
	int *h_Hist_PSFWidth;

	// mean value of current data
	float MeanTotalPhoton;
	float MeanLocPrecisionXY; // sqrt(LocPrec_x^2+LocPrec_y^2)
	float MeanPeakSNR;
	float MeanBackground;
	float MeanPSFWidth;
	float MeanPSFWidth_Ctl; // h_MeanPSFWidth1 = MeanPSFWidth  * 10.0f / MeanPeakSNR1; // lower is better for feedback control
	float MeanLocDensity2D; // cur activation density

	// time variation curve data 

	vector<float> TimeVary_TotalPhoton;
	vector<float> TimeVary_LocPrecisionXY;
	vector<float> TimeVary_PeakSNR;
	vector<float> TimeVary_Background;
	vector<float> TimeVary_PSFWidth;
	vector<float> TimeVary_PSFWidth_Ctl;
	vector<float> TimeVary_LocDensity2D;


	// come from localization
	float *h_OntimeRatio; // for display and time vatiation
	float MeanOntime;


	vector<float> TimeVary_OntimeF1Ratio;
	vector<float> TimeVary_OntimeF2Ratio;
	vector<float> TimeVary_OntimeF3Ratio;
	vector<float> TimeVary_OntimeF4Ratio;


private:

	float *d_LocArry;

	// Distribution statistic info in gpu
	int *d_Hist_TotalPhoton; // total photon
	int *d_Hist_LocPrecisionXY; // theoretical precision
	int *d_Hist_PeakSNR; // peak PSF signal to PSF and background induced shot noise
	int *d_Hist_PSFWidth; // PSF width in pixel
	

	// get average SNR and PSF width
	float *h_ValidNum;
	float *h_LocPrecSum;
	float *h_SNRSum;
	float *h_PSFWSum;
	float *h_BgSum;

	float *d_ValidNum;
	float *d_LocPrecSum;
	float *d_SNRSum;
	float *d_PSFWSum;
	float *d_BgSum;

	// overlap molecular calculate to calculate localization density

	int *h_MaxFrame; // curremt max frame id

	int *h_FilteredFluoNum;
	int *h_FluoNum_n0; // 0 neighbor
	int *h_FluoNum_n1; // 1 neighbor

	int *d_FilteredFluoNum;
	float *d_FilteredLocArry;

	int *d_FluoNum_n0; // 0 neighbor
	int *d_FluoNum_n1; // 1 neighbor


public:

	void Init();
	void Deinit();

	// reset all distributions and time variation data
	void ResetAllDat(cudaStream_t cstream); // should be called befor first use GetStatisticalInf


	void GetStatisticalInf(float *h_LocArry, LocalizationPara & LocPara, int FluoNum, cudaStream_t cstream);

	void UpdateOntimeRatio(float *i_h_OntimeRatio);


public:
	static float GetActivationDensity(float Neibhgor0_Ratio, float Neibhgor1_Ratio, float RadiusTh_um);


	static float GetHistogramMeanData(int *HistData, int DatLen, float PercentTh);
	static int GetHistogramMaxData(int *HistData, int DatLen);
	static int GetHistogramMaxDataPos(int *HistData, int DatLen);
	static float GetHistogramWidth(int *HistData, int MaxPos, int DatLen);

	static float GetTimeVaryMean(vector<float>& iTimeVaryData);
	static float GetTimeVaryMax(vector<float>& iTimeVaryData);
	static float GetTimeVaryMin(vector<float>& iTimeVaryData);


private:
	// reset all distributions to avoid accumulate, keep time variation data,
	void ResetDistribDat(cudaStream_t cstream); // already be called in GetStatisticalInf()

	void GetOverlapMolecules(int*oNeighbor0_Num, int *oNeighbor1_Num, int *oTotalFluo, float RadiusTh_pixel, int CurFrame, LocalizationPara & LocPara, int FluoNum, cudaStream_t cstream);

	void UpdateStatDat(float *h_LocArry, int FluoNum);

	void CalcLocalizationDensity2D(float *h_LocArry, LocalizationPara & LocPara, int FluoNum, cudaStream_t cstream);

};



////

__global__ void CalcLocResultStatInf(float *d_LocArry, int *d_Hist_TotalPhoton, int *d_Hist_LocPrecisionXY, int *d_Hist_PeakSNR, int *d_Hist_PSFWidth, float *d_ValidNum, float *d_LocPrecSum, float *d_SNRSum, float *d_PSFWSum, float *d_BgSum, int FluoNum);

__global__ void FillerMoleculesInCurFrame(float *d_LocArry, float*d_FilteredLocArry, int *d_FilteredFluoNum, int FilterFrame, int LocType, int ImageWidth, int ImageHigh, int FluoNum);
__global__ void CalcOverlapFluoNum(float* d_LocArry, int *d_FluoNum_n0, int *d_FluoNum_n1, float RadiusTh_pixel, int FluoNum);





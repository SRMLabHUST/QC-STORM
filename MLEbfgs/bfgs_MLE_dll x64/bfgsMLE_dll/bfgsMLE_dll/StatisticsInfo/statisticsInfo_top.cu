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


#include "CurveFitting.h"
#include "cudaWrapper.h"

#include "statisticsInfo.h"


void FluoStatisticData_TypeDef::UpdateOntimeRatio(float *i_h_OntimeRatio)
{
	int cnt;

	memcpy(h_OntimeRatio, i_h_OntimeRatio, MaxOnTimeConsecutiveNum * sizeof(float));


	BFGSOptimizer< float, 2, 4, 11> ExpCurveFit(ExpFit_PreFitting, ExpFit_TargerF);
	ExpCurveFit.BFGSOptimize(NULL, h_OntimeRatio, MaxOnTimeConsecutiveNum);

	float k1 = -logf(1 - ExpCurveFit.FitPara[0]);
	float k2 = -ExpCurveFit.FitPara[1];
	MeanOntime = 2 / (k1 + k2); // frame



	// for cut off length data display, only show 4
	for (cnt = 4; cnt < MaxOnTimeConsecutiveNum; cnt++)
	{
		h_OntimeRatio[3] += h_OntimeRatio[cnt];
	}

	for (cnt = 0; cnt < 4; cnt++)
	{
		h_OntimeRatio[cnt] = max(h_OntimeRatio[cnt], 0.0f);
		h_OntimeRatio[cnt] = min(h_OntimeRatio[cnt], 1.0f);

	}

	// update  Variation  for display
	TimeVary_OntimeF1Ratio.push_back(h_OntimeRatio[0]);
	TimeVary_OntimeF2Ratio.push_back(h_OntimeRatio[1]);
	TimeVary_OntimeF3Ratio.push_back(h_OntimeRatio[2]);
	TimeVary_OntimeF4Ratio.push_back(h_OntimeRatio[3]);

}


void FluoStatisticData_TypeDef::GetStatisticalInf(float *h_LocArry, LocalizationPara & LocPara, int FluoNum, cudaStream_t cstream)
{

	if (FluoNum <= 0)return;


	// calculate stastics information
	cudaMemsetAsync(d_ValidNum, 0, sizeof(float), cstream);
	cudaMemsetAsync(d_LocPrecSum, 0, sizeof(float), cstream);
	cudaMemsetAsync(d_SNRSum, 0, sizeof(float), cstream);
	cudaMemsetAsync(d_PSFWSum, 0, sizeof(float), cstream);
	cudaMemsetAsync(d_BgSum, 0, sizeof(float), cstream);


	cudaMemcpyAsync(d_LocArry, h_LocArry, FluoNum*OutParaNumGS2D * sizeof(float), cudaMemcpyHostToDevice, cstream);

	int BlockDim = ThreadsPerBlock;
	int BlockNum = ((FluoNum + ThreadsPerBlock - 1) / ThreadsPerBlock);

	CalcLocResultStatInf << <BlockNum, BlockDim, 0, cstream >> >(d_LocArry, d_Hist_TotalPhoton, d_Hist_LocPrecisionXY, d_Hist_PeakSNR, d_Hist_PSFWidth, d_ValidNum, d_LocPrecSum, d_SNRSum, d_PSFWSum, d_BgSum, FluoNum);

	cudaMemcpyAsync(h_Hist_TotalPhoton, d_Hist_TotalPhoton, StatInf_Hist_DatLenMax * sizeof(int), cudaMemcpyDeviceToHost, cstream);
	cudaMemcpyAsync(h_Hist_LocPrecisionXY, d_Hist_LocPrecisionXY, StatInf_Hist_DatLenMax * sizeof(int), cudaMemcpyDeviceToHost, cstream);
	cudaMemcpyAsync(h_Hist_PeakSNR, d_Hist_PeakSNR, StatInf_Hist_DatLenMax * sizeof(int), cudaMemcpyDeviceToHost, cstream);
	cudaMemcpyAsync(h_Hist_PSFWidth, d_Hist_PSFWidth, StatInf_Hist_DatLenMax * sizeof(int), cudaMemcpyDeviceToHost, cstream);

	cudaMemcpyAsync(h_ValidNum, d_ValidNum, sizeof(float), cudaMemcpyDeviceToHost, cstream);
	cudaMemcpyAsync(h_LocPrecSum, d_LocPrecSum, sizeof(float), cudaMemcpyDeviceToHost, cstream);
	cudaMemcpyAsync(h_SNRSum, d_SNRSum, sizeof(float), cudaMemcpyDeviceToHost, cstream);
	cudaMemcpyAsync(h_PSFWSum, d_PSFWSum, sizeof(float), cudaMemcpyDeviceToHost, cstream);
	cudaMemcpyAsync(h_BgSum, d_BgSum, sizeof(float), cudaMemcpyDeviceToHost, cstream);

	cudaStreamSynchronize(cstream);

	if ((*h_ValidNum) <= 0)*h_ValidNum = 1;

	// avoid histgram accumulate
	ResetDistribDat(cstream);


	CalcLocalizationDensity2D(h_LocArry, LocPara, FluoNum, cstream);

	// update to display data
	UpdateStatDat(h_LocArry, FluoNum);

}

void FluoStatisticData_TypeDef::CalcLocalizationDensity2D(float *h_LocArry, LocalizationPara & LocPara, int FluoNum, cudaStream_t cstream)
{

	// calculate current activation density
	float RadiusTh_pixel;
	float RadiusTh_um;

	if (LocPara.LocType == LocType_GS2D)
	{
		RadiusTh_pixel = LocDensityCalc_RadiusTh_pixel_2D; // 2d localization
	}
	else
	{
		RadiusTh_pixel = LocDensityCalc_RadiusTh_pixel_3D; // 3d localization, commonly sparser than 2d
	}


	RadiusTh_um = RadiusTh_pixel*LocPara.PixelSize / 1000; // distance threshold in um


	int StartFrame = LDLocData_TypeDef::GetFirstFrame(h_LocArry, FluoNum);
	int EndFrame = LDLocData_TypeDef::GetLastFrame(h_LocArry, FluoNum);

	int TotalFrame = EndFrame - StartFrame + 1;

	int FilteredStartFrame = StartFrame;
	int FilteredEndFrame = EndFrame;

	if (TotalFrame > 3)
	{
		 FilteredStartFrame = StartFrame + 1;
		 FilteredEndFrame = EndFrame - 1;
	}


	float Neibhgor0_Ratio = 0;
	float Neibhgor1_Ratio = 0;
	int TotalNum = 0;
	int CurFillerFrame = 0;



	// get molular loc arry only in this frame
	// accumulate several images to avlid small images, single image can also be used
	for (CurFillerFrame = FilteredStartFrame; CurFillerFrame <= FilteredEndFrame; CurFillerFrame++)
	{
		int oNeighbor0_Num, oNeighbor1_Num;
		int oTotalFluo;


		GetOverlapMolecules(&oNeighbor0_Num, &oNeighbor1_Num, &oTotalFluo, RadiusTh_pixel, CurFillerFrame, LocPara, FluoNum, cstream);

		Neibhgor0_Ratio += oNeighbor0_Num;
		Neibhgor1_Ratio += oNeighbor1_Num;

		TotalNum += oTotalFluo;
	}

	if (TotalNum > 50)
	{
		// if there are enough molecular
		Neibhgor0_Ratio /= TotalNum;
		Neibhgor1_Ratio /= TotalNum;

		MeanLocDensity2D = GetActivationDensity(Neibhgor0_Ratio, Neibhgor1_Ratio, RadiusTh_um);

		// a compensation between practical and simulation

		if (LocPara.MultiEmitterFitEn)
		{
			// second order polynomial compensation
			MeanLocDensity2D = MeanLocDensity2D*MeanLocDensity2D* 0.5098f + MeanLocDensity2D* 0.6364f; 
		}
		else
		{
			MeanLocDensity2D = MeanLocDensity2D*MeanLocDensity2D* 0.8082f + MeanLocDensity2D* 0.519f;

		}

	}


}


void FluoStatisticData_TypeDef::UpdateStatDat(float *h_LocArry, int FluoNum)
{

	// get distribution mean value by center of mass
	// only calculate pos larger than 5% of highest distribution 
	MeanTotalPhoton = Hist_PhotonGap* GetHistogramMeanData(h_Hist_TotalPhoton, StatInf_Hist_DatLenMax, 0.1f);


	MeanLocPrecisionXY = *h_LocPrecSum / *h_ValidNum;
	MeanPeakSNR = *h_SNRSum / *h_ValidNum;
	MeanBackground = *h_BgSum / *h_ValidNum;

	float MeanPSFWidth1 = (*h_PSFWSum) / (*h_ValidNum);

	MeanPSFWidth = MeanPSFWidth1;
	MeanPSFWidth_Ctl = MeanPSFWidth1 * (-MeanPeakSNR); // lower is better for feedback control


	// update  Variation  for display
	TimeVary_TotalPhoton.push_back(MeanTotalPhoton);
	TimeVary_LocPrecisionXY.push_back(MeanLocPrecisionXY);
	TimeVary_PeakSNR.push_back(MeanPeakSNR);
	TimeVary_Background.push_back(MeanBackground);

	TimeVary_PSFWidth.push_back(MeanPSFWidth);
	TimeVary_PSFWidth_Ctl.push_back(MeanPSFWidth_Ctl);

	TimeVary_LocDensity2D.push_back(MeanLocDensity2D);

}


void FluoStatisticData_TypeDef::GetOverlapMolecules(int *oNeighbor0_Num, int *oNeighbor1_Num, int *oTotalFluo, float RadiusTh_pixel, int CurFrame, LocalizationPara & LocPara, int FluoNum, cudaStream_t cstream)
{

	// analysis activation density based on overlap molular ratio in current frame
	cudaMemsetAsync(d_FilteredFluoNum, 0, sizeof(int), cstream);

	cudaMemsetAsync(d_FluoNum_n0, 0, sizeof(int), cstream);
	cudaMemsetAsync(d_FluoNum_n1, 0, sizeof(int), cstream);

	int BlockDim = ThreadsPerBlock;
	int BlockNum = ((FluoNum + ThreadsPerBlock - 1) / ThreadsPerBlock);

	FillerMoleculesInCurFrame << <BlockNum, BlockDim, 0, cstream >> >(d_LocArry, d_FilteredLocArry, d_FilteredFluoNum, CurFrame, LocPara.LocType, LocPara.ImageWidth, LocPara.ImageHigh, FluoNum);
	cudaMemcpyAsync(h_FilteredFluoNum, d_FilteredFluoNum, sizeof(int), cudaMemcpyDeviceToHost, cstream);

	cudaStreamSynchronize(cstream);

	//	printf("filter frame:%d %d %d \n", StartFrame, CurFillerFrame, *h_FilteredFluoNum);

	int TotalFluoNum = *h_FilteredFluoNum;

	BlockNum = ((TotalFluoNum + ThreadsPerBlock - 1) / ThreadsPerBlock);

	CalcOverlapFluoNum << <BlockNum, BlockDim, 0, cstream >> >(d_FilteredLocArry, d_FluoNum_n0, d_FluoNum_n1, RadiusTh_pixel, TotalFluoNum);

	cudaMemcpyAsync(h_FluoNum_n0, d_FluoNum_n0, sizeof(int), cudaMemcpyDeviceToHost, cstream);
	cudaMemcpyAsync(h_FluoNum_n1, d_FluoNum_n1, sizeof(int), cudaMemcpyDeviceToHost, cstream);

	cudaStreamSynchronize(cstream);


	*oNeighbor0_Num = (*h_FluoNum_n0);
	*oNeighbor1_Num = (*h_FluoNum_n1);

	*oTotalFluo = TotalFluoNum;
}


void FluoStatisticData_TypeDef::Init()
{

	// host and gpu
	cudaError_t err;

	// malloc gpu and cpu resource

	err = cudaMalloc((void **)&d_LocArry, MaxPointNum*OutParaNumGS2D * sizeof(float));

	// Distribution statistic info in gpu
	err = cudaMalloc((void **)&d_Hist_TotalPhoton, StatInf_Hist_DatLenMax * sizeof(int));
	err = cudaMalloc((void **)&d_Hist_LocPrecisionXY, StatInf_Hist_DatLenMax * sizeof(int));
	err = cudaMalloc((void **)&d_Hist_PeakSNR, StatInf_Hist_DatLenMax * sizeof(int));
	err = cudaMalloc((void **)&d_Hist_PSFWidth, StatInf_Hist_DatLenMax * sizeof(int));

	HandleErr(err, "cudaMalloc d_Hist_PSFWidth");

	// Distribution statistic info in cpu
	err = cudaMallocHost((void **)&h_Hist_TotalPhoton, StatInf_Hist_DatLenMax * sizeof(int));
	err = cudaMallocHost((void **)&h_Hist_LocPrecisionXY, StatInf_Hist_DatLenMax * sizeof(int));
	err = cudaMallocHost((void **)&h_Hist_PeakSNR, StatInf_Hist_DatLenMax * sizeof(int));
	err = cudaMallocHost((void **)&h_Hist_PSFWidth, StatInf_Hist_DatLenMax * sizeof(int));


	// on time variation
	h_OntimeRatio = new float[MaxOnTimeConsecutiveNum];


	// overlap molecular calculate

	err = cudaMallocHost((void **)&h_MaxFrame, sizeof(int));

	err = cudaMallocHost((void **)&h_FilteredFluoNum, sizeof(int));
	err = cudaMallocHost((void **)&h_FluoNum_n0, sizeof(int));
	err = cudaMallocHost((void **)&h_FluoNum_n1, sizeof(int));

	err = cudaMalloc((void **)&d_FilteredFluoNum, sizeof(int));
	err = cudaMalloc((void **)&d_FilteredLocArry, MaxPointNum*OutParaNumGS2D * sizeof(float));

	err = cudaMalloc((void **)&d_FluoNum_n0, sizeof(int));
	err = cudaMalloc((void **)&d_FluoNum_n1, sizeof(int));

	// get average SNR and PSF width
	err = cudaMallocHost((void **)&h_ValidNum, sizeof(float));
	err = cudaMallocHost((void **)&h_LocPrecSum, sizeof(float));
	err = cudaMallocHost((void **)&h_SNRSum, sizeof(float));
	err = cudaMallocHost((void **)&h_PSFWSum, sizeof(float));
	err = cudaMallocHost((void **)&h_BgSum, sizeof(float));

	err = cudaMalloc((void **)&d_ValidNum, sizeof(float));
	err = cudaMalloc((void **)&d_LocPrecSum, sizeof(float));
	err = cudaMalloc((void **)&d_SNRSum, sizeof(float));
	err = cudaMalloc((void **)&d_PSFWSum, sizeof(float));
	err = cudaMalloc((void **)&d_BgSum, sizeof(float));



}

void FluoStatisticData_TypeDef::Deinit()
{
	cudaError_t err;

	cudaFree(d_LocArry);

	err = cudaFree(d_Hist_TotalPhoton);
	err = cudaFree(d_Hist_LocPrecisionXY);
	err = cudaFree(d_Hist_PeakSNR);
	err = cudaFree(d_Hist_PSFWidth);

	HandleErr(err, "cudaFree d_Hist_PSFWidth");

	//
	cudaFreeHost(h_Hist_TotalPhoton);
	cudaFreeHost(h_Hist_LocPrecisionXY);
	cudaFreeHost(h_Hist_PeakSNR);
	cudaFreeHost(h_Hist_PSFWidth);



	// on time variation
	delete[] h_OntimeRatio;


	// overlap molecular calculate
	cudaFreeHost(h_MaxFrame);

	cudaFreeHost(h_FilteredFluoNum);
	cudaFreeHost(h_FluoNum_n0);
	cudaFreeHost(h_FluoNum_n1);

	cudaFree(d_FilteredFluoNum);
	cudaFree(d_FilteredLocArry);

	cudaFree(d_FluoNum_n0);
	cudaFree(d_FluoNum_n1);


	// get average SNR and PSF width
	cudaFreeHost(h_ValidNum);
	cudaFreeHost(h_SNRSum);
	cudaFreeHost(h_PSFWSum);
	cudaFreeHost(h_BgSum);

	cudaFree(d_ValidNum);
	cudaFree(d_SNRSum);
	cudaFree(d_PSFWSum);
	cudaFree(d_BgSum);


}


void FluoStatisticData_TypeDef::ResetDistribDat(cudaStream_t cstream)
{

	cudaMemsetAsync(d_Hist_TotalPhoton, 0, StatInf_Hist_DatLenMax * sizeof(int), cstream);
	cudaMemsetAsync(d_Hist_LocPrecisionXY, 0, StatInf_Hist_DatLenMax * sizeof(int), cstream);
	cudaMemsetAsync(d_Hist_PeakSNR, 0, StatInf_Hist_DatLenMax * sizeof(int), cstream);
	cudaMemsetAsync(d_Hist_PSFWidth, 0, StatInf_Hist_DatLenMax * sizeof(int), cstream);

	cudaStreamSynchronize(cstream);

}



void FluoStatisticData_TypeDef::ResetAllDat(cudaStream_t cstream)
{
	ResetDistribDat(cstream);

	// time variation curve data 
	TimeVary_TotalPhoton.clear();
	TimeVary_LocPrecisionXY.clear();
	TimeVary_PeakSNR.clear();
	TimeVary_Background.clear();
	TimeVary_PSFWidth.clear();
	TimeVary_PSFWidth_Ctl.clear();
	TimeVary_LocDensity2D.clear();


	TimeVary_OntimeF1Ratio.clear();
	TimeVary_OntimeF2Ratio.clear();
	TimeVary_OntimeF3Ratio.clear();
	TimeVary_OntimeF4Ratio.clear();


	// mean value of current data
	MeanTotalPhoton = 0;
	MeanLocPrecisionXY = 0;
	MeanPeakSNR = 0;
	MeanBackground = 0;
	MeanPSFWidth = 0;
	MeanPSFWidth_Ctl = 0;
	MeanLocDensity2D = 0;


}


float FluoStatisticData_TypeDef::GetTimeVaryMean(vector<float> &iTimeVaryData)
{
	float MeanData = 0;

	if (iTimeVaryData.size() > 2)
	{
		float TotalData = 0;
		int ValidNum = 0;

		for (int cnt = 1; cnt < iTimeVaryData.size() - 1; cnt++)
		{
			TotalData += iTimeVaryData[cnt];
			ValidNum++;
		}

		if (ValidNum <= 0)ValidNum = 1;

		MeanData = TotalData / ValidNum;
	}

	return MeanData;
}

float FluoStatisticData_TypeDef::GetTimeVaryMax(vector<float> &iTimeVaryData)
{
	float MaxDat = 0;

	if (iTimeVaryData.size() > 0)
	{
		MaxDat = iTimeVaryData[0];

		for (int cnt = 0; cnt < iTimeVaryData.size(); cnt++)
		{
			MaxDat = max(MaxDat, iTimeVaryData[cnt]);
		}
	}

	return MaxDat;
}

float FluoStatisticData_TypeDef::GetTimeVaryMin(vector<float> &iTimeVaryData)
{
	float MinDat = 0;

	if (iTimeVaryData.size() > 0)
	{
		 MinDat = iTimeVaryData[0];

		for (int cnt = 0; cnt < iTimeVaryData.size(); cnt++)
		{
			MinDat = min(MinDat, iTimeVaryData[cnt]);
		}
	}

	return MinDat;
}


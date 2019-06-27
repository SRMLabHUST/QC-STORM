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

#include "bfgs_MLE_dll.h"

extern LocalizationPara LocPara_Global;


extern int GPUID_1Best;
extern int GPUID_2Best;


// molecule finding, MLE localization,rendering,statistical processing
class ZDriftMeasureData_TypeDef
{

public:
	LDROIExtractData_TypeDef LDROIExtractData_;
	LDLocData_TypeDef LDLocData_;

	ZeroLocalizationsRemovel_TypeDef ZeroLocRemovel_;

	FluoStatisticData_TypeDef FluoStatData_;

	cudaStream_t loc_stream1_;


	int MLELocalization(unsigned short *h_RawImg, int BatchedImgNum)
	{
		cudaSetDevice(GPUID_1Best);

		FluoStatData_.ResetAllDat(loc_stream1_);



		// subregion extraction for current image
		LDROIExtractData_.ExtractMolecules(LDROIExtractData_.h_RawImg, ImageSource_CPU_Pinned, LocPara_Global, 1, BatchedImgNum, loc_stream1_);

		int FluoNum = LDROIExtractData_.GetAccumulatedROINum();
		LDROIExtractData_.ResetROINum();


		// localization
		LDLocData_.BFGS_MLELocalization(LDROIExtractData_.h_ImageROI, LDROIExtractData_.Get_h_WLEPara(), LocPara_Global, FluoNum, loc_stream1_);

		// remove invalid molecules and sort frame, frame is disordered by LDROIExtractData and LDLocData
		ZeroLocRemovel_.RemoveZeroLocalizations(LDLocData_.h_LocArry, LDLocData_.oValidFluoNum, 1, 1, BatchedImgNum, loc_stream1_);


		// get statistic information
		FluoStatData_.GetStatisticalInf(ZeroLocRemovel_.h_LocArry, LocPara_Global, ZeroLocRemovel_.ValidFluoNum, loc_stream1_);

		FluoNum = ZeroLocRemovel_.ValidFluoNum;

//		printf("mean psf:%f %f\n", FluoStatData_.h_MeanPSFWidth, FluoStatData_.h_MeanPSFWidth_ctl);

		return FluoNum;
	}

	unsigned short* GetRawImageMem()
	{
		return LDROIExtractData_.h_RawImg;
	}

	void Init(LocalizationPara &LocPara_Global)
	{
		cudaSetDevice(GPUID_1Best);

		LDROIExtractData_.Init(LocPara_Global);
		LDLocData_.Init(LocPara_Global);

		ZeroLocRemovel_.Init();

		FluoStatData_.Init();

		CreatStream(&loc_stream1_);

		int BatchedImgNum = 32;
		int BatchedImgSize = BatchedImgNum*LocPara_Global.ImageWidth*LocPara_Global.ImageHigh;

	}

	void Deinit(LocalizationPara &LocPara_Global)
	{
		cudaSetDevice(GPUID_1Best);

		LDROIExtractData_.Deinit();
		LDLocData_.Deinit(LocPara_Global);

		ZeroLocRemovel_.Deinit();

		FluoStatData_.Deinit();

		// free stream
		FreeStream(loc_stream1_);

	}

};



extern ZDriftMeasureData_TypeDef ZDriftCtl;



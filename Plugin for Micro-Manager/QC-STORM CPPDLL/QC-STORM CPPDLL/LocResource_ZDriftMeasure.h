#pragma once


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

#include "bfgs_MLE_dll.h"

#include "LocResource.h"


// molecule finding, MLE localization,rendering,statistical processing
class ZDriftMeasureData_TypeDef
{

public:
	LDROIExtractData_TypeDef LDROIExtractData_;
	LDLocData_TypeDef LDLocData_;
	FluoStatisticData_TypeDef FluoStatData_;

	cudaStream_t loc_stream1_;



	int MLELocalization(unsigned short *h_RawImg, int BatchedImgNum)
	{

		FluoStatData_.ResetAllDat(loc_stream1_);


		memcpy(LDROIExtractData_.h_RawImg, h_RawImg, BatchedImgNum*LocPara_Global.ImageWidth*LocPara_Global.ImageHigh * sizeof(short));

		// subregion extraction for current image
		LDROIExtractData_.ExtractMolecules(LDROIExtractData_.h_RawImg, ImageSource_CPU_Pinned, LocPara_Global, 1, BatchedImgNum, loc_stream1_);

		int FluoNum = LDROIExtractData_.GetAccumulatedROINum();
		LDROIExtractData_.ResetROINum();


		// localization
		LDLocData_.BFGS_MLELocalization(LDROIExtractData_.h_ImageROI, LDROIExtractData.Get_h_WLEPara(), LocPara_Global, FluoNum, loc_stream1_);

		// get statistic information
		FluoStatData_.GetStatisticalInf(LDLocData_.h_LocArry, LocPara_Global, LDLocData_.oValidFluoNum, loc_stream1_);

		FluoNum = LDLocData_.oValidFluoNum;

//		printf("mean psf:%f %f\n", FluoStatData_.h_MeanPSFWidth, FluoStatData_.h_MeanPSFWidth_ctl);

		printf("batch loc molecule num:%d\n", FluoNum);


		return FluoNum;
	}


	void Init(LocalizationPara &LocPara_Global)
	{
		LDROIExtractData_.Init(LocPara_Global);
		LDLocData_.Init(LocPara_Global);

		FluoStatData_.Init();

		CreatStream(&loc_stream1_);
	}

	void Deinit(LocalizationPara &LocPara_Global)
	{
		LDROIExtractData_.Deinit();
		LDLocData_.Deinit(LocPara_Global);

		FluoStatData_.Deinit();

		// free stream
		FreeStream(loc_stream1_);
	}

};



extern ZDriftMeasureData_TypeDef ZDriftCtl;



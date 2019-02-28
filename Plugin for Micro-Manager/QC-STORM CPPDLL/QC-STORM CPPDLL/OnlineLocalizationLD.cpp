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

#include "stdafx.h"
#include "OnlineLocalizationLD.h"

#include <time.h>

#include "OnlineFeedback_rc.h"




#include <io.h>  
#include <fcntl.h>  
void OpenConsole()
{
	// create a console
	AllocConsole();

	// 
	FILE* stream;
	freopen_s(&stream, "CON", "r", stdin);//
	freopen_s(&stream, "CON", "w", stdout);//

	SetConsoleTitleA("Information output");//

	HANDLE _handleOutput;
	_handleOutput = GetStdHandle(STD_OUTPUT_HANDLE);

	// FreeConsole();
}



UINT th_OnlineLocalizationLD(LPVOID params)
{
	IsLocRunning = true;

	//	OpenConsole();
	printf("image inf:%d %d %d %d \n", LocPara_Global.ImageWidth, LocPara_Global.ImageHigh, LocPara_Global.ROISize, LocPara_Global.LocType);
	printf("image inf:%f %f %f \n", LocPara_Global.Offset, LocPara_Global.KAdc, LocPara_Global.QE);
	printf("structure size:%f\n", LocPara_Global.StrucuteSize_2D);


	int FluoNum = 0;
	int CurFrame = 0;
	int BatchFrameNum = 0;
	int LastFrame = 0;



	// accumulate a group of molecular from multiple images to accelerate localization
	int RecFluoNumTh = PointNumTh/2;


	if (LocPara_Global.SpatialResolutionCalcEn)
	{
		RecFluoNumTh = PointNumTh*(LocPara_Global.ImageWidth + LocPara_Global.ImageHigh) / 2 / 2048;
		RecFluoNumTh = min(RecFluoNumTh, PointNumTh);
		RecFluoNumTh = max(RecFluoNumTh, 4096);
		RecFluoNumTh = RecFluoNumTh / 32 * 32;
	}



	CString LocFileName;
	LocFileName.Format(_T("loc_result%dD%d_%s.txt"), LocPara_Global.LocType + 2, LocPara_Global.ROISize, CreateTimeIdxStr);
	LocFileName = ResultSavePath + LocFileName;


	// write to file
	float *WriteLocArry = NULL;
	int WriteLocNum = 0;

	CFile LocFile;

	bool suc = LocFile.Open(LocFileName, CFile::modeCreate | CFile::modeWrite); //CFile::modeNoTruncate |
	if (!suc)
	{
		printf("invalid result save path\n");
	}



	// time counting
	unsigned int time1;
	unsigned int WholeProc_StartTime, WholeProc_EndTime;

	unsigned int ExtractTime = 0;
	unsigned int LocTime = 0;
	unsigned int RendTime = 0;
	unsigned int StatTime = 0;
	unsigned int OnTimeCalcTime = 0;
	unsigned int ResolutionTime = 0;

	unsigned int TotalFluoNum = 0;

	//
	CurDispPosX = LocPara_Global.ImageWidth / 2;
	CurDispPosY = LocPara_Global.ImageHigh / 2;
	CurDispZoom = min(1.0f *SR_IMAGE_DISPLAY_WIDTH / (LocPara_Global.SRImageWidth), 1.0f *SR_IMAGE_DISPLAY_HIGH / (LocPara_Global.SRImageHigh));



	CurFrame = 0;
	LastFrame = 0;

	bool IsBreak = 0;


	// reset rendered image
	RendData.ResetFillImgTop(LocPara_Global);
	FluoStatData.ResetAllDat(loc_stream1);

	// spatial resolution calculation
	DimensionDensityCalc.ResetAccumulatedData();



	SpatialResolutionCalc.ResetData();
	SpatialResolutionCalc.SetStructureSize(LocPara_Global.StrucuteSize_2D);


	RenderingState.MakeAProcess(); // display window

	qImgData RecImg;

	WholeProc_StartTime = clock();

	while (1)
	{
		
		IsBreak = !(OnlineLocAlive || (!ImgDataQueue.empty()));

		// get images and perform molecular detection
		if (ImgDataQueue.try_pop(RecImg))
		{
			// group many images to process them at the same time to inprove performance
			BatchFrameNum = RecImg.BatchFrameNum;


			unsigned short * ImageMemoryBuf = NULL;
			int ImageMemorySource = ImageSource_ERR;

			//get image from queue
			if (RecImg.ImageSource == ImageSource_GPU)
			{
				// GPU type image memory
				ImageMemoryBuf = RecImg.pImgData;
				ImageMemorySource = ImageSource_GPU;
			}
			else if (RecImg.ImageSource == ImageSource_CPU_Normal)
			{
				// LDROIExtractData.h_RawImg is pinned
				ImageMemoryBuf = LDROIExtractData.h_RawImg;
				ImageMemorySource = ImageSource_CPU_Pinned;

				// CPU type image memory
				memcpy(ImageMemoryBuf, RecImg.pImgData, BatchFrameNum*LocPara_Global.ImageWidth*LocPara_Global.ImageHigh * 2);
			}


			time1 = clock();
			// subregion extraction for current image
			LDROIExtractData.ExtractMolecules(ImageMemoryBuf, ImageMemorySource, LocPara_Global, CurFrame + 1, BatchFrameNum, loc_stream1);

			ExtractTime += (clock() - time1);


			CurFrame += BatchFrameNum;


			// after use the image
			if (RecImg.ImageSource == ImageSource_CPU_Normal)
			{
				delete[] RecImg.pImgData;
			}
			else if (RecImg.ImageSource == ImageSource_GPU)
			{
				cudaFree(RecImg.pImgData);
			}


		}

		// accumulate enough molecule for fast localization
		if ((LDROIExtractData.GetRegionNum() >= RecFluoNumTh) || IsBreak)
		{
			FluoNum = LDROIExtractData.GetRegionNum();
			LDROIExtractData.ResetRegionNum();

			TotalFluoNum += FluoNum;


			time1 = clock();

			// localization
			LDLocData.BFGS_MLELocalization(LDROIExtractData.h_RegionMem, LocPara_Global, FluoNum, loc_stream1);

			LocTime += (clock() - time1);
			

			// write localization data into file
			WriteLocArry = LDLocData.h_LocArry;
			WriteLocNum = LDLocData.oValidFluoNum;



			time1 = clock();

			// consecutive fit
			if (LocPara_Global.ConsecFitEn)
			{
				ConsecutiveFitData.FitConsecutiveFluo(LDLocData.d_LocArry, LocPara_Global, LDLocData.oValidFluoNum, loc_stream1, IsBreak);

				if (!IsBreak)
				{
					// get avalible data
					ConsecutiveFitData.GetAvailableData(loc_stream1);
				}
				else
				{
					// use at the last time
					ConsecutiveFitData.GetResidualData(loc_stream1);
				}

				WriteLocArry = ConsecutiveFitData.h_OutLocArry;
				WriteLocNum = ConsecutiveFitData.OutFluoNum;
			}

			// remove invalid molecules
			ZeroLocRemovel.RemoveZeroLocalizations(WriteLocArry, WriteLocNum, loc_stream1);

			WriteLocArry = ZeroLocRemovel.h_LocArry;
			WriteLocNum = ZeroLocRemovel.ValidFluoNum;

			LocTime += (clock() - time1);


			// write localization data into file
			LocFile.Write(WriteLocArry, WriteLocNum * OutParaNumGS2D * sizeof(float));

			time1 = clock();
			// super resolution image rendering
			RendData.FluoRenderTop(WriteLocArry, LocPara_Global, RenderingMode_FittedPhoton_CalculatedLocPrec, 0, WriteLocNum, loc_stream1);

			RendTime += (clock() - time1);


			time1 = clock();

			if (LocPara_Global.OnTimeCalcEn)
			{
				// ontime stastics, should be call after LDLocData.BFGS_MLELocalization()
				LDLocData.OntimeCalc(LocPara_Global, FluoNum, loc_stream1);
				FluoStatData.UpdateOntimeRatio(LDLocData.h_OntimeRatio);
			}

			OnTimeCalcTime += (clock() - time1);

			time1 = clock();

			// get statistic information
			FluoStatData.GetStatisticalInf(WriteLocArry, LocPara_Global, WriteLocNum, loc_stream1);

			StatTime += (clock() - time1);


			if (!IsBreak)
			{
				// must after stastical data are avalible
				OnlineFeedbackState.MakeAProcess(); // update feedback
			}
			

			time1 = clock();

			if (LocPara_Global.SpatialResolutionCalcEn)
			{
				// filter margin data, to avoid overlap effect in multi ROI acquisition
				MarginDataFilter.FilterMarginData(WriteLocArry, WriteLocNum, (FOVOverlapPercent + 0.2f) / 100, LocPara_Global.ImageWidth, LocPara_Global.ImageHigh, loc_stream1);


				int IsEnough = DimensionDensityCalc.AddLocArray_FewFrames(WriteLocArry, WriteLocNum, LocPara_Global.ConsecFit_DistanceTh_nm / LocPara_Global.PixelSize, IsBreak, loc_stream1);

				if (IsEnough)
				{
					DimensionDensityCalc.GetDimensionLocDensity_AGroup(DimensionDensityCalc.ImagesPerGroup_Valid * 1 / 10, DimensionDensityCalc.ImagesPerGroup_Valid * 2 / 10, DimensionDensityCalc.ImagesPerGroup_Valid, LocPara_Global.PixelSize, RendType_Is3D(LocPara_Global.LocType), loc_stream1);

					// add current group's dimension and density to calculate spatial resolution
					SpatialResolutionCalc.AddDimensionLocDensity_AGroup(DimensionDensityCalc.DimensionFD, DimensionDensityCalc.LocDensityFD, DimensionDensityCalc.ImagesPerGroup_Valid);

					// get spatial resolution vary data
					bool Is3DImaging = LocPara_Global.LocType == LocType_AS3D;

					float Mean_LocPrecisionXY = FluoStatisticData_TypeDef::GetTimeVaryMean(FluoStatData.TimeVary_LocPrecisionXY);
					float Mean_LocPrecisionZ = Mean_LocPrecisionXY / 1.414f * 2.0f;

					SpatialResolutionCalc.GetSpatialResolutionVary(Is3DImaging, LocPara_Global.IsHollowTube, Mean_LocPrecisionXY, Mean_LocPrecisionZ, NYQUIST_RESOLUTION_OVERSAMPLING);


					DimensionDensityCalc.ResetAccumulatedData();

				}
			}

			ResolutionTime += (clock() - time1);



			// at least accumulate 10 frame to refresh display
			// it's not madatory for accumulate several images, but for high efficient and not update too frequently
			if (CurFrame - LastFrame >= 20)
			{
				RenderingState.MakeAProcess(); // update super-resolution image display and stastical information
				LastFrame = CurFrame;
			}
		}

		if (IsBreak)break;
	}


	RenderingState.MakeAProcess();
	LocFile.Close();

	while (RenderingState.HaveProcessWait()); // wait rend finish

	RendSaveImage();


	ClearAllOnlineThreadAlive();

	WholeProc_EndTime = clock();


	printf("out loc th\n");



	if (LocTime <= 0)LocTime = 1;
	float LocSpeed = TotalFluoNum*1.0f / LocTime*1000.0f;

	printf("total time : %d ms\n", WholeProc_EndTime - WholeProc_StartTime);

	printf("Molecular Extraction time: %d ms\n", ExtractTime);
	printf("Localization time : %d ms\n", LocTime);
	printf("Rendering time : %d ms\n", RendTime);
	printf("Stastics calc time : %d ms\n", StatTime);
	printf("OnTime calc time : %d ms\n", OnTimeCalcTime);

	printf("Resolution calc time : %d ms\n", ResolutionTime);

	printf("TotalFluoNum : %d\n", TotalFluoNum);
	printf("Localization speed : %.0f/s\n", LocSpeed);


	IsLocRunning = false;


	return 0;
}

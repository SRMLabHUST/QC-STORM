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
#include "PostProcess.h"

#include "CppFileWrapper.h"

CString RerendDataPath;

CWinThread *Thread_Rerend = NULL;
CWinThread *Thread_ConvertBin = NULL;


volatile int RerendProgress = 0;

int IsDriftCorrection = 0;
int DriftCorrGroupFrameNum = 800;

UINT th_RerendImage(LPVOID params)
{
	IsLocRunning = true;

	OpenConsole(); // open console window to display printf

	RerendProgress = 0; // start

	int TotalFluoNum;
	int CurFluoNum = PointNumTh;
	int cnt;

	int ImageWidth = 0;
	int ImageHigh = 0;


	int TotalFrame = 0;
	

	CString CorrectedDataPath = RerendDataPath;

	if (IsDriftCorrection)
	{
		CString PostFixStr;
		PostFixStr.Format(L"_DriftCorrected_g%df.txt", DriftCorrGroupFrameNum);

		CorrectedDataPath.TrimRight(L".txt");
		CorrectedDataPath += PostFixStr;
	}


	// find image size from loc arry first
	
	wstring iFileName_ws = RerendDataPath.GetBuffer();
	wstring oFileName_ws = CorrectedDataPath.GetBuffer();

	string iFileName = EncodeConvert::ws2s(iFileName_ws);
	string oFileName = EncodeConvert::ws2s(oFileName_ws);



	// get image size
	SRDriftCorrData_TypeDef::GetMaxImgSize(iFileName, &ImageWidth, &ImageHigh);
	TotalFrame = SRDriftCorrData_TypeDef::GetTotalFrame(iFileName);

	printf("find img size:%d %d %d\n", ImageWidth, ImageHigh, TotalFrame);

	LocPara_Global.ImageWidth = ImageWidth;
	LocPara_Global.ImageHigh = ImageHigh;
	LocPara_Global.UpdateSRImageSize();




	if (IsDriftCorrection)
	{
		printf("begin shift correction\n");

		cudaStream_t CurStream1;
		cudaStreamCreate(&CurStream1);

		SRDriftCorrData_TypeDef h_SRDriftCorrData;
		h_SRDriftCorrData.Init(ImageWidth, ImageHigh, TotalFrame);

		// 500 frame as a rendering group for cross-correlation


		int time1, time2;

		time1 = clock();


		h_SRDriftCorrData.CorrectSampleShift(iFileName, oFileName, LocPara_Global, DriftCorrGroupFrameNum, CurStream1);


		time2 = clock();

		printf("cor time:%d\n", time2 - time1);

		h_SRDriftCorrData.Deinit();
		cudaStreamDestroy(CurStream1);

		printf("shift correction finish\n");
	}

	// begin rendering super-resolution image

	// rend image in progress


	InitAllLocResource(1);

	RendData.ResetFillImgTop(LocPara_Global);
	FluoStatData.ResetAllDat(loc_stream1);

	SpatialResolutionCalc.ResetData();
	SpatialResolutionCalc.SetStructureSize(LocPara_Global.StrucuteSize_2D);


	RerendProgress = 1; // get img size and super-resolution image resource are allocated
	
	CFile LocFile;
	LocFile.Open(CorrectedDataPath, CFile::modeRead);



	TotalFluoNum = LocFile.GetLength() / OutParaNumGS2D / sizeof(float);
	printf("total fluo:%d\n", TotalFluoNum);


	int RecFluoNumTh = PointNumTh;

	if (LocPara_Global.SpatialResolutionCalcEn)
	{
		// avoid long processing time,which result in error
		RecFluoNumTh = PointNumTh*(LocPara_Global.ImageWidth + LocPara_Global.ImageHigh) / 2 / 2048;
		RecFluoNumTh = min(RecFluoNumTh, PointNumTh);
		RecFluoNumTh = max(RecFluoNumTh, 4096);
		RecFluoNumTh = RecFluoNumTh / 32 * 32;
	}


	int ProcNum = (TotalFluoNum + RecFluoNumTh - 1) / RecFluoNumTh;




	for (cnt = 0; cnt < ProcNum; cnt++)
	{
		if (cnt == ProcNum - 1)CurFluoNum = TotalFluoNum - cnt*RecFluoNumTh;
		else CurFluoNum = RecFluoNumTh;

		bool IsBreak = (cnt == ProcNum - 1);

		LocFile.Read(LDLocData.h_LocArry, CurFluoNum*OutParaNumGS2D*sizeof(float));


		// write to file
		float *WriteLocArry = LDLocData.h_LocArry;
		int WriteLocNum = CurFluoNum;



		// super resolution image rendering
		RendData.FluoRenderTop(WriteLocArry, LocPara_Global, 0, 0, WriteLocNum, loc_stream1);



		// get statistic information
		FluoStatData.GetStatisticalInf(WriteLocArry, LocPara_Global, WriteLocNum, loc_stream1);


		if (LocPara_Global.SpatialResolutionCalcEn)
		{
			// dimension and density calculation, only calculate after enough frames are accumulated

			int IsEnough = DimensionDensityCalc.AddLocArray_FewFrames(WriteLocArry, WriteLocNum, LocPara_Global.ConsecFit_DistanceTh_nm / LocPara_Global.PixelSize, IsBreak, loc_stream1);


			if (IsEnough)
			{
				bool Is3DImaging = LocPara_Global.LocType == LocType_AS3D;

				DimensionDensityCalc.GetDimensionLocDensity_AGroup(DimensionDensityCalc.ImagesPerGroup_Valid * 1 / 10, DimensionDensityCalc.ImagesPerGroup_Valid * 2 / 10, DimensionDensityCalc.ImagesPerGroup_Valid - 1, LocPara_Global.PixelSize, Is3DImaging, loc_stream1);

//				printf("get dimension:%d %.2f %.2f %.2f\n", DimensionDensityCalc.ImagesPerGroup_Valid, DimensionDensityCalc.DimensionFD, DimensionDensityCalc.LocDensityFD, LocPara_Global.PixelSize);

				// add current group's dimension and density to calculate spatial resolution
				SpatialResolutionCalc.AddDimensionLocDensity_AGroup(DimensionDensityCalc.DimensionFD, DimensionDensityCalc.LocDensityFD, DimensionDensityCalc.ImagesPerGroup_Valid);

				// get spatial resolution vary data

				float Mean_LocPrecisionXY = FluoStatisticData_TypeDef::GetTimeVaryMean(FluoStatData.TimeVary_LocPrecisionXY);
				float Mean_LocPrecisionZ = Mean_LocPrecisionXY / 1.414f * 2.0f;

				SpatialResolutionCalc.GetSpatialResolutionVary(Is3DImaging, LocPara_Global.IsHollowTube, Mean_LocPrecisionXY, Mean_LocPrecisionZ, NYQUIST_RESOLUTION_OVERSAMPLING);


				DimensionDensityCalc.ResetAccumulatedData();
			}
	
		}

		if((cnt + 1) % 3 ==0)
		{
			RenderingState.MakeAProcess();
			while (RenderingState.HaveProcessWait()); // wait rend finish
		}
	}


	RenderingState.MakeAProcess();

	while (RenderingState.HaveProcessWait()); // wait rend finish
	OnlineRendAlive = false; // stop rend image

	LocFile.Close();




	printf("rend finish\n");

	RerendProgress = 2; // rend finish

	IsLocRunning = false;

	return 0;
}


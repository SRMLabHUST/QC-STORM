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

#include "OnlineSpatialResolutionCalc.h"

#include <time.h>

#include "OnlineFeedback_rc.h"

#include "DepthMapCalc.h"


#define DEPTH_MAP_DISP_EN		1


#if(DEPTH_MAP_DISP_EN)

CImgDisplay DepthMap_Display;

#endif


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
	cudaSetDevice(GPUID_1Best);

	IsLocRunning = true;

	//	OpenConsole();
	printf("image inf:%d %d %d %d \n", LocPara_Global.ImageWidth, LocPara_Global.ImageHigh, LocPara_Global.ROISize, LocPara_Global.LocType);
	printf("image inf:%f %f %f \n", LocPara_Global.Offset, LocPara_Global.KAdc, LocPara_Global.QE);
	printf("structure size:%f\n", LocPara_Global.StrucuteSize_2D);


	int CurFrame = 0;
	int BatchFrameNum = 0;
	int LastFrame = 0;



	// accumulate a group of molecular from multiple images to accelerate localization


	int RecFluoNumTh = PointNumTh*(LocPara_Global.ImageWidth + LocPara_Global.ImageHigh) / 2 / 1024;

	if (LocPara_Global.LocType == LocType_AS3D)
	{
		RecFluoNumTh = RecFluoNumTh / 2;
	}

	RecFluoNumTh = min(RecFluoNumTh, PointNumTh);
	RecFluoNumTh = max(RecFluoNumTh, 4096);
	RecFluoNumTh = RecFluoNumTh / 32 * 32;



	if (DepthMapDispEn)
	{
		RecFluoNumTh = 2048;
	}



	CString MultiFitStr;
	if (LocPara_Global.MultiEmitterFitEn) MultiFitStr = L"_M";
	else MultiFitStr = L"_S";

	CString ConsecFitStr;
	if (LocPara_Global.ConsecFitEn) ConsecFitStr.Format(L"_Consec%.0fnm", LocPara_Global.ConsecFit_DistanceTh_nm);
	else ConsecFitStr = L"";


	CString LocFileName;
	LocFileName.Format(_T("loc_result%dD%d_%s%s%s.txt"), LocPara_Global.LocType + 2, LocPara_Global.ROISize, CreateTimeIdxStr, MultiFitStr, ConsecFitStr);
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

	unsigned int TotalFluoNum = 0;

	//
	CurDispPosX = LocPara_Global.ImageWidth / 2;
	CurDispPosY = LocPara_Global.ImageHigh / 2;
	CurDispZoom = min(1.0f *SR_IMAGE_DISPLAY_WIDTH / (LocPara_Global.SRImageWidth), 1.0f *SR_IMAGE_DISPLAY_HIGH / (LocPara_Global.SRImageHigh));



	CurFrame = 0;
	LastFrame = 0;

	bool IsBreak = 0;

	// reset state
	LDROIExtractData.ResetROINum();

	RendData.ResetFillImgTop(LocPara_Global); // reset rendered image

	FluoStatData.ResetAllDat(loc_stream1);




	RenderingState.MakeAProcess(); // display window


	// receive image data send from plug-in GUI
	qImgData RecImg;

	// send localizations for resolution calc and multi-GPU acceleration
	qLocArray LocArray_Send;


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
			if (RecImg.ImageSource == ImageSource_CPU_Pinned)
			{
				// GPU type image memory
				ImageMemoryBuf = RecImg.pImgData;
				ImageMemorySource = ImageSource_CPU_Pinned;
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
			if (RecImg.ImageSource == ImageSource_CPU_Pinned)
			{
				cudaFreeHost(RecImg.pImgData);
			}
			if (RecImg.ImageSource == ImageSource_CPU_Normal)
			{
				delete[] RecImg.pImgData;
			}
		}

		// accumulate enough molecule for fast localization
		if ((LDROIExtractData.GetAccumulatedROINum() >= RecFluoNumTh) || IsBreak)
		{
			// get first and last frame of a batch
			int FirstFrame = LDROIExtractData.FirstFrame;
			int EndFrame = LDROIExtractData.EndFrame;

			int FluoNum = LDROIExtractData.GetAccumulatedROINum();
			LDROIExtractData.ResetROINum();



			time1 = clock();

			// localization
			LDLocData.BFGS_MLELocalization(LDROIExtractData.h_ImageROI, LDROIExtractData.Get_h_WLEPara(), LocPara_Global, FluoNum, loc_stream1);

			LocTime += (clock() - time1);
			
			TotalFluoNum += LDLocData.oValidFluoNum;



			// remove invalid molecules and sort frame, frame is disordered by LDROIExtractData and LDLocData
			ZeroLocRemovel.RemoveZeroLocalizations(LDLocData.h_LocArry, LDLocData.oValidFluoNum, 1, FirstFrame, EndFrame, loc_stream1);

			// write localization data into file
			WriteLocArry = ZeroLocRemovel.h_LocArry;
			WriteLocNum = ZeroLocRemovel.ValidFluoNum;


#if(DEPTH_MAP_DISP_EN)
			// only for beads localization
			if (DepthMapDispEn)
			{
				DepthMapCalc_TypeDef *DepthMapCalc = new DepthMapCalc_TypeDef(LocPara_Global.ImageWidth, LocPara_Global.ImageHigh);

				DepthMapCalc->FitZPlane(WriteLocArry, WriteLocNum);

				CImg<unsigned char> *DepthMapImage = new CImg<unsigned char>(DepthMapCalc->DepthMap_ImageWidth, DepthMapCalc->DepthMap_ImageHigh, 1, 3, 0);

				unsigned char *pImage = DepthMapImage->data();

				ConvertRGBToCImg(pImage, DepthMapCalc->DepthMap, DepthMapCalc->DepthMap_ImageWidth, DepthMapCalc->DepthMap_ImageHigh);


				DepthMap_Display.resize().display(*DepthMapImage);
				DepthMap_Display.set_title("depth map");
				if (DepthMap_Display.is_closed())
				{
					DepthMap_Display.show();
				}

				delete DepthMapCalc;
				delete DepthMapImage;

			}				

#endif


			time1 = clock();
			// consecutive fit
			if (LocPara_Global.ConsecFitEn)
			{

				ConsecutiveFitData.ConsecutiveFit_WeightedAvg(WriteLocArry, WriteLocNum, IsBreak, LocPara_Global, loc_stream1); // d_iLocArry come from localization data 

				// frame is not disordered by ConsecutiveFitData
				ZeroLocRemovel.RemoveZeroLocalizations(ConsecutiveFitData.h_OutLocArry, ConsecutiveFitData.OutFluoNum, 0, 0, 0, loc_stream1);

				// write localization data into file
				WriteLocArry = ZeroLocRemovel.h_LocArry;
				WriteLocNum = ZeroLocRemovel.ValidFluoNum;


				// ontime stastics
				FluoStatData.UpdateOntimeRatio(ConsecutiveFitData.h_OntimeRatio);

			}

			StatTime += (clock() - time1);


			// send for resolution calculation
			if (LocPara_Global.SpatialResolutionCalcEn)
			{

				LocArray_Send.h_LocArray = new float[WriteLocNum * OutParaNumGS2D];
				LocArray_Send.FluoNum = WriteLocNum;
				LocArray_Send.IsEndCalc = IsBreak;

				memcpy(LocArray_Send.h_LocArray, WriteLocArry, WriteLocNum * OutParaNumGS2D * sizeof(float));

				LocArray_Resolution_Queue.push(LocArray_Send);
			}


			// write localization data into file
			LocFile.Write(WriteLocArry, WriteLocNum * OutParaNumGS2D * sizeof(float));

			time1 = clock();
			// super resolution image rendering
			RendData.FluoRenderTop(WriteLocArry, LocPara_Global, RenderingMode_FittedPhoton_CalculatedLocPrec, 0, WriteLocNum, loc_stream1);

			RendTime += (clock() - time1);


			time1 = clock();

			// get statistic information
			FluoStatData.GetStatisticalInf(WriteLocArry, LocPara_Global, WriteLocNum, loc_stream1);

			StatTime += (clock() - time1);


			if (!IsBreak)
			{
				// must after stastical data are avalible
				OnlineFeedbackState.MakeAProcess(); // update feedback
			}



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


	if (LocPara_Global.SpatialResolutionCalcEn)
	{
		while (LocArray_Resolution_Queue.unsafe_size()>0);
	}


	ClearAllOnlineThreadAlive();

	WholeProc_EndTime = clock();


	printf("out localization thread\n");



	if (LocTime <= 0)LocTime = 1;
	float LocSpeed = TotalFluoNum*1.0f / LocTime*1000.0f;

	printf("total time : %d ms\n", WholeProc_EndTime - WholeProc_StartTime);

	printf("Molecular Extraction time : %d ms\n", ExtractTime);
	printf("Localization time : %d ms\n", LocTime);
	printf("Rendering time : %d ms\n", RendTime);
	printf("Stastics calc time : %d ms\n", StatTime);
	printf("OnTime calc time : %d ms\n", OnTimeCalcTime);

	

	printf("TotalFluoNum : %d\n", TotalFluoNum);
	printf("Localization speed : %.0f/s\n", LocSpeed);


	IsLocRunning = false;


	return 0;
}

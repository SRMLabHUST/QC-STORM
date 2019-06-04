#include <iostream>
using namespace std;

#include "bfgs_MLE_dll.h"

#include <afx.h>

#include "ListFiles_x64.h"
#include "tinytiffreader.h"

int main(int argc, char* argv[])
{

	int cnt = 0;	

	CStdioFile FileToProc(L"filepath to proc.txt", CFile::modeRead);


	int ImageWidth = 256;
	int ImageHigh = 256;

	int FrameNum = 200;



	float PixelSize_nm = 100.0f;


	int ROISize = 7;
	int LocType = LocType_GS2D;

	int MultiEmitterFitEn = 1;
	int WLEEn = 1;



	// camera info
	LocalizationPara LocPara_Global;

	// camera parameter
	LocPara_Global.QE = 0.77; 

	LocPara_Global.Offset = 100;
	LocPara_Global.KAdc = 0.45;



	LocPara_Global.ROISize = ROISize;
	LocPara_Global.LocType = LocType;

	LocPara_Global.MultiEmitterFitEn = MultiEmitterFitEn;
	LocPara_Global.WLEEn = WLEEn;



	LocPara_Global.ConsecFit_DistanceTh_nm = 80.0; // nm
	LocPara_Global.ConsecFitEn = 0; // better to enable

	LocPara_Global.PixelSize = PixelSize_nm;
	LocPara_Global.PixelZoom = 5;

	LocPara_Global.SNR_th = 5;


	// 3D imaging
	LocPara_Global.ZDepthCorrFactor=1.0f;

	// calibration of sigma X >= sigma Y
	LocPara_Global.p4_XGY = 0.0287;
	LocPara_Global.p3_XGY = -0.3324;
	LocPara_Global.p2_XGY = 2.7762;
	LocPara_Global.p1_XGY = -59.0184;
	LocPara_Global.p0_XGY = 3.2860;


	// calibration of sigma X < sigma Y
	LocPara_Global.p4_XLY =  0.0056;
	LocPara_Global.p3_XLY =  0.2122;
	LocPara_Global.p2_XLY =  -0.4926;
	LocPara_Global.p1_XLY =  -57.1170;
	LocPara_Global.p0_XLY =  3.5634;



	int fcnt = 0;


	CString ImgName;
	TinyTIFFReaderFile* tiffr = NULL;
	unsigned short *h_RawImg = new unsigned short[ImageWidth * ImageHigh];

	cudaStream_t loc_stream1;

	LDROIExtractData_TypeDef LDROIExtractData;
	LDLocData_TypeDef LDLocData;
	ConsecutiveFit_TypeDef ConsecutiveFitData;

	ZeroLocalizationsRemovel_TypeDef ZeroLocRemovel;

	FluoStatisticData_TypeDef FluoStatData;


	int leastPriority, greatestPriority;
	cudaDeviceGetStreamPriorityRange(&leastPriority, &greatestPriority);

	// loc_stream1 declear at bfgs lib
	CreatStreamWithPriority(&loc_stream1, leastPriority);



	int RecFluoNumTh = PointNumTh;


	char pFilePath[1024];
	char CurFileName[1024];


	int GroupCnt = 0;

	CString curPath;
	while (FileToProc.ReadString(curPath))
	{
		if (curPath.GetLength() < 5)continue;

		GroupCnt++;

		wprintf(L"CurFile:%s\n", curPath);


		CString MultiFitStr;
		if (LocPara_Global.MultiEmitterFitEn) MultiFitStr = L"_M";
		else MultiFitStr = L"_S";

		CString ConsecFitStr;
		if (LocPara_Global.ConsecFitEn) ConsecFitStr.Format(L"_Consec%.0fnm", LocPara_Global.ConsecFit_DistanceTh_nm);
		else ConsecFitStr = L"";


		CString PostFix;
		PostFix.Format(_T("_result%dD%d%s%s.txt"), LocPara_Global.LocType + 2, LocPara_Global.ROISize, MultiFitStr, ConsecFitStr);

		CString oFileName = curPath;
		oFileName = oFileName + PostFix;

		int FluoNum=0;
		int TotalFluoNum = 0;
		int CurFrameNum;


		vector<string> files;
		WideCharToMultiByte(CP_ACP, 0, curPath.GetBuffer(), -1, pFilePath, 1024, NULL, NULL);



		// write to file
		float *WriteLocArry = NULL;
		int WriteLocNum = 0;

		CFile LocFile(oFileName, CFile::modeCreate | CFile::modeWrite);

		int IsBreak = 0;
		fcnt = 0;

//		float SumRatio_WLE = 0;
		float SumRatio_1 = 0;
		float SumRatio_2 = 0;
		float SumRatio_3 = 0;
		float SumNum = 0;

		// ratio of ROI by last fitting modality to all detected ROIs
		float FitRatio_Final_1E= 0; // single molecule fitting ratio
		float FitRatio_Final_2E= 0; // two emitter fitting ratio
		float FitRatio_Final_3E= 0; // three emitter fitting ratio
		float FitRatio_Final_4E= 0; // four or more emitter fitting ratio


		float SumDensity = 0;

		while (1)
		{

			sprintf(CurFileName, "%s\\img%d.tif", pFilePath, fcnt + 1);
			fcnt++;


//			printf("%s\n", CurFileName);

			IsBreak = (fcnt == FrameNum);


			if (fcnt % 100 == 0)printf("f%d ", fcnt);

			tiffr = TinyTIFFReader_open(CurFileName);

			if (fcnt == 1)
			{
				ImageWidth = TinyTIFFReader_getWidth(tiffr);
				ImageHigh = TinyTIFFReader_getHeight(tiffr);

				LocPara_Global.ImageWidth = ImageWidth;
				LocPara_Global.ImageHigh = ImageHigh;

				LocPara_Global.UpdateSRImageSize();


				LDROIExtractData.Init(LocPara_Global);

				LDLocData.Init(LocPara_Global);
				ConsecutiveFitData.Init();

				ZeroLocRemovel.Init();

				FluoStatData.Init();

				printf("image size:%d %d %d\n", ImageWidth, ImageHigh, FrameNum);
			}


			if (tiffr)
			{
				TinyTIFFReader_getSampleData(tiffr, h_RawImg, 0); // get image data
			}
			else
			{
				printf("read error:%d\n", fcnt);
			}
			TinyTIFFReader_close(tiffr);
			
			// processing

			memcpy(LDROIExtractData.h_RawImg, h_RawImg, ImageWidth*ImageHigh * 2);

			// subregion extraction for current image
			LDROIExtractData.ExtractMolecules(LDROIExtractData.h_RawImg, ImageSource_CPU_Pinned, LocPara_Global, fcnt, 1, loc_stream1);



			// 
			// accumulate enough molecule for fast localization
			if ((LDROIExtractData.GetAccumulatedROINum() >= RecFluoNumTh) || IsBreak)
			{
				int FirstFrame = LDROIExtractData.FirstFrame;
				int EndFrame = LDROIExtractData.EndFrame;

				FluoNum = LDROIExtractData.GetAccumulatedROINum();
				LDROIExtractData.ResetROINum();

				TotalFluoNum += FluoNum;

				// localization
				LDLocData.BFGS_MLELocalization(LDROIExtractData.h_ImageROI, LDROIExtractData.Get_h_WLEPara(), LocPara_Global, FluoNum, loc_stream1);


				// write localization data into file
				WriteLocArry = LDLocData.h_LocArry;
				WriteLocNum = LDLocData.oValidFluoNum;



				// remove invalid molecules
				ZeroLocRemovel.RemoveZeroLocalizations(LDLocData.h_LocArry, LDLocData.oValidFluoNum, 1, FirstFrame, EndFrame, loc_stream1);

				WriteLocArry = ZeroLocRemovel.h_LocArry;
				WriteLocNum = ZeroLocRemovel.ValidFluoNum;

				// get statistic information
				FluoStatData.GetStatisticalInf(WriteLocArry, LocPara_Global, WriteLocNum, loc_stream1);


				// write localization data into file
				LocFile.Write(WriteLocArry, WriteLocNum * OutParaNumGS2D * sizeof(float));



				if (WriteLocNum > 4000)
				{

//					SumRatio_WLE += LDLocData.FitRatio_WLE; // single molecule fit ratio

					SumRatio_1 += LDLocData.FitRatio_1E;
					SumRatio_2 += LDLocData.FitRatio_2E;
					SumRatio_3 += LDLocData.FitRatio_3E;


					FitRatio_Final_1E += LDLocData.FitRatio_Final_1E; // single molecule fitting ratio
					FitRatio_Final_2E += LDLocData.FitRatio_Final_2E; // two emitter fitting ratio
					FitRatio_Final_3E += LDLocData.FitRatio_Final_3E; // three emitter fitting ratio
					FitRatio_Final_4E += LDLocData.FitRatio_Final_4E; // four or more emitter fitting ratio


					SumDensity += FluoStatData.MeanLocDensity2D;

					SumNum += 1;

				}


			}


			if (IsBreak)break;


		}


//		SumRatio_WLE /= SumNum;
		SumRatio_1 /= SumNum;
		SumRatio_2 /= SumNum;
		SumRatio_3 /= SumNum;
		SumDensity /= SumNum;

		FitRatio_Final_1E /= SumNum;
		FitRatio_Final_2E /= SumNum;
		FitRatio_Final_3E /= SumNum;
		FitRatio_Final_4E /= SumNum;

//		printf("\nWLE Ratio: %.4f\n", SumRatio_WLE);
		printf("\nMFitRatio_1E 2E 3E: %.4f %.4f %.4f\n", SumRatio_1, SumRatio_2, SumRatio_3);

		printf("\nFinal FitRatio_1E 2E 3E 4E: %.4f %.4f %.4f %.4f\n", FitRatio_Final_1E, FitRatio_Final_2E, FitRatio_Final_3E, FitRatio_Final_4E);


		printf("Locdensity: %.4f\n", SumDensity);


		LocFile.Close();

		printf("total fluo:%d\n", TotalFluoNum);


		LDROIExtractData.Deinit();
		LDLocData.Deinit(LocPara_Global);
		ConsecutiveFitData.Deinit();
		ZeroLocRemovel.Deinit();

		FluoStatData.Deinit();
	}



	// free stream
	FreeStream(loc_stream1);
	//	FreeStream(loc_stream2);

	delete[]h_RawImg;

	FileToProc.Close();

	system("pause");

	return 0;
}


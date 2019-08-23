#include <iostream>
using namespace std;

#include "bfgs_MLE_dll.h"

#include <afx.h>

#include "tinytiffreader.h"

#include <time.h>


int main(int argc, char* argv[])
{

	float PixelSize_nm = 100.0f;


	int ROISize = 7;
	int LocType = LocType_GS2D;

	int MultiEmitterFitEn = 1;
	int WLEEn = 1;



	// camera info
	LocalizationPara LocPara_Global;

	// camera parameter
	LocPara_Global.Offset = 100;
	LocPara_Global.KAdc = 0.45;
	LocPara_Global.QE = 0.72;



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
	LocPara_Global.ZDepthCorrFactor = 1.0f;

	// calibration of sigma X >= sigma Y
	LocPara_Global.p4_XGY = 0.0174;
	LocPara_Global.p3_XGY = -0.2339;
	LocPara_Global.p2_XGY = 2.6260;
	LocPara_Global.p1_XGY = -59.2057;
	LocPara_Global.p0_XGY = 3.5972;


	// calibration of sigma X < sigma Y
	LocPara_Global.p4_XLY = -0.0476;
	LocPara_Global.p3_XLY = -0.3630;
	LocPara_Global.p2_XLY = -2.2169;
	LocPara_Global.p1_XLY = -58.2213;
	LocPara_Global.p0_XLY = 3.9583;




	LDROIExtractData_TypeDef LDROIExtractData;
	LDLocData_TypeDef LDLocData;
	ConsecutiveFit_TypeDef ConsecutiveFitData;

	ZeroLocalizationsRemovel_TypeDef ZeroLocRemovel;



	cudaStream_t loc_stream1;
	
	// loc_stream1 declear at bfgs lib
	CreatStream(&loc_stream1);


	int FluoNum;
	int TotalFluoNum;


	int RecFluoNumTh = PointNumTh;


	CString curPath;
	CString oFileName;

	CStdioFile FileToProc(L"filepath to proc.txt", CFile::modeRead);

	char buf[1024];

	while (FileToProc.ReadString(curPath))
	{
		if (curPath.GetLength() < 5)continue;

		if (curPath.Find(L"#") >= 0) continue;

		oFileName = curPath;
		oFileName.TrimRight(L".tif");
		oFileName = oFileName + L"_LocArray.txt";

		CString ImgName = curPath; //+L"\\MMStack_Pos0.ome.tif"


		wprintf(L"CurFile:%s\n", curPath);


		WideCharToMultiByte(CP_ACP, 0, ImgName.GetBuffer(), -1, buf, 1024, NULL, NULL);


		TinyTIFFReaderFile* tiffr = TinyTIFFReader_open(buf);

		if (!tiffr)
		{
			printf("image open error\n");
			TinyTIFFReader_close(tiffr);
	
			continue;

		}

		// write to file
		float *WriteLocArry = NULL;
		int WriteLocNum = 0;

		CFile LocFile(oFileName, CFile::modeCreate | CFile::modeWrite);


		TotalFluoNum = 0;


		int TotalFrame = TinyTIFFReader_countFrames(tiffr);
		printf("total image frame:%d\n", TotalFrame);

		int ImageWidth = TinyTIFFReader_getWidth(tiffr);
		int ImageHigh = TinyTIFFReader_getHeight(tiffr);
		printf("raw image size:%d x %d\n", ImageWidth, ImageHigh);


		LocPara_Global.ImageWidth = ImageWidth;
		LocPara_Global.ImageHigh = ImageHigh;

		LocPara_Global.UpdateSRImageSize();



		LDROIExtractData.Init(LocPara_Global);

		LDLocData.Init(LocPara_Global);
		ConsecutiveFitData.Init();

		ZeroLocRemovel.Init();



		bool IsBreak = 0;


		int fcnt = 0;
		

		while (1)
		{
			fcnt++;

			IsBreak = (fcnt == TotalFrame);

			if (fcnt % (TotalFrame / 10) == 0)printf("progress: %.2f%%\n", fcnt * 100.0f / TotalFrame);


			TinyTIFFReader_getSampleData(tiffr, LDROIExtractData.h_RawImg, 0); // get image data
			TinyTIFFReader_readNext(tiffr);


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


				// write localization data into file
				LocFile.Write(WriteLocArry, WriteLocNum * OutParaNumGS2D * sizeof(float));


			}

			if (IsBreak)break;

		}



		printf("total fluo:%d\n\n\n", TotalFluoNum);



		TinyTIFFReader_close(tiffr);

		LDROIExtractData.Deinit();
		LDLocData.Deinit(LocPara_Global);
		ConsecutiveFitData.Deinit();
		ZeroLocRemovel.Deinit();


		LocFile.Close();

	}



	// free stream
	FreeStream(loc_stream1);


	FileToProc.Close();

	printf("all process finished\n");

	system("pause");

	return 0;
}


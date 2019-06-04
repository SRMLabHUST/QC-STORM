#include <iostream>
using namespace std;
#include <afx.h>


#include "bfgs_MLE_dll.h"

#include "tinytiffreader.h"
#include "tinytiffwriter.h"
#include "CppFileWrapper.h"

void SRImgRendering(CString FilePath, ImageRenderData_TypeDef &RendData, LocalizationPara& LocPara_Global, cudaStream_t loc_stream1);


#define RenderingMode	RenderingMode_FittedPhoton_1Pixel

int main()
{


	// in multi ROI, set false, and set real image size
	bool AutoFindImgSize = false;
	 
	int ImageWidth = 1024;
	int ImageHigh = 1024;


	float RenderingPixelSize = 35; // nm

	float PixelSize_nm = 108.3f; // 97 108.3f


	// camera info
	LocalizationPara LocPara_Global;

	// camera parameter
	LocPara_Global.Offset = 100;
	LocPara_Global.KAdc = 0.45;
	LocPara_Global.QE = 0.72;  // be careful
	LocPara_Global.ReadNoise_e = 1.3;// e-

	LocPara_Global.ROISize = 7;
	LocPara_Global.LocType = 0;



	LocPara_Global.ConsecFit_DistanceTh_nm = 80.0; // nm
	LocPara_Global.ConsecFitEn = 0; // better to enable

	LocPara_Global.PixelSize = PixelSize_nm; // raw pixel size
	LocPara_Global.PixelZoom = LocPara_Global.PixelSize / RenderingPixelSize; //

	LocPara_Global.SNR_th = 5;



	CStdioFile FileToProc(L"filepath to proc.txt", CFile::modeRead);

	cudaStream_t loc_stream1;

	int leastPriority, greatestPriority;
	cudaDeviceGetStreamPriorityRange(&leastPriority, &greatestPriority);

	// loc_stream1 declear at bfgs lib
	CreatStreamWithPriority(&loc_stream1, leastPriority);




	CString curPath;

	int GroupCnt = 0;

	while (FileToProc.ReadString(curPath))
	{
		if (curPath.GetLength() < 5)continue;
		wprintf(L"CurFile:%s\n", curPath);

		GroupCnt++;

		CString SRFileName = curPath;
		SRFileName.TrimRight(L".txt");
		SRFileName.Format(L"%s_Rerend_%.2fnm.tif", SRFileName, RenderingPixelSize);



		wstring curPath_ws = curPath.GetBuffer();
		string curPath_s = EncodeConvert::ws2s(curPath_ws);

		cout << "curPath_s: " << curPath_s << endl;


		// get image size
		if (AutoFindImgSize)
		{
			SRDriftCorrData_TypeDef::GetMaxImgSize(curPath_s, &ImageWidth, &ImageHigh);

		}

		int TotalFrame = SRDriftCorrData_TypeDef::GetTotalFrame(curPath_s);


		printf("find img size:%d x %d x %d f\n", ImageWidth, ImageHigh, TotalFrame);


		LocPara_Global.ImageWidth = ImageWidth;
		LocPara_Global.ImageHigh = ImageHigh;
		LocPara_Global.UpdateSRImageSize();

		ImageRenderData_TypeDef RendData;

		RendData.Init(LocPara_Global, 4, 4); // no sr display in c++
		float* h_RendFloatImage2D = new float[LocPara_Global.SRImageWidth*LocPara_Global.SRImageHigh];



		SRImgRendering(curPath, RendData, LocPara_Global, loc_stream1);



		// create 2d image 
		cudaMemcpyAsync(h_RendFloatImage2D, RendData.d_SRIntensityImg, LocPara_Global.SRImageWidth* LocPara_Global.SRImageHigh*sizeof(float), cudaMemcpyDeviceToHost, loc_stream1);
		WaitGPUStream(loc_stream1);


		char SRImgFileName[360];
//		sprintf_s(SRImgFileName, "ReRendered sr image_%d-%dnm.tif", GroupCnt, (int)RenderingPixelSize);

		WideCharToMultiByte(CP_OEMCP, NULL, SRFileName.GetBuffer(), -1, SRImgFileName, 360, NULL, FALSE);

		TinyTIFFFile* tif = TinyTIFFWriter_open(SRImgFileName, 32, LocPara_Global.SRImageWidth, LocPara_Global.SRImageHigh); // 32 bit float image
		TinyTIFFWriter_writeImage(tif, h_RendFloatImage2D);

		TinyTIFFWriter_close(tif);


		RendData.Deinit(LocPara_Global);
		delete[]h_RendFloatImage2D;

	}


	return 0;
}

void SRImgRendering(CString FilePath, ImageRenderData_TypeDef &RendData, LocalizationPara& LocPara_Global, cudaStream_t loc_stream1)
{
	// clear old image
	RendData.ResetFillImgTop(LocPara_Global);

	CFile LocFile(FilePath, CFile::modeRead);

	// image rendering
	int TotalFluoNum = LocFile.GetLength() / sizeof(float) / OutParaNumGS2D; //
	printf("fluo num:%d\n", TotalFluoNum);

	const int BatchProcNum = PointNumTh;
	int ProcNum = (TotalFluoNum + BatchProcNum - 1) / BatchProcNum;

	int CurFluoNum = 0;
	for (int pcnt = 0; pcnt < ProcNum; pcnt++)
	{
		CurFluoNum = min(BatchProcNum, TotalFluoNum - pcnt*BatchProcNum);

		LocFile.Read(RendData.h_LocArry, CurFluoNum*OutParaNumGS2D*sizeof(float));

		RendData.FluoRenderTop(RendData.h_LocArry, LocPara_Global, RenderingMode, 0, CurFluoNum, loc_stream1);

	}

	LocFile.Close();
}


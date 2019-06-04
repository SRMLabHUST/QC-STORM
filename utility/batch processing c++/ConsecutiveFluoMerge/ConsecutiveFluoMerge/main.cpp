#include <iostream>
using namespace std;

#include <stdio.h>
#include <afx.h>


#include "cuda_runtime.h"
#include "device_launch_parameters.h"

#include "ConsecutiveFluoMerge.h"




CString GetFatherPath(CString ipath);
CString GetFileName(CString ipath);
CString GetFileName1(CString ipath, CString PostFix);


int main()
{

	// consecutive filter


	float PixelSize_nm = 97; //97 108.3 nm
	float ConsecFit_DistanceTh_nm = 50; // nm

	float QE = 0.77;

	float Offset = 100;
	float KAdc = 0.45;


	float ReadNoise_e = 0.1;
	



	float ConsecFilter_DistanceTh_pixel = ConsecFit_DistanceTh_nm / PixelSize_nm;

	// camera info
	LocalizationPara LocPara_Global;
	// camera parameter
	LocPara_Global.Offset = Offset;
	LocPara_Global.KAdc = KAdc;
	LocPara_Global.QE = QE;  // be careful
	LocPara_Global.ReadNoise_e = ReadNoise_e;

	LocPara_Global.ROISize = 7;
	LocPara_Global.LocType = 0;



	LocPara_Global.ConsecFit_DistanceTh_nm = ConsecFit_DistanceTh_nm; // nm
	LocPara_Global.ConsecFitEn = 0; // better to enable

	LocPara_Global.PixelSize = PixelSize_nm;
	LocPara_Global.PixelZoom = 5;

	LocPara_Global.SNR_th = 5;

	//

	CStdioFile FileToProc(L"filepath to proc.txt", CFile::modeRead);


	CString wStr1;


	cudaStream_t loc_stream1;

	cudaStreamCreate(&loc_stream1);


	// mean distance calculation for a group of localizations
	ConsecutiveFluoMerger_TypeDef NyqConsecutiveFilter;


	CString FileName;

	int ValidFileNum = 0;


	int fcnt = 0;

	// for each file, only a single file supported
	while (FileToProc.ReadString(FileName))
	{
		if (FileName.GetLength() < 5)continue;
		if (FileName.Find(L"#") >= 0) continue;

		wprintf(L"file name:%s\n", FileName);

		ValidFileNum++;

		CFile LocFile;
		bool suc = LocFile.Open(FileName, CFile::modeRead);
		if (!suc)
		{
			printf("file open error\n");
			break;
		}

		CString InfFileName;
		InfFileName.Format(L"_ConsecMerge_Th%.2fnm.txt", ConsecFit_DistanceTh_nm);

		InfFileName = GetFatherPath(FileName) + GetFileName1(FileName, L".txt") + InfFileName;

		CFile oLocFile(InfFileName, CFile::modeCreate | CFile::modeWrite);



		int TotalFluoNum = LocFile.GetLength() / OutParaNumGS2D / sizeof(float);
		printf("TotalFluoNum:%d\n", TotalFluoNum);


		NyqConsecutiveFilter.Init(TotalFluoNum);


		// get data
		LocFile.Read(NyqConsecutiveFilter.h_LocArry, LocFile.GetLength());


		NyqConsecutiveFilter.MergeConsecutiveFluo(NyqConsecutiveFilter.h_LocArry, TotalFluoNum, LocPara_Global, 0, ConsecFilter_DistanceTh_pixel, loc_stream1);

		oLocFile.Write(NyqConsecutiveFilter.h_LocArry, TotalFluoNum * OutParaNumGS2D * sizeof(float));



		LocFile.Close();
		oLocFile.Close();

		//

		NyqConsecutiveFilter.DeInit();


	}


	cudaStreamDestroy(loc_stream1);

	return 0;

}


CString GetFatherPath(CString ipath)
{

	int id = ipath.ReverseFind(L'\\');

	CString opath = ipath.Left(id + 1);

	return opath;
}

CString GetFileName(CString ipath)
{
	int id = ipath.ReverseFind(L'\\');
	int selLen = ipath.GetLength() - id - 1;

	CString fname = ipath.Right(selLen);


	return fname;

}

CString GetFileName1(CString ipath, CString PostFix)
{
	int id = ipath.ReverseFind(L'\\');
	int selLen = ipath.GetLength() - id - 1;

	CString fname = ipath.Right(selLen);

	fname.TrimRight(PostFix);

	return fname;

}

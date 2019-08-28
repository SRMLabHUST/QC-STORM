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

#include "LocResource.h"
#include "OnlineLocalizationLD.h"
#include "StatInfDisplay.h"

#include "BatchLocalization.h"

// control whether thread running
volatile bool OnlineLocAlive = false;
volatile bool OnlineRendAlive = false;


// thread state, 0:idle, 1:running, 2:is busy
volatile bool IsLocRunning = false;
volatile bool IsRendRunning = false;
volatile bool IsBatchLocRunning = false;


// image for display and save in ImageJ
float *h_RendFloatImage2D = NULL;



// image info
CString ImageName; // image path of current processing
CString LocFileName; // localization data file name

cudaStream_t loc_stream1;

cudaStream_t render_stream1;



concurrent_queue<qImgData> ImgDataQueue;


LocalizationPara LocPara_Global;

LDROIExtractData_TypeDef LDROIExtractData;

LDLocData_TypeDef LDLocData;

ZeroLocalizationsRemovel_TypeDef ZeroLocRemovel;

ConsecutiveFit_TypeDef ConsecutiveFitData;

DH3D_MoleculePair_TypeDef DH3D_MoleculePair;

FluoStatisticData_TypeDef FluoStatData;
ImageRenderData_TypeDef RendData;


int GPUID_1Best = 0;
int GPUID_2Best = 0;



cudaStream_t Resolution_stream1;

NyqDimensionDensityCalc_TypeDef DimensionDensityCalc;
SpatialResolutionCalc_TypeDef SpatialResolutionCalc;


volatile bool IsLocResourceAllocated = false;

LocalizationPara LocPara_Last;


void StartLocalizationThread()
{

	OnlineLocAlive = true;
	OnlineRendAlive = true;


	InitAllLocResource(0);


	AfxBeginThread(th_OnlineLocalizationLD, NULL);

	AfxBeginThread(th_OnlineRendDispLD, NULL);


	if (LocPara_Global.SpatialResolutionCalcEn)
	{
		AfxBeginThread(th_OnlineSpatialResolutionCalc, NULL);
	}
}

void StopLocThread()
{
	OnlineLocAlive = false;

}

void FinishLocalizationThread()
{
	OnlineLocAlive = false;

	DeinitAllLocResource(0);

}




void CreateFeedImgMemory(qImgData & CurImgInf, unsigned short* pImg, int ImageWidth, int ImageHigh, int FrameNum)
{

	CurImgInf.BatchFrameNum = FrameNum;


	cudaError_t err = cudaErrorMemoryAllocation;

	int BatchedImgSize = ImageWidth * ImageHigh * FrameNum;

	// try allocate CPU memory
	err = cudaMallocHost((void **)&CurImgInf.pImgData, BatchedImgSize * sizeof(short));


	if (err == cudaSuccess)
	{
		//	printf("cuda suc:%s\n", str);
		CurImgInf.ImageSource = ImageSource_CPU_Pinned;
	}
	else
	{
		while (1)
		{
			// if GPU memory is allocate error, then try CPU memory
			CurImgInf.pImgData = new unsigned short[BatchedImgSize];
			CurImgInf.ImageSource = ImageSource_CPU_Normal;

			if (CurImgInf.pImgData != NULL)break;
		}
	}


	memcpy(CurImgInf.pImgData, pImg, BatchedImgSize * sizeof(short));

}

// for nromal imageJ process
void SetLocDataFileName()
{
	// set localization file name
	LocFileName = ImageName;
	LocFileName.TrimRight(_T(".tif"));

	CString MultiFitStr;
	if (LocPara_Global.MultiEmitterFitEn) MultiFitStr = L"_M";
	else MultiFitStr = L"_S";

	CString ConsecFitStr;
	if (LocPara_Global.ConsecFitEn) ConsecFitStr.Format(L"_Consec%.0fnm", LocPara_Global.ConsecFit_DistanceTh_nm);
	else ConsecFitStr = L"";


	CString PostFix;
	PostFix.Format(_T("_result%dD%d%s%s.txt"), LocPara_Global.LocType + 2, LocPara_Global.ROISize, MultiFitStr, ConsecFitStr);

	LocFileName = LocFileName + PostFix;

}

// for batch processing
void SetLocDataFileName_BatchProc()
{
	// set localization file name
	LocFileName = ImageName;
	LocFileName.TrimRight(_T(".tif"));

	CString MultiFitStr;
	if (LocPara_Global.MultiEmitterFitEn) MultiFitStr = L"_M";
	else MultiFitStr = L"_S";

	CString ConsecFitStr;
	if (LocPara_Global.ConsecFitEn) ConsecFitStr.Format(L"_Consec%.0fnm", LocPara_Global.ConsecFit_DistanceTh_nm);
	else ConsecFitStr = L"";


	CString PostFix;
	PostFix.Format(_T("_result%dD%d%s%s.txt"), LocPara_Global.LocType + 2, LocPara_Global.ROISize, MultiFitStr, ConsecFitStr);

	LocFileName = LocFileName + PostFix;

	// modify save path
	LocFileName.TrimLeft(BatchProc_FolderName.c_str());
	LocFileName.TrimLeft(L"\\");
	LocFileName.Replace(L"\\", L"__");

	CString SavePath = BatchProc_SavePath.c_str();
	SavePath = SavePath + L"\\";
	LocFileName = SavePath + LocFileName;

}


void InitAllLocResource(int IsPostprocess)
{
	SelectBestGPU();

	if ((IsLocResourceAllocated == false) || (!LocPara_Global.IsEqual(LocPara_Last)))
	{
		DeinitAllLocResource(IsPostprocess);

		LocPara_Last = LocPara_Global;


		// the GPU 1
		cudaSetDevice(GPUID_1Best);


		// alloc stream with priority
		int leastPriority, greatestPriority;
		cudaDeviceGetStreamPriorityRange(&leastPriority, &greatestPriority);

		cudaStreamCreate(&loc_stream1);
		cudaStreamCreate(&render_stream1);

		if (IsPostprocess == 0)
		{
			// extraction and localization
			LDROIExtractData.Init(LocPara_Global);

		}

		LDLocData.Init(LocPara_Global);

		ZeroLocRemovel.Init();

		DH3D_MoleculePair.Init();

		if (LocPara_Global.ConsecFitEn)
		{
			ConsecutiveFitData.Init(); 
		}


		// rendering and stastical
		FluoStatData.Init();
		RendData.Init(LocPara_Global, 4, 4); // no sr display in c++


		//
		ConvertInfImageBuf = new unsigned char[AxesImgWidth*AxesImgHigh * 4];
		
		if (RendType_Is2D(LocPara_Global.LocType))
		{
			h_RendFloatImage2D = new float[LocPara_Global.SRImageWidth*LocPara_Global.SRImageHigh];
		}
		else
		{
		}

		if (LocPara_Global.SpatialResolutionCalcEn)
		{
			// the GPU 2
			cudaSetDevice(GPUID_2Best);
			
			cudaStreamCreate(&Resolution_stream1);

			// spatial resolution calculation
			DimensionDensityCalc.Init();
			SpatialResolutionCalc.Init();
		}

		IsLocResourceAllocated = true;
	}
}


void DeinitAllLocResource(int IsPostprocess)
{
	SelectBestGPU();

	// for both 2d and 3d


	if (IsLocResourceAllocated == true)
	{

		cudaSetDevice(GPUID_1Best);

		// free stream
		cudaStreamDestroy(loc_stream1);
		cudaStreamDestroy(render_stream1);

		if (IsPostprocess == 0)
		{
			// extraction and localization
			LDROIExtractData.Deinit();

		}

		LDLocData.Deinit(LocPara_Global);

		ZeroLocRemovel.Deinit();

		DH3D_MoleculePair.Deinit();


		if (LocPara_Global.ConsecFitEn)
		{
			ConsecutiveFitData.Deinit();

		}
		// rendering and stastical

		FluoStatData.Deinit();
		RendData.Deinit(LocPara_Global);


		//
		delete[] ConvertInfImageBuf;


		//
		if (RendType_Is2D(LocPara_Global.LocType))
		{
			delete[] h_RendFloatImage2D;
		}
		else
		{
		}


		if (LocPara_Global.SpatialResolutionCalcEn)
		{
			cudaSetDevice(GPUID_2Best);

			cudaStreamDestroy(Resolution_stream1);


			// spatial resolution calculation
			DimensionDensityCalc.DeInit();
			SpatialResolutionCalc.DeInit();
		}

		//
		cudaSetDevice(GPUID_2Best);
		cudaDeviceReset();
		cudaSetDevice(GPUID_1Best);
		cudaDeviceReset();


		IsLocResourceAllocated = false;
	}
}

void SelectBestGPU()
{
	GPUID_1Best = 0;
	GPUID_2Best = 0;


	int DevNum = 0;

	int cnt = 0;
	cudaDeviceProp iProp;


	cudaGetDeviceCount(&DevNum);
	printf("GPU num:%d\n", DevNum);

	if ((DevNum < 1) || (DevNum > 20))
	{
		printf("error: No GPU find\n\n");

		return ;
	}


	int *ProcessorNumArray = new int[DevNum];

	// get processor number per GPU
	for (cnt = 0; cnt < DevNum; cnt++)
	{
		cudaGetDeviceProperties(&iProp, cnt);
		ProcessorNumArray[cnt] = iProp.multiProcessorCount;
	}


	// find the GPU with largest processor number

	for (cnt = 0; cnt < DevNum; cnt++)
	{
		if (ProcessorNumArray[cnt] > ProcessorNumArray[GPUID_1Best])
		{
			GPUID_1Best = cnt;
		}
	}

	// find the send best GPU
	GPUID_2Best = GPUID_1Best;
	ProcessorNumArray[GPUID_1Best] = 0;

	if (DevNum > 1)
	{
		for (cnt = 0; cnt < DevNum; cnt++)
		{
			if (ProcessorNumArray[cnt] > ProcessorNumArray[GPUID_2Best])
			{
				GPUID_2Best = cnt;
			}
		}
	}


	delete[] ProcessorNumArray;

	// print information
	cudaGetDeviceProperties(&iProp, GPUID_1Best);
	printf("1st BestGPUId %d:%s\n", GPUID_1Best, iProp.name);
	printf("Number of multiprocessors: %d\n", iProp.multiProcessorCount);

	if (DevNum > 1)
	{
		cudaGetDeviceProperties(&iProp, GPUID_2Best);
		printf("2nd BestGPUId %d:%s\n", GPUID_2Best, iProp.name);
		printf("Number of multiprocessors: %d\n", iProp.multiProcessorCount);
	}
}

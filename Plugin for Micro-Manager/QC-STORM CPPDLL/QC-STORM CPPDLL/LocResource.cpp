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
#include "StatInfDisplay.h"

#include "LocResource_ZDriftMeasure.h"

int ProcessorID = 0;

volatile bool OnlineLocAlive = false;
volatile bool OnlineRendAlive = false;
volatile bool OnlineFeedbackAlive = false;


volatile bool IsLocRunning = false;
volatile bool IsRendRunning = false;


int DepthMapDispEn = 0;

// image info


CString ResultSavePath;
CString CreateTimeIdxStr;


cudaStream_t loc_stream1;

cudaStream_t render_stream1;



int GPUID_1Best = 0;
int GPUID_2Best = 0;


// image queue from ImageJ/micro-manager to c++ thread
tbb::concurrent_queue<qImgData> ImgDataQueue;

LocalizationPara LocPara_Global;

LDROIExtractData_TypeDef LDROIExtractData;

LDLocData_TypeDef LDLocData;

ZeroLocalizationsRemovel_TypeDef ZeroLocRemovel;

ConsecutiveFit_TypeDef ConsecutiveFitData;


FluoStatisticData_TypeDef FluoStatData;
ImageRenderData_TypeDef RendData;


cudaStream_t Resolution_stream1;

NyqDimensionDensityCalc_TypeDef DimensionDensityCalc;
SpatialResolutionCalc_TypeDef SpatialResolutionCalc;

MarginDataFilter_TypeDef MarginDataFilter;

float FOVOverlapPercent = 0;


void SetAllOnlineThreadAlive()
{
	OnlineLocAlive = true;

	OnlineRendAlive = true;
	OnlineFeedbackAlive = true;

}

void ClearAllOnlineThreadAlive()
{
	OnlineLocAlive = false;

	OnlineRendAlive = false;
	OnlineFeedbackAlive = false;

}


volatile bool IsLocResourceAllocated = false;

LocalizationPara LocPara_Last;



void InitAllLocResource()
{
	SelectBestGPU();

	if ((IsLocResourceAllocated == false) || (!LocPara_Global.IsEqual(LocPara_Last)))
	{
		DeinitAllLocResource();

		LocPara_Last = LocPara_Global;

		// the GPU 1
		cudaSetDevice(GPUID_1Best);

		// alloc stream with priority
		int leastPriority, greatestPriority;
		cudaDeviceGetStreamPriorityRange(&leastPriority, &greatestPriority);

		CreatStreamWithPriority(&render_stream1, leastPriority + 1);
		CreatStreamWithPriority(&loc_stream1, leastPriority);


		// extraction and localization
		LDROIExtractData.Init(LocPara_Global);
		LDLocData.Init(LocPara_Global);

		ZeroLocRemovel.Init();


		if (LocPara_Global.ConsecFitEn)
		{
			ConsecutiveFitData.Init();
		}


		// rendering and stastical
		FluoStatData.Init();
		RendData.Init(LocPara_Global, SR_IMAGE_DISPLAY_WIDTH, SR_IMAGE_DISPLAY_HIGH); // not the same with ImageJ plugin

		
		if (RendType_Is2D(LocPara_Global.LocType))
		{
			h_RendFloatImage2D = new float[LocPara_Global.SRImageWidth*LocPara_Global.SRImageHigh];
		}


		if (LocPara_Global.SpatialResolutionCalcEn)
		{
			cudaSetDevice(GPUID_2Best);
			
			cudaStreamCreate(&Resolution_stream1);
			
			// spatial resolution calculation
			DimensionDensityCalc.Init();
			SpatialResolutionCalc.Init();

			MarginDataFilter.Init();

		}

		// z drift measure by seperated localization

		ZDriftCtl.Init(LocPara_Global);


		IsLocResourceAllocated = true;
	}

}


void DeinitAllLocResource()
{
	SelectBestGPU();

	if (IsLocResourceAllocated == true)
	{
		// the GPU 1
		cudaSetDevice(GPUID_1Best);


		// free stream
		FreeStream(render_stream1);
		FreeStream(loc_stream1);

		// extraction and localization
		LDROIExtractData.Deinit();
		LDLocData.Deinit(LocPara_Global);

		ZeroLocRemovel.Deinit();


		if (LocPara_Global.ConsecFitEn)
		{
			ConsecutiveFitData.Deinit();
		}


		// rendering and stastical

		FluoStatData.Deinit();
		RendData.Deinit(LocPara_Global);


		//
		if (RendType_Is2D(LocPara_Global.LocType))
		{
			delete[] h_RendFloatImage2D;
		}


		// z drift measure by seperated localization
		ZDriftCtl.Deinit(LocPara_Global);


		if (LocPara_Global.SpatialResolutionCalcEn)
		{
			cudaSetDevice(GPUID_2Best);
			
			cudaStreamDestroy(Resolution_stream1);
			
			// spatial resolution calculation
			DimensionDensityCalc.DeInit();
			SpatialResolutionCalc.DeInit();

			MarginDataFilter.DeInit();

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

	if ((DevNum < 1) || (DevNum > 10))
	{
		printf("error: No GPU find\n\n");

		return;
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


	cudaGetDeviceProperties(&iProp, GPUID_2Best);
	printf("2nd BestGPUId %d:%s\n", GPUID_2Best, iProp.name);
	printf("Number of multiprocessors: %d\n", iProp.multiProcessorCount);
}

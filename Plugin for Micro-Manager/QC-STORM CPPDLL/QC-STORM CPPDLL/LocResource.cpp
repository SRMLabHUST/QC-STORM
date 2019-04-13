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


volatile bool IsLocResourceAllocated = false;


// image info


CString ResultSavePath;
CString CreateTimeIdxStr;


cudaStream_t loc_stream1;
cudaStream_t loc_stream2;

cudaStream_t render_stream1;
cudaStream_t render_stream2;

int SelectedGPUID;

// image queue from ImageJ/micro-manager to c++ thread
tbb::concurrent_queue<qImgData> ImgDataQueue;

LocalizationPara LocPara_Global;

LDROIExtractData_TypeDef LDROIExtractData;

LDLocData_TypeDef LDLocData;

ZeroLocalizationsRemovel_TypeDef ZeroLocRemovel;

ConsecutiveFit_TypeDef ConsecutiveFitData;


FluoStatisticData_TypeDef FluoStatData;
ImageRenderData_TypeDef RendData;



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


void InitAllLocResource()
{

	if (IsLocResourceAllocated == false)
	{
		SelectedGPUID = SelectBestGPU();

		// alloc stream with priority
		int leastPriority, greatestPriority;
		cudaDeviceGetStreamPriorityRange(&leastPriority, &greatestPriority);

		CreatStreamWithPriority(&render_stream1, leastPriority + 1);
		CreatStreamWithPriority(&render_stream2, leastPriority + 1);

		CreatStreamWithPriority(&loc_stream1, leastPriority);
		CreatStreamWithPriority(&loc_stream2, leastPriority);


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
		else
		{
		}

		if (LocPara_Global.SpatialResolutionCalcEn)
		{
			// spatial resolution calculation
			DimensionDensityCalc.Init(MAX_FLUO_NUM_PER_GROUP, LocPara_Global.ImagesPerGroup);
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

	if (IsLocResourceAllocated == true)
	{
		// for both 2d and 3d
		cudaSetDevice(SelectedGPUID);

		// free stream
		FreeStream(render_stream1);
		FreeStream(render_stream2);

		FreeStream(loc_stream1);
		FreeStream(loc_stream2);

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


		if (LocPara_Global.SpatialResolutionCalcEn)
		{
			// spatial resolution calculation
			DimensionDensityCalc.DeInit();
			SpatialResolutionCalc.DeInit();

			MarginDataFilter.DeInit();

		}

		// z drift measure by seperated localization
		ZDriftCtl.Deinit(LocPara_Global);


//		cudaDeviceReset();

		IsLocResourceAllocated = false;
	}



}

int SelectBestGPU()
{
	int DevCnt;

	int cnt = 0;
	cudaDeviceProp iProp;


	int BestGPUId = 0;

	cudaGetDeviceCount(&DevCnt);
	printf("GPU num:%d\n", DevCnt);

	if ((DevCnt < 1) || (DevCnt > 10))
	{
		printf("GPU load error\n");

		return -1;
	}


	int *ProcNumArray = new int[DevCnt];
	for (cnt = 0; cnt < DevCnt; cnt++)
	{
		cudaGetDeviceProperties(&iProp, cnt);
		ProcNumArray[cnt] = iProp.multiProcessorCount;
	}


	// find the device with largest multi processor count
	BestGPUId = 0;
	for (cnt = 0; cnt < DevCnt; cnt++)
	{
		if (ProcNumArray[cnt] > ProcNumArray[BestGPUId])
		{
			BestGPUId = cnt;
		}
	}


	delete[] ProcNumArray;

	cudaGetDeviceProperties(&iProp, BestGPUId);
	printf("BestGPUId %d:%s\n", BestGPUId, iProp.name);
	printf("Number of multiprocessors: %d\n", iProp.multiProcessorCount);

	cudaSetDevice(BestGPUId);

	return BestGPUId;
}

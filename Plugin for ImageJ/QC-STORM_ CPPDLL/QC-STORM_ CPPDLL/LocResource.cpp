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



volatile bool OnlineLocAlive = false;
volatile bool OnlineRendAlive = false;

// 0:idle, 1:running, 2:is busy
volatile bool IsLocRunning = false;
volatile bool IsRendRunning = false;


bool IsLocResourceAllocated = false;

// image for display and save in ImageJ
float *h_RendFloatImage2D = NULL;

// image info

CString ImageName;


cudaStream_t loc_stream1;
cudaStream_t loc_stream2;

cudaStream_t render_stream1;
cudaStream_t render_stream2;



concurrent_queue<qImgData> ImgDataQueue;

LocalizationPara LocPara_Global;
LDROIExtractData_TypeDef LDROIExtractData;

LDLocData_TypeDef LDLocData;

ZeroLocalizationsRemovel_TypeDef ZeroLocRemovel;

ConsecutiveFit_TypeDef ConsecutiveFitData;


FluoStatisticData_TypeDef FluoStatData;
ImageRenderData_TypeDef RendData;


int SelectedGPUID;


int ImagesPerGroup;
float StrucuteSize_2D;
float StrucuteSize_3D;


NyqDimensionDensityCalc_TypeDef DimensionDensityCalc;
SpatialResolutionCalc_TypeDef SpatialResolutionCalc;



void InitAllLocResource(int IsPostprocess)
{
	SelectedGPUID = SelectBestGPU();

	if (IsLocResourceAllocated == false)
	{
		// alloc stream with priority
		int leastPriority, greatestPriority;
		cudaDeviceGetStreamPriorityRange(&leastPriority, &greatestPriority);

		CreatStreamWithPriority(&render_stream1, leastPriority + 1);
		CreatStreamWithPriority(&render_stream2, leastPriority + 1);

		CreatStreamWithPriority(&loc_stream1, leastPriority);
		CreatStreamWithPriority(&loc_stream2, leastPriority);

		if (IsPostprocess == 0)
		{
			// extraction and localization
			LDROIExtractData.Init(LocPara_Global);

		}

		LDLocData.Init(LocPara_Global);

		ZeroLocRemovel.Init();


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
			// spatial resolution calculation
			DimensionDensityCalc.Init(MAX_FLUO_NUM_PER_GROUP, LocPara_Global.ImagesPerGroup);
			SpatialResolutionCalc.Init();
		}

		IsLocResourceAllocated = true;
	}

}


void DeinitAllLocResource(int IsPostprocess)
{
	// for both 2d and 3d
	cudaSetDevice(SelectedGPUID);


	if (IsLocResourceAllocated == true)
	{

		// free stream
		FreeStream(render_stream1);
		FreeStream(render_stream2);

		FreeStream(loc_stream1);
		FreeStream(loc_stream2);

		if (IsPostprocess == 0)
		{
			// extraction and localization
			LDROIExtractData.Deinit();

		}

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
			// spatial resolution calculation
			DimensionDensityCalc.DeInit();
			SpatialResolutionCalc.DeInit();

		}


		//
		cudaDeviceReset();

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

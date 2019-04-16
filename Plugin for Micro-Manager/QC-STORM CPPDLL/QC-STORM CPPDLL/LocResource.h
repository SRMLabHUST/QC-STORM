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

#pragma once


// molecule finding, MLE localization,rendering,statistical processing
#include "bfgs_MLE_dll.h"

// nyquist and spatial resolution calculation
#include "SpatialResolutionCalc.h"


// for multi-thread program useage
#include "tbb.h"
//using namespace tbb;




struct qImgData
{
	unsigned short *pImgData;
	int ImageSource;
	int BatchFrameNum;
};

// to receive feedback cmd from master to slave
class ThreadCmdProcessState
{
private:
	volatile unsigned int CmdCnt; // number for asked commands
	volatile unsigned int ExeCnt; // number of executed commands

public:

	ThreadCmdProcessState()
	{
		Reset();
	}

	void Reset()
	{
		CmdCnt = 0;
		ExeCnt = 0;
	}

	void MakeAProcess()
	{
		CmdCnt++;
	}

	int HaveProcessWait()
	{
		return (ExeCnt != CmdCnt);
	}

	void ProcessFinish()
	{
		 ExeCnt = CmdCnt;
	}
};


extern int ProcessorID;

extern volatile bool OnlineLocAlive;
extern volatile bool OnlineRendAlive;
extern volatile bool OnlineFeedbackAlive;


// proc status

extern volatile  bool IsLocRunning;
extern volatile  bool IsRendRunning;

extern volatile bool IsLocResourceAllocated;


// image info


extern CString ResultSavePath;
extern CString CreateTimeIdxStr;


extern cudaStream_t loc_stream1;

extern cudaStream_t render_stream1;

extern int SelectedGPUID;

extern tbb::concurrent_queue<qImgData> ImgDataQueue;


extern LocalizationPara LocPara_Global;

extern LDROIExtractData_TypeDef LDROIExtractData;

extern LDLocData_TypeDef LDLocData;

extern ZeroLocalizationsRemovel_TypeDef ZeroLocRemovel;


extern ConsecutiveFit_TypeDef ConsecutiveFitData;


extern FluoStatisticData_TypeDef FluoStatData;
extern ImageRenderData_TypeDef RendData;



extern NyqDimensionDensityCalc_TypeDef DimensionDensityCalc;
extern SpatialResolutionCalc_TypeDef SpatialResolutionCalc;

extern MarginDataFilter_TypeDef MarginDataFilter;
extern float FOVOverlapPercent;


void SetAllOnlineThreadAlive();
void ClearAllOnlineThreadAlive();


void InitAllLocResource();
void DeinitAllLocResource();

int SelectBestGPU();


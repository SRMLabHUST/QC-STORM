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
using namespace tbb;




struct qImgData
{
	unsigned short *pImgData;
	int ImageSource;
	int BatchFrameNum;

};



class ThreadCmdProcessState
{
private:
	volatile unsigned int CmdCnt; // number for asked commands
	volatile unsigned int ExeCnt; // number of executed commands

public:
	ThreadCmdProcessState()
	{
		ResetProcessState();
	}
	void ResetProcessState()
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
		int valid = (ExeCnt != CmdCnt);
		return valid;
	}
	void ProcessFinish()
	{
		ExeCnt = CmdCnt;
	}
};


extern volatile bool OnlineLocAlive;
extern volatile bool OnlineRendAlive;


extern volatile bool IsLocRunning;
extern volatile bool IsRendRunning;


extern volatile bool IsLocResourceAllocated;

extern float *h_RendFloatImage2D;


// image info
extern CString ImageName;


extern cudaStream_t loc_stream1;

extern cudaStream_t render_stream1;



extern concurrent_queue<qImgData> ImgDataQueue;


extern LocalizationPara LocPara_Global;
extern LDROIExtractData_TypeDef LDROIExtractData;

extern LDLocData_TypeDef LDLocData;

extern ZeroLocalizationsRemovel_TypeDef ZeroLocRemovel;


extern ConsecutiveFit_TypeDef ConsecutiveFitData;


extern FluoStatisticData_TypeDef FluoStatData;
extern ImageRenderData_TypeDef RendData;



extern cudaStream_t Resolution_stream1;

extern NyqDimensionDensityCalc_TypeDef DimensionDensityCalc;
extern SpatialResolutionCalc_TypeDef SpatialResolutionCalc;



extern int GPUID_1Best;
extern int GPUID_2Best;


void InitAllLocResource(int IsPostprocess);
void DeinitAllLocResource(int IsPostprocess);


void SelectBestGPU();


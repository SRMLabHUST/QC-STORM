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

#include "cuda_runtime.h"
#include "device_launch_parameters.h"


#include "DimensionLocDensityCalc_Para.h"

#include "ConsecutiveFilter.h"

#include "CurveFitting.h"

#include "cudaWrapper.h"

#include <iostream>
#include <vector>

using namespace std;



#define ThreadsPerBlock				32	//Threads Per Block




// calculate dimension and localization density, frame id is start from 1
class NyqDimensionDensityCalc_TypeDef
{
public:
	// for a group of localization data with dozens of frames

	float * h_LocArry;
	float * d_LocArry;


	// add a fraction of loc data, may be only several frames
	int AccumulatedFluoNum;
	int ImagesPerGroup_Valid; 
	
	
	float DimensionFD;
	float LocDensityFD; 

private:


	// fluo pos and number for each frame, frame id is start from 0
	vector<int> FrameStartPosArry;
	vector<int> FrameFluoNumArry;


	// valid frame number in a group to calculate calculate dimension and density
	int ValidAccuFrameVaryDataNum; 

	vector<float> AccumulatedFrameNumVary; // frame number
	vector<float> AccumulatedMeanDistanceVary; // unit is um

	
	// images number for a group to calculate dimension and localization density
	int ImagesPerGroup_Ideal;
	
	// consecutive filter for a group of localizations
	NyqConsecutiveFilter_TypeDef ConsecutiveFilter;


	// calculate min and mean nearest-neighbor distance
	float *h_ValidNum;
	float *h_TotalValue;
	float *d_ValidNum;
	float *d_TotalValue;

	float *d_MinDistance;

public:
	void Init(int MaxFluoNumPerGroup, int ImagesPerGroup);
	void DeInit();

	
	void ResetAccumulatedData();

	// two methods of add localization data: AddLocArray_FewFrames, AddLocArray_AGroup
	// only add localization data of few frames, when enough frames are accumulated, notification will be return
	int AddLocArray_FewFrames(float * h_iLocArry, int FluoNum, float ConsecFilter_DistanceTh_pixel, bool IsEndCalc, cudaStream_t cstream);
	
	// entire group of localizations are add, dimension and density can be calculated immediately
	void AddLocArray_AGroup(float * h_iLocArry, int FluoNum, float ConsecFilter_DistanceTh_pixel, cudaStream_t cstream);

	
	// get dimension, localization density fitting, use OffsetFrame as offset, calculate data point from StartFrame to EndFrame (not from the frist frame to EndFrame)
	void GetDimensionLocDensity_AGroup(int OffsetFrame, int StartFrame, int EndFrame, float PixelSize_nm, bool Is3DImaging, cudaStream_t cstream);

private:
	
	// calculate the mean Neighboring distance from(include) OffsetFrame to CurFrame
	float GetMeanMinNeighborDistance_AFrameGap(float PixelSize_nm, int OffsetFrame, int CurFrame, bool Is3DImaging, cudaStream_t cstream);


	void CalcMinNeighborDistance(float *d_LocArry, float PixelSize, int CalcStartPos, int SelDatLen, bool Is3DImaging, int FluoNum, cudaStream_t cstream);

	float CalcMeanDistance(float DistanceTh_pixel, int SelDatLen, cudaStream_t cstream);

	void GetDimensionDensity();


	void FindAllFramePos(float * h_LocArry, int TotalFluoNum);

	// can only be used after FindAllFramePos, frame id is start from 1
	int GetAccumulatedFluoNum(int OffsetFrame, int CurFrame);
	int GetFrameStartPos(int Frame);
	int GetFrameEndPos(int Frame);

public:
	static int GetFirstFrame(float * h_LocArry, int FluoNum);
	static int GetLastFrame(float * h_LocArry, int FluoNum);
	static int GetFluoNumPerFrame(float * h_LocArry, int FluoNum, int CurFrame, int StartPos);

};



__global__ void gpuMinDistanceCalc_2D(float *d_LocArry, float *d_MinDistance, float PixelSize_nm, int CalcStartPos, int SelDatLen, int FluoNum);
__global__ void gpuMinDistanceCalc_3D(float *d_LocArry, float *d_MinDistance, float PixelSize_nm, int CalcStartPos, int SelDatLen, int FluoNum);

__global__ void gpuMeanDistanceCalc(float *d_MinDistance, float *d_ValidNum, float *d_TotalValue, float DistanceTh_pixel, int FluoNum);


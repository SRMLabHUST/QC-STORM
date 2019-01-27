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

#include "bfgs_MLE_para.h"


#include <string>
#include <iostream>
#include <vector>
using namespace std;


#pragma comment(lib,"cudart.lib")
#pragma comment(lib,"SpatialResolutionCalc.lib")



// max fluo mum for a group with 50 frame 2048*2048 images
#define MAX_FLUO_NUM_PER_GROUP								(10240*3*50)

// max images in a group
#define MAX_FRAME_NUM_PER_GROUP								2000


// calculate only min neighboring distance of some molecules, calculation of all molecules is not necessary and time consuming
#define NEIGHBOR_DISTANCE_CALC_DATA_SELECT_NUMBER			22000





// filter molecules in consecutive frames by a radius threshold, only keep the molecule in the first frame
class NyqConsecutiveFilter_TypeDef
{
public:

	//
	int *d_ForwardLinkID;
	int *d_BackwardLinkID;
	int *d_MaxSNRID;

public:
	void Init(int TotalFluoNum); // create CPU&GPU memory
	void DeInit(); // release CPU&GPU memory

	void FilterConsecutiveFluo(float * d_LocArry, int FluoNum, float Distance_th_pixel, cudaStream_t cstream);

};





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





// get Nyquist and spatial resolution variation by each group's localization density and dimension
class SpatialResolutionCalc_TypeDef
{
public:

	// dimension and density of each group
	vector<int> FrameNum_Vary_Group;
	vector<int> AccumulatedFrameNum_Vary_Group;

	// dimension and localization density vary with group
	vector<float> Dimension_Vary_Group_FD; 
	vector<float> LocDensity_Vary_Group_FD; 
	vector<float> LocDensity_Vary_Group_2D; 
	vector<float> LocDensity_Vary_Group_3D; 

	// accumulated localization density
	vector<float> AccumulatedLocDensity_Vary_FD;
	vector<float> AccumulatedLocDensity_Vary_2D;
	vector<float> AccumulatedLocDensity_Vary_3D;


	// each data point is seperated by 10 frames
	vector<float> NyqResolutionVary_FD_10f; // mean Nyquist resolution of FD
	vector<float> NyqResolutionVary_2D_10f; // mean Nyquist resolution of 2D
	vector<float> NyqResolutionVary_3D_10f; // mean Nyquist resolution of 3D

	vector<float> NyquistResolutionVary_10f; // corrected mean Nyquist resolution
	vector<float> SpatialResolutionVary_10f; // lateral convolved spatial resolution

private:

	float StructureSize_2D; // nm

	// Nyquist and spatial resolution vary with frame
	int MetStructureSize2D_Pos; 
	int MetStructureSize3D_Pos; 
	float MetStructureSize2D_Value;
	float MetStructureSize3D_Value;

public:

	void Init();
	void DeInit();

	void SetStructureSize(float StructureSize_2D); // nm
	
	void AddDimensionLocDensity_AGroup(float DimensionFD, float LocDensityFD, int FrameNum_CurGroup);

	// calculate spatial resolution by Nyquist resolution and localization precision, each data is 10 frame increasing
	// MeanLocPrecXY is sqrt(loc prec X^2 + loc prec Y^2)
	// IsHollowTube : StructureSize_2D significantly larget than MeanLocPrec*2.35
	void GetSpatialResolutionVary(bool Is3DImaging, bool IsHollowTube, float MeanLocPrecXY, float MeanLocPrecZ, float OverSamplingRatio);
	
	
	float GetCurSpatialResolution();
	float GetCurNyquistResolution();

	int GetEstimatedAcquisitionFrame(float RSCResolution_th, float DimensionFD, float LocDensityFD, bool Is3DImaging, bool IsHollowTube, float MeanLocPrecXY, float MeanLocPrecZ, float OverSamplingRatio);

	void ResetData(); // already called in Init(), should be called if a new acquisition set will be calculated

	// debug
	void PrintData();

	float GetMeanVectorValue(vector<float> &iData_Vary_Group);

private:

	void GetDimensionConversionWidth(float *oWidth2D, float *oWidth3D, bool Is3DImaging, bool IsHollowTube, float MeanLocPrecXY, float MeanLocPrecZ);
	int GetAccumulatedFrameNum();

	int GetLocDensityGroupNum();

	float GetAccumulatedDensity(vector<float> &iLocDensity_Vary_Group, vector<int> &iFrameNum_Vary_Group);

	int GetGroupId_CurFrame(int CurFrame);

	float GetMeanDimensionFD(int CurGroup);


	int GetAccumulatedValue(vector<int> & iData_Vary);
	float GetAccumulatedValue(vector<float> & iData_Vary);

};





// use only center part of image, for multi-ROI stitching
class MarginDataFilter_TypeDef
{
public:
	float * h_LocArry;
	float * d_LocArry;

	float ValidFluoNum;


	void FilterMarginData(float* ih_LocArry, int FluoNum, float MarginPercentage, int ImageWidth, int ImageHigh, cudaStream_t cstream);

	void Init();
	void DeInit();

};



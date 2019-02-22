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

#include <iostream>
#include <vector>

using namespace std;


#include "DimensionLocDensityCalc_Para.h"


#define	SPATIAL_RESOLUTION_FRAME_GAP					10

#define NYQUIST_RESOLUTION_OVERSAMPLING					4



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


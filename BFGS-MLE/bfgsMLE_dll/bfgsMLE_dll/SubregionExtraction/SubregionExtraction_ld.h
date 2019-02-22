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

#include <stdio.h>
#include "cuda_runtime.h"
#include "device_launch_parameters.h"

#include "bfgs_CommonPara.h"

#include "cudaWrapper.h"

#include "bfgs_base.h"


#define ThreadsPerBlock				32 //Threads Per Block



#define BatchLineNum				8

#define ValidProcCore_7x7			(ThreadsPerBlock-6)
#define MarginSize_7x7				3

#define ValidProcCore_3x3			(ThreadsPerBlock-2)
#define MarginSize_3x3				1


#define Size32MB					(32 * 1024 * 1024)

#define MaxBatchedImageSize			Size32MB


#define RegionPosInfNum				2	// x,y



// molecular subregion detection and extraction for 2d and 3d low density
class LDROIExtractData_TypeDef
{
public:
	unsigned short *h_RawImg;
	unsigned short *d_RawImg;

	unsigned short *d_GaussImg;
	unsigned short *d_AnnuImg;
	unsigned short *d_StdImg;

	// extracted molecular ROI
	unsigned short * h_RegionMem;
	unsigned short * d_RegionMem;


	int *h_RegionNum; // region number array for batched frames
//	int *h_RegionPosMem;

	int h_TotalRegionCount;

	int *d_RegionNum; // region number array for batched frames
	int *d_RegionPosMem;

public:
	void Init(LocalizationPara & LocPara);
	void Deinit();

	// enable extract several images, the image memory must be continuously stored in the pImgData, will be very high efficient for small images, limited by MaxBatchedImageSize
	void ExtractMolecules(unsigned short *pImgData, int ImageSource, LocalizationPara & LocPara, int StartFrame, int BatchFrameNum, cudaStream_t cstream);

	int GetRegionNum();
	void ResetRegionNum();

};



 __global__ void gpuImgFilleringLD1(unsigned short *d_RawImg, unsigned short *d_AnnuImg, unsigned short *d_StdImg, int ImageWidth, int ImageHigh);
 __global__ void gpuImgFilleringLD2(unsigned short *d_GaussImg, unsigned short *d_AnnuImg, int ImageWidth, int ImageHigh);

 __global__ void gpuSubregionFindingLD(unsigned short *d_AnnuImg, unsigned short *d_StdImg, int *d_RegionPosMem, int *d_RegionNum, int RawImageWidth, int GroupImageHigh, int RawImgHigh);
 __global__ void gpuSubregionExtractionLD(unsigned short *d_RawImg, unsigned short *d_RegionMem, int *d_RegionPosMem, int FluoNum, int ROISize, int ImageWidth, int GroupImgHigh, int RawImgHigh, int StartFrame);


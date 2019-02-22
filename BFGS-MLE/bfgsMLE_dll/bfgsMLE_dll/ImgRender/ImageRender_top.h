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


#include "cudaWrapper.h"

#include "bfgs_CommonPara.h"

#include "bfgs_base.h"

#include "ImageRender_base.h"



// both 2d and 3d 
class ImageRenderData_TypeDef
{
public:
	float *h_LocArry;
	float *d_LocArry;

	float *d_SRIntensityImg; // gray-scale super-resolution image
	float *d_SRColorMapImg; // color weight for 3d color super-resolution image rendering

	int *h_MaxImageVal;
	int *d_MaxImageVal;

	int *h_HistMaxDat;
	int *d_HistMaxDat;

	char *h_DispRendImg; // croped and down-sampled super-resolution image for display
	char *d_DispRendImg; // croped and down-sampled super-resolution image for display

	char *h_SaveRendImg; // whole super-resolution image for save
	char *d_SaveRendImg; // whole super-resolution image for save

private:
	int tRendFluoNum;

public:

	void Init(LocalizationPara & LocPara, int MaxDispImgWidth, int MaxDispImgHigh);
	void Deinit(LocalizationPara & LocPara);

	// render super-resolution image by fill a gaussian ROI to each localizaed point
	void FluoRenderTop(float *h_LocArry, LocalizationPara & LocPara, int RenderingMode, float FixedlocPrec, int FluoNum, cudaStream_t cstream);

	 
	// get whole rendered image
	void GetDispImgTop(LocalizationPara & LocPara, float BrightRatio, int oImgWidth, int oImgHigh, int cposX, int cposY, float DispZoom, cudaStream_t cstream);
	void GetSaveImgTop(LocalizationPara & LocPara, float BrightRatio, int RGBImageEncodeMode, cudaStream_t cstream);
	void ResetFillImgTop(LocalizationPara & LocPara);

	int GetDispMaxVal();


	void ResetFillMaxVal(int Mode);

	// image render without knowing the image size, find it first
	static void GetMaxImgSizeFromLocArry(float *h_LocArry, float *d_LocArry, int *MaxImgWidth, int *MaxImgHigh, int FluoNum, cudaStream_t cstream);

};




__global__ void FindMaxImgSize(float *d_LocArry, int *d_MaxImageWidth, int *d_MaxImageHigh, int FluoNum, int FillParaNum);

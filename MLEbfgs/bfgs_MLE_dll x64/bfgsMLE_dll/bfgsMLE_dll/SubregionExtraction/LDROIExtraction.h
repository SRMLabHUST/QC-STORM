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

#include <stdio.h>


#include "bfgs_base.h"
#include "cudaWrapper.h"

#include "WLEParaEstimation.h"

#include "LDROIExtraction_Param.h"

#include <time.h>


// molecular subregion detection and extraction for 2d and 3d low density
class LDROIExtractData_TypeDef
{
public:
	// raw image
	unsigned short *h_RawImg;
	unsigned short *d_RawImg;

	// extracted molecular ROI data
	unsigned short * h_ImageROI;
	unsigned short * d_ImageROI;

	int FirstFrame;
	int EndFrame;

	WLEParameterEstimation_TypeDef *WLEParameterEstimator;

private:

	unsigned short * d_MoleculePosImage;

	// region number for batched frames
	int *h_ValidROINum;
	int *d_ValidROINum;

	
	int *d_ROIMarkInfArray; // molecule position

	int TotalROINumber;

private:
	// raw imaging filtering
	unsigned short *d_RawImg_Smoothed;
	unsigned short *d_BackgroundImage;

	unsigned short *d_LineFilterImage_t;
	unsigned short *d_LineFilterImage_t1;
	
	// image threshold calculate
	float *h_MeanDataX;
	float *h_MeanDataX2;
	float *d_MeanDataX;
	float *d_MeanDataX2;

	// threshold  = 10*sqrt(ImageVariance)
	float ImageVariance;

	// seperable image filter kernel
	float *h_LineFilterH_Signal;


public:
	void Init(LocalizationPara & LocPara);
	void Deinit();

	// BatchedImageNum is better an even number
	void ExtractMolecules(unsigned short *pImgData, int ImageSource, LocalizationPara & LocPara, int StartFrame_CurBatch, int BatchedImageNum, cudaStream_t cstream);

	int GetAccumulatedROINum();
	void ResetROINum();

	float * Get_h_WLEPara();
	float * Get_d_WLEPara();


public:
	static int GetMaxBatchedNumForCurrentImageSize(int ImageWidth, int ImageHigh);

private:

	void FilterInit(int ROISize, int LocType, int MultiEmitterFitEn);

	void ImageVarianceCalc(unsigned short *d_iRawImg, int ImageWidth, int ImageHigh, cudaStream_t cstream);

	void ImageFiltering(int LocType, int MultiEmitterFitEn, int ImageWidth, int ImageHigh, int BatchedImageNum, cudaStream_t cstream);

	void ROIExtraction(int ROISize, int LocType, int MultiEmitterFitEn, int ImageWidth, int ImageHigh, int BatchedImageNum, int StartFrame, cudaStream_t cstream);
};



// image filtering

// calculate background intensity
void BackgroundIntensityRemove(unsigned short *d_RawImg, unsigned short *d_BackgroundImage, int ImageWidth, int ImageHigh, cudaStream_t cstream);


// seperable convolution horizontal
void ImageConv_Seperable_H(unsigned short *d_iImage, unsigned short *d_oImage, int ImageWidth, int ImageHigh, cudaStream_t cstream);

// seperable convolution vertical
void ImageConv_Seperable_V(unsigned short *d_iImage, unsigned short *d_oImage, int ImageWidth, int ImageHigh, cudaStream_t cstream);

// ImageSubstract, d_ImageC = d_ImageA - d_ImageB
void ImageSubstract(unsigned short *d_ImageA, unsigned short *d_ImageB, unsigned short *d_ImageC, int ImageWidth, int ImageHigh, cudaStream_t cstream);

void ImageVarianceCalcOnSelectRegion(unsigned short *d_RawImg, float *d_MeanDataX, float *d_MeanDataX2, int ImageCenterX, int ImageCenterY, int StdImageSize, int ImageWidth, int ImageHigh, cudaStream_t cstream);

void ImageSharpen(unsigned short *d_RawImg, unsigned short *d_FilteredImg, int ImageWidth, int ImageHigh, cudaStream_t cstream);



// ROI extraction

void ROIFindingLD(unsigned short *d_RawImg_Smoothed, unsigned short * d_MoleculePosImage, int *d_ValidROINum, int *d_ROIMarkInfArray, float ImageVariance, int LocType, int MultiEmitterFitEn, int ROISize, int ImageWidth, int BatchedImageHigh, int ImageHigh, cudaStream_t cstream);


void DetectedROIClasify(unsigned short * d_MoleculePosImage, int *d_ROIMarkInfArray, int ROINumber, int MultiEmitterFitEn, int ROISize, int ImageWidth, cudaStream_t cstream);


void SubregionExtractionLD(unsigned short *d_RawImg, unsigned short *d_ImageROI, int *d_ROIMarkInfArray, int ROISize, int ROINumber, int ImageWidth, int BatchedImageHigh, int ImageHigh, int StartFrame, cudaStream_t cstream);




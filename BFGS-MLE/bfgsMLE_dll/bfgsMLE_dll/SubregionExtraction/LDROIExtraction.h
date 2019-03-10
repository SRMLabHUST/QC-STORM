#pragma once

#include "cuda_runtime.h"
#include "device_launch_parameters.h"

#include <stdio.h>


#include "bfgs_base.h"
#include "cudaWrapper.h"

#include "WLEParaEstimation.h"



#define ThreadsPerBlock				32 //Threads Per Block




#define MaxBatchedImageSize			(2048*2048*2*8)



#define ImageFilter_KernelSize		5


// gaurantee pixels in a thread block are in the same line
#define GetImageWidth_Thread(x)			((x + ThreadsPerBlock - 1) / ThreadsPerBlock * ThreadsPerBlock)

#define GetMaxFluoNumPerImage(x)		((int)(MaxPointNum / x))


#define RegionPosInfNum				2	// x,y



// molecular subregion detection and extraction for 2d and 3d low density
class LDROIExtractData_TypeDef
{
public:
	// raw image
	unsigned short *h_RawImg;
	unsigned short *d_RawImg;

	// extracted molecular ROI data
	unsigned short * h_ROIMem;
	unsigned short * d_ROIMem;

	WLEParameterEstimation_TypeDef *WLEParameterEstimator;

private:
	// region number for batched frames
	int *h_ROINumPerImage;
	int *d_ROINumPerImage;

	// region position array
	int *d_ROIPosArray;

	int TotalROINumber;

private:
	unsigned short *d_RawImg_Smoothed;
	unsigned short *d_BackgroundImage;
	unsigned short *d_LineFilterImage_t;
	
	// image threshold calculate
	float *h_MeanDataX;
	float *h_MeanDataX2;
	float *d_MeanDataX;
	float *d_MeanDataX2;

	// threshold  = 10*sqrt(ImageVariance)
	float ImageVariance;

	// image filter
	float *h_LineFilterH_Bkg;
	float *d_LineFilterH_Bkg;

	float *h_LineFilterH_Signal;
	float *d_LineFilterH_Signal;


	// consecutive fitting, to merge ROI by assign each consecutive ROI the average
	int *d_ForwardLinkID;
	int *d_BackwardLinkID;
	int *d_ConsecutiveNum;


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

	void FilterInit();

	void ImageVarianceCalc(unsigned short *d_iRawImg, int ImageWidth, int ImageHigh, cudaStream_t cstream);

	void ImageFiltering(int ImageWidth, int ImageHigh, int BatchedImageNum, cudaStream_t cstream);

	void ROIExtraction(int ROISize, int ImageWidth, int ImageHigh, int BatchedImageNum, int StartFrame, cudaStream_t cstream);
};



// image filtering

// calculate background intensity
void BackgroundIntensityCalc(unsigned short *d_RawImg, unsigned short *d_BackgroundImage, int ImageWidth, int ImageHigh, cudaStream_t cstream);


// seperable convolution horizontal
void ImageConv_Seperable_H(unsigned short *d_iImage, unsigned short *d_oImage, float *d_LineFilter, int ImageWidth, int ImageHigh, cudaStream_t cstream);

// seperable convolution vertical
void ImageConv_Seperable_V(unsigned short *d_iImage, unsigned short *d_oImage, float *d_LineFilter, int ImageWidth, int ImageHigh, cudaStream_t cstream);

// ImageSubstract, d_ImageC = d_ImageA - d_ImageB
void ImageSubstract(unsigned short *d_ImageA, unsigned short *d_ImageB, unsigned short *d_ImageC, int ImageWidth, int ImageHigh, cudaStream_t cstream);

void ImageVarianceCalcOnSelectRegion(unsigned short *d_RawImg, float *d_MeanDataX, float *d_MeanDataX2, int ImageCenterX, int ImageCenterY, int StdImageSize, int ImageWidth, int ImageHigh, cudaStream_t cstream);



// ROI extraction

void ROIFindingLD(unsigned short *d_RawImg_Smoothed, int *d_ROIPosArray, int *d_ROINumPerImage, float ImageVariance, int ROISize, int ImageWidth, int BatchedImageHigh, int ImageHigh, cudaStream_t cstream);


void SubregionExtractionLD(unsigned short *d_RawImg, unsigned short *d_ROIMem, int *d_ROIPosArray, int ROISize, int RegionNum_CurImage, int ImageWidth, int BatchedImageHigh, int ImageHigh, int StartFrame, cudaStream_t cstream);




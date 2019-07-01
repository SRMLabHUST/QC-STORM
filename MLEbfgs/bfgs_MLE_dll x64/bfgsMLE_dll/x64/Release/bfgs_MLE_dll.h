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

#ifndef __BFGS_MLE_GPU_H
#define __BFGS_MLE_GPU_H

#include "cuda_runtime.h"
#include "device_launch_parameters.h"

#include "bfgs_MLE_para.h"

#include <string>
#include <iostream>
#include <vector>
using namespace std;


#pragma comment(lib,"cudart.lib")
#pragma comment(lib,"bfgsMLE_dll.lib")




#ifdef BFGS_EXPORT
// export api
#define BFGS_MLE_API		__declspec(dllexport)

#else
// import api
// extern: static lib, __declspec(dllimport): DLL
#define BFGS_MLE_API		extern

#endif





// basic parameters for localization for both 2d and 3d
class LocalizationPara
{
public:
	// camera parameter
	float Offset; // DN
	float KAdc; // e-/DN
	float QE; // e-/pn
	float ReadNoise_e; // e-

	// localization and rendering
	int ROISize;
	int LocType;

	int MultiEmitterFitEn;

	int WLEEn;

	// consecutive molecules filter or fit
	float ConsecFit_DistanceTh_nm; // pixel
	int ConsecFitEn;

	// rendering
	float PixelSize; // raw image pixel size
	float PixelZoom; // the super-resolution image size is PixelSize/PixelZoom

	float SNR_th; // minimum SNR to render
	int ColorMode_3D;


	// raw image size in pixel
	int ImageWidth;
	int ImageHigh;

	int TotalFrameNum; // only valid for ImageJ plugin

	int SRImageWidth; // super-resolution image size in pixel
	int SRImageHigh;

	// 3d localization calibration curve para

	float MinZDepth; // min z depth of 3d imaging
	float MaxZDepth; // max z depth of 3d imaging

	float ZDepthCorrFactor;

	// calibration of sigma X >= sigma Y
	float p4_XGY;
	float p3_XGY;
	float p2_XGY;
	float p1_XGY;
	float p0_XGY;

	// calibration of sigma X < sigma Y
	float p4_XLY;
	float p3_XLY;
	float p2_XLY;
	float p1_XLY;
	float p0_XLY;


	// double-helix special
	int DH_RotateType;
	float DH_MeanDistance;
	float DH_DistanceTh;


	// spatial resolution calculation
	int ImagesPerGroup;
	int IsHollowTube; // tube width is significantly larger than localization precision

	float StrucuteSize_2D;

	// calculation control
	int SpatialResolutionCalcEn;


public:
	// construction
	LocalizationPara();

	void UpdateSRImageSize(); // SRImageWidth = ((int)(ImageWidth*PixelZoom) + 3) / 4 * 4; same for ImageHigh

	bool IsEqual(LocalizationPara & iPara);

};




struct CoreFittingPara
{
	int MultiEmitterFitEn;
	int WLEEn;

	// camera parameter
	float Offset; // DN
	float KAdc; // e-/DN
	float QE; // e-/pn
	float ReadNoise_e; // e-

	// 3D imaging
	float ZDepthCorrFactor;

	// calibration of sigma X >= sigma Y
	float p4_XGY;
	float p3_XGY;
	float p2_XGY;
	float p1_XGY;
	float p0_XGY;

	// calibration of sigma X < sigma Y
	float p4_XLY;
	float p3_XLY;
	float p2_XLY;
	float p1_XLY;
	float p0_XLY;

};


struct FitPosInf_TypeDef
{

	// single molecule fitting
	int * h_SingleFitFluoNum;
	int * d_SingleFitFluoNum; 
	int * d_SingleFitFluoPos;


	// two emitter fitting
	int * h_MultiFitFluoNum_2E;
	int * d_MultiFitFluoNum_2E;
	int * d_MultiFitFluoPos_2E;

	// three emitter fitting
	int * h_MultiFitFluoNum_3E;
	int * d_MultiFitFluoNum_3E;
	int * d_MultiFitFluoPos_3E;

	int * h_RejectedFluoNum;
	int * d_RejectedFluoNum;


	int *h_MultiFit_AddedFluoNum;
	int *d_MultiFit_AddedFluoNum;

};



// estimage WLE para, contained in the ROI extraction
class WLEParameterEstimation_TypeDef
{
public:

	float *h_WLEPara;
	float *d_WLEPara;

public:

	void Init(LocalizationPara & LocPara);
	void Deinit();

	void WLEParameterEstimate(unsigned short * d_ImageROI, int LocType, int MultiEmitterFitEn, int ROISize, int FluoNum, cudaStream_t cstream);
};




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








// both 2d and 3d localization data structure
class LDLocData_TypeDef
{
public:
	unsigned short * h_ImageROI;
	unsigned short * d_ImageROI;

	float *d_WLEPara; // copy data from roi extraction


	float * h_LocArry;
	float * d_LocArry;


	// valid number after localization
	int oValidFluoNum;


public:

	// ratio of ROI to all detected ROIs, include molecules fitted by multi modality
	float FitRatio_1E; // single molecule fitting ratio
	float FitRatio_2E; // two emitter fitting ratio
	float FitRatio_3E; // three emitter fitting ratio

	// ratio of ROI by last fitting modality to all detected ROIs
	float FitRatio_Final_1E; // single molecule fitting ratio
	float FitRatio_Final_2E; // two emitter fitting ratio
	float FitRatio_Final_3E; // three emitter fitting ratio
	float FitRatio_Final_4E; // four or more emitter fitting ratio

private:

	// for convinient parameter transimission
	CoreFittingPara* h_FitPara;
	CoreFittingPara* d_FitPara;

	FitPosInf_TypeDef* h_FitPosInf;
	FitPosInf_TypeDef* d_FitPosInf;



public:
	void Init(LocalizationPara & LocPara); // create CPU&GPU memory
	void Deinit(LocalizationPara & LocPara); // release CPU&GPU memory

	void BFGS_MLELocalization(unsigned short * h_ImageROI, float *h_WLEPara, LocalizationPara & LocPara, int FluoNum, cudaStream_t cstream);

private:

	void CopyFittingPara(LocalizationPara & LocPara, cudaStream_t cstream);


	void MoleculePreFitClasify(int ROISize, int MultiEmitterFitEn, int FluoNum, cudaStream_t cstream);

	void ResetNumbers(cudaStream_t cstream);

public:
	// note the frame must be sorted by ZeroLocalizationsRemove class
	static int GetFirstFrame(float * h_LocArry, int FluoNum);
	static int GetLastFrame(float * h_LocArry, int FluoNum);


	// two optional localization precision method, only suitable for 2d localization with symmetric Gaussian PSF
	static void LocPrecCalc_GaussianCRLB(float* d_LocArry, LocalizationPara & LocPara, int FluoNum, cudaStream_t cstream);

};







// consecutive fit to merge molecules in consecutive frames
class ConsecutiveFit_TypeDef
{
public:

	int OutFluoNum;
	float * h_OutLocArry;
	
	// on time calculate

	float *h_OntimeRatio;
private:

	int FluoNum_LastGroup;

	float * d_LocArry_LastGroup;
	float * d_LocArry_ConsecFit;


	// consecutive molecule find
	int *d_ForwardLinkID;
	int *d_BackwardLinkID;


	// on time calculate
	int *h_OntimeDistrib;
	int *d_OntimeDistrib;

	int *h_ValidFluoNum; // for ontime distribution 
	int *d_ValidFluoNum; // for ontime distribution 

public:
	void Init(); 
	void Deinit();

	void ConsecutiveFit_WeightedAvg(float * h_iLocArry, int FluoNum_CurGrop, int IsEndProc, LocalizationPara & LocPara, cudaStream_t cstream); // d_iLocArry come from localization data 

};





// remove invalid localizations with zero value
class ZeroLocalizationsRemovel_TypeDef {
public:
	int ValidFluoNum;

	float *h_LocArry;

public:
	// sort frame is necessary for multi emitter fitting that extra localizations are added on existing LocArry
	void RemoveZeroLocalizations(float *ih_LocArry,int iFluoNum, int SortFrameEn, int FirstFrame, int EndFrame, cudaStream_t cstream);

	void Init();

	void Deinit();


private:
	void RemoveZero_SortFrame(float *ih_LocArry, int iFluoNum, int FirstFrame, int EndFrame, cudaStream_t cstream);

	void RemoveZeroWithoutSortFrame(float *ih_LocArry, int iFluoNum, cudaStream_t cstream);
	void FindCopyID(float *ih_LocArry, int iFluoNum);

private:
	float *d_LocArry_Raw;
	float *d_LocArry_Valid;


	int *h_FluoID_Valid;
	int *d_FluoID_Valid;

	int *h_FluoNum_Valid;
	int *d_FluoNum_Valid;

};





// based on localized or even consecutive fitted position from bfgs_top
class DH3D_MoleculePair_TypeDef
{

public:
	float *h_LocArry;
	float *d_LocArry;

	float *h_oLocArry;
	float *d_oLocArry;

	int oValidFluoNum;


private:
	// paired molecule id for the same molecule
	int *h_PairID;
	int *d_PairID;

	int * h_ValidoFluoNum;
	int * d_ValidoFluoNum;

public:

	void Init();

	void Deinit();

	void MoleculePair(float *h_iLocArry, int FluoNum, LocalizationPara & LocPara, cudaStream_t cstream);

};




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
	
	char *h_SaveRendImg; // whole super-resolution image for save, mainly for 3D
	char *d_SaveRendImg; // whole super-resolution image for save, mainly for 3D


	// sr image histgram for 3d rendering threshold to correct save image
	float *h_SRImageHist;
	float *d_SRImageHist;

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




// stastical information class for both 2d and 3d
class FluoStatisticData_TypeDef
{
public:
	int FirstFrame;
	int EndFrame;

	// Distribution of current data
	int *h_Hist_TotalPhoton; 
	int *h_Hist_LocPrecisionXY;	
	int *h_Hist_PeakSNR;
	int *h_Hist_PSFWidth;

	// mean value of current data
	float MeanTotalPhoton;
	float MeanLocPrecisionXY; // sqrt(LocPrec_x^2+LocPrec_y^2)
	float MeanPeakSNR;
	float MeanBackground;
	float MeanPSFWidth;
	float MeanPSFWidth_Ctl; // h_MeanPSFWidth1 = MeanPSFWidth  * 10.0f / MeanPeakSNR1; // lower is better for feedback control
	float MeanLocDensity2D; // cur activation density

	// time variation curve data 

	vector<float> TimeVary_TotalPhoton;
	vector<float> TimeVary_LocPrecisionXY;
	vector<float> TimeVary_PeakSNR;
	vector<float> TimeVary_Background;
	vector<float> TimeVary_PSFWidth;
	vector<float> TimeVary_PSFWidth_Ctl;
	vector<float> TimeVary_LocDensity2D;


	// come from localization
	float *h_OntimeRatio; // for display and time vatiation
	float MeanOntime;


	vector<float> TimeVary_OntimeF1Ratio;
	vector<float> TimeVary_OntimeF2Ratio;
	vector<float> TimeVary_OntimeF3Ratio;
	vector<float> TimeVary_OntimeF4Ratio;


private:

	float *d_LocArry;

	// Distribution statistic info in gpu
	int *d_Hist_TotalPhoton; // total photon
	int *d_Hist_LocPrecisionXY; // theoretical precision
	int *d_Hist_PeakSNR; // peak PSF signal to PSF and background induced shot noise
	int *d_Hist_PSFWidth; // PSF width in pixel
	

	// get average SNR and PSF width
	float *h_ValidNum;
	float *h_LocPrecSum;
	float *h_SNRSum;
	float *h_PSFWSum;
	float *h_BgSum;

	float *d_ValidNum;
	float *d_LocPrecSum;
	float *d_SNRSum;
	float *d_PSFWSum;
	float *d_BgSum;

	// overlap molecular calculate to calculate localization density

	int *h_MaxFrame; // curremt max frame id

	int *h_FilteredFluoNum;
	int *h_FluoNum_n0; // 0 neighbor
	int *h_FluoNum_n1; // 1 neighbor

	int *d_FilteredFluoNum;
	float *d_FilteredLocArry;

	int *d_FluoNum_n0; // 0 neighbor
	int *d_FluoNum_n1; // 1 neighbor


public:

	void Init();
	void Deinit();

	// reset all distributions and time variation data
	void ResetAllDat(cudaStream_t cstream); // should be called befor first use GetStatisticalInf


	void GetStatisticalInf(float *h_LocArry, LocalizationPara & LocPara, int FluoNum, cudaStream_t cstream);

	void UpdateOntimeRatio(float *i_h_OntimeRatio);


public:
	static float GetActivationDensity(float Neibhgor0_Ratio, float Neibhgor1_Ratio, float RadiusTh_um);


	static float GetHistogramMeanData(int *HistData, int DatLen, float PercentTh);
	static int GetHistogramMaxData(int *HistData, int DatLen);
	static int GetHistogramMaxDataPos(int *HistData, int DatLen);
	static float GetHistogramWidth(int *HistData, int MaxPos, int DatLen);

	static float GetTimeVaryMean(vector<float>& iTimeVaryData);
	static float GetTimeVaryMax(vector<float>& iTimeVaryData);
	static float GetTimeVaryMin(vector<float>& iTimeVaryData);


private:
	// reset all distributions to avoid accumulate, keep time variation data,
	void ResetDistribDat(cudaStream_t cstream); // already be called in GetStatisticalInf()

	void GetOverlapMolecules(int*oNeighbor0_Num, int *oNeighbor1_Num, int *oTotalFluo, float RadiusTh_pixel, int CurFrame, LocalizationPara & LocPara, int FluoNum, cudaStream_t cstream);

	void UpdateStatDat(float *h_LocArry, int FluoNum);

	void CalcLocalizationDensity2D(float *h_LocArry, LocalizationPara & LocPara, int FluoNum, cudaStream_t cstream);

};




// drift correction by super-resoluction image cross-correlation 
class SRDriftCorrData_TypeDef
{

public:

	void Init(int RawImgWidth, int RawImgHigh, int MaxFrameNum);
	void Deinit();

	void CorrectSampleShift(string FileName, string oFileName, LocalizationPara & LocPara, int CorrFrameNum, cudaStream_t cstream);


public:
	static int GetTotalFrame(string FileName);
	static void GetMaxImgSize(string FileName, int *ImageWidth, int *ImageHigh);


private:
	int ImageWidth; // raw image width
	int ImageHigh; // raw image high
	int SRImageWidth; // raw image width
	int SRImageHigh; // raw image high


	int TotalFrame;

	int CorrGroupNum;

	int *GroupStartFrame; // start frame of each group
	int *GroupEndFrame; // start frame of each group

	int *GroupFrameStartPos; // start fluo position of each group
	int *GroupFrameEndPos; // start fluo position of each group


	float *XSliceShift;
	float *YSliceShift;

	float *XFrameShift;
	float *YFrameShift;

	float *d_XFrameShift;
	float *d_YFrameShift;

	// super-resolution image rendering
	float *h_LocArry; // render super-resolution image for corss-correlation
	float *d_LocArry; // render super-resolution image for corss-correlation

	float *h_FillImg1; // reference image
	float *h_FillImg2; // shifted image

	float *h_SumLine; // temporal use for correlation

	float *d_FillImg1; // reference image
	float *d_FillImg2; // shifted image

	float *d_MulImg;  // temporal use for correlation
	float *d_SumLine; // temporal use for correlation


private:

	// low level

	void CalcGroupNum(int TotalFrame, int CorrFrameNum);

	int GetAFrameEndFluoPos(string FileName, int OffsetFluoNum, int FindFrame);
	void GetCorrGroupFluoPos(string FileName);


	void ShiftInterpolation();
	void ApplyShiftTop(string FileName, string oFileName, int ShiftCorrEnable, cudaStream_t cstream);


	void RenderSlice(string FileName, int RendGroup, float PixelSize, float QE, float SNR_th, cudaStream_t cstream);

	void GetSliceShift(float *ShiftX, float *ShiftY, int CorrShiftBiasX, int CorrShiftBiasY, cudaStream_t cstream);

	void ResetFillImage(float *d_SRIntensityImg, int SRImageWidth, int SRImageHigh, cudaStream_t cstream);

	double CrossCorrelation(int ShiftX, int ShiftY, int CorrShiftBiasX, int CorrShiftBiasY, cudaStream_t cstream);

	// render 20480 molecules
	void ImageRender(float *h_LocArry, float *d_LocArry, float *d_SRIntensityImg, float QE, float SNR_th, float PixelSize, float PixelZoom, int SRImageWidth, int SRImageHigh, int FluoNum, cudaStream_t cstream);


	void gpuApplyShift(float *h_LocArry, int ShiftCorrEnable, int FluoNum, cudaStream_t cstream);


};



class ImageRender3DStackData_TypeDef
{
public:

	float *h_LocArry;
	float *d_LocArry;

	float *h_SRIntensityImg;
	float *d_SRIntensityImg;


public:


	void Init(int SRImageWidth, int SRImageHigh, int TotalFluoNum);
	void Deinit();

	// put all molecules data once for each z depth
	// FixedlocPrec_z is double of FixedlocPrec_x
	void GetSaveImgTop(float *ih_LocArry, LocalizationPara & LocPara, float RenderZDepth, int RGBImageEncodeMode, float FixedlocPrec_x, int TotalFluoNum, cudaStream_t cstream);

};



// wrapper for cuda function
BFGS_MLE_API char AllocHostMemory(void **ptr, long size);
BFGS_MLE_API char FreeHostMemory(void * ptr);
BFGS_MLE_API char AllocGPUMemory(void **ptr, long size);
BFGS_MLE_API char FreeGPUMemory(void * ptr);
BFGS_MLE_API char WaitGPUStream(cudaStream_t cstream);

BFGS_MLE_API char CreatStream(cudaStream_t*pstream);
BFGS_MLE_API char CreatStreamWithPriority(cudaStream_t *pstream, int prio);

BFGS_MLE_API char FreeStream(cudaStream_t cstream);

BFGS_MLE_API void HandleErr(cudaError_t err, const char * str );
BFGS_MLE_API int CudaDeviceNum();
BFGS_MLE_API char CudaSetDevice(int id);



#endif

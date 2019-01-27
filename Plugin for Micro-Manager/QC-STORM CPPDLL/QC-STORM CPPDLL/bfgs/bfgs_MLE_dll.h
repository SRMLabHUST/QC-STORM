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
	float ReadNoise_e;// e-

	// localization and rendering
	int ROISize;
	int LocType;
	int BadFitFilterWithAutoThEn; // filter fit result with automatic snr threshold: th = mean(SNR>4)/2

	int BackgroundNoiseFilterEn;

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

	float p4;
	float p3;
	float p2;
	float p1;
	float p0;


	// double-helix 3d loc para
	float MeanDistance;
	float DistanceTh;
	int RotateType;

	// spatial resolution calculation
	int ImagesPerGroup;
	int IsHollowTube; // tube width is significantly larger than localization precision
	float StrucuteSize_2D;
	float RSCResolutionTh; // expected resolution threshold

	// calculation control
	int OnTimeCalcEn;
	int SpatialResolutionCalcEn;


public:
	// construction
	LocalizationPara();

	void UpdateSRImageSize(); // SRImageWidth = ((int)(ImageWidth*PixelZoom) + 3) / 4 * 4; same for ImageHigh
};



// double-helix 3d pair data define
class DH3DPairData_TypeDef
{
public:
	float * h_PairedLocArry;
	float * d_PairedLocArry;

	float * h_PairedPosArry;
	float * d_PairedPosArry;

	int *h_ValidFluoNum;
	int *d_ValidFluoNum;

	int oValidFluoNum;

	int *h_DistanceDistrib;
	int *d_DistanceDistrib;


	float MeanDistance;
	float DistanceHistWidth;


public:
	void PairMolecules(float *d_LocArry, LocalizationPara & LocPara, int FluoNum, cudaStream_t cstream);

	void Init();
	void Deinit();


};







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





// both 2d and 3d localization data structure
class LDLocData_TypeDef
{
public:
	unsigned short * h_SubRegion;
	float * h_LocArry;

	unsigned short * d_SubRegion;
	float * d_LocArry;


	// Consecutive finding from adjecent frames
	// same with consecutive fitting method
	int *d_ForwardLinkID;
	int *d_BackwardLinkID;
	int *d_ConsecutiveNum;

	int *h_OntimeDistrib; 
	int *d_OntimeDistrib; 

	int *h_ValidFluoNum; // for ontime distribution 
	int *d_ValidFluoNum; // for ontime distribution 

	float *h_OntimeRatio; // for display and time vatiation

	// valid number after localization, still include filtered molecule number
	int oValidFluoNum;

	// double-helix 3d localization, just pair molecules after 2d localization

	DH3DPairData_TypeDef DH3DPairData;
	
private:
	// for loc filter
	float *h_SNRSumUp;
	int *h_ValidNum;

	float *d_SNRSumUp;
	int *d_ValidNum;


public:
	void Init(LocalizationPara & LocPara); // create CPU&GPU memory
	void Deinit(LocalizationPara & LocPara); // release CPU&GPU memory

	void BFGS_MLELocalization(unsigned short * h_SubRegion, LocalizationPara & LocPara, int FluoNum, cudaStream_t cstream);

	
	void OntimeCalc(LocalizationPara & LocPara, int FluoNum, cudaStream_t cstream);

	void CopyDataToGPU(float * h_LocArry, int FluoNum, cudaStream_t cstream);

public:

	static int GetFirstFrame(float * h_LocArry, int FluoNum);
	static int GetLastFrame(float * h_LocArry, int FluoNum);

	static int GetFirstFrameFromROI(unsigned short * h_SubRegion, int ROISize, int FluoNum);
	static int GetLastFrameFromROI(unsigned short * h_SubRegion, int ROISize, int FluoNum);


	// two optional localization precision method, only suitable for 2d localization with symmetric Gaussian PSF
	static void LocPrecCalc_GaussianCRLB(float* d_LocArry, LocalizationPara & LocPara, int FluoNum, cudaStream_t cstream);

private:
	void FilterBadFit(LocalizationPara & LocPara, int FluoNum, cudaStream_t cstream);

};







// filter background noised by automatic threshold
class BackgroundNoiseFilter_TypeDef
{
public:
	float * h_LocArry;
	float * d_LocArry;

private:
	// find distance threshold
	float *d_MinDistance;

	float *h_ValidNum;
	float *h_TotalValue;

	float *d_ValidNum;
	float *d_TotalValue;

	// filter noise based on neighbor numbe in threshold
	int* d_NeighborNum_Th1;
	int* d_NeighborNum_Th2;
	int* d_NeighborNum_Th3;

	int* d_IsNoise;

public:
	void Init(int TotalFluoNum); // create CPU&GPU memory
	void DeInit(); // release CPU&GPU memory

	// if the data source is cpu, data is copied out, otherwise, it's not
	void FilterBackgroundNoise(float * h_iLocArry, int FluoNum, int DataSource, cudaStream_t cstream);

private:
	void BackgroundNoiseRemove(int FluoNum, cudaStream_t cstream);
	void GetMinNearestNeighborDistance(int FluoNum, int SelFluoNum, cudaStream_t cstream);

	float GetMeanMinNearestNeighborDistance(int FluoNum, cudaStream_t cstream);
	void NeighborNumCalc(int FluoNum, cudaStream_t cstream);

};




// consecutive fit to group molecules scattered in consecutive frames into one
class ConsecutiveFit_TypeDef
{
public:

	int OutFluoNum;
	float * h_OutLocArry;


public:
	void Init(); // create CPU&GPU memory
	void Deinit(); // release CPU&GPU memory

	void FitConsecutiveFluo(float * d_iLocArry, LocalizationPara & LocPara, int FluoNum, cudaStream_t cstream, int IsLastProc); // d_iLocArry come from localization data 

	void GetAvailableData(cudaStream_t cstream); // get last data
	void GetResidualData(cudaStream_t cstream); // get last data+cur data,unsed at the last processing


	void ResetData();

private:

	int CurInputID;
	int BufFluoNum[2];

	float * d_BufLocArry[2];

	//
	int *d_ForwardLinkID;
	int *d_BackwardLinkID;


};



// remove invalid localizations with zero value
class ZeroLocalizationsRemovel_TypeDef {
public:
	int ValidFluoNum;

	float *h_LocArry;

public:
	void RemoveZeroLocalizations(float *ih_LocArry,int iFluoNum, cudaStream_t cstream);

	void Init();

	void Deinit();


private:
	void FindCopyID(float *ih_LocArry, int iFluoNum);

private:
	float *d_LocArry_Raw;
	float *d_LocArry_Valid;


	int *h_FluoID_Valid;
	int *d_FluoID_Valid;

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






// stastical information class for both 2d and 3d
class FluoStatisticData_TypeDef
{
public:

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
	static float GetActivationDensity(float Ov1MoleculesRatio, float RadiusTh_um);


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

	void GetOverlapMolecules(int*oNeighbor0_Num, int *oTotalFluo, float RadiusTh_pixel, int CurFrame, LocalizationPara & LocPara, int FluoNum, cudaStream_t cstream);

	void UpdateStatDat(float *h_LocArry, int FluoNum);

	void CalcLocalizationDensity2D(float *h_LocArry, LocalizationPara & LocPara, int FluoNum, cudaStream_t cstream);


};




// super-resoluction image cross-correlation drift correction

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
	void GetSaveImgTop(float *ih_LocArry, LocalizationPara & LocPara, float RenderZDepth, int RendMode, float FixedlocPrec_x, int TotalFluoNum, cudaStream_t cstream);



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

void HandleErr(cudaError_t err, const char * str);
int  CudaDeviceNum();
char CudaSetDevice(int id);




#endif

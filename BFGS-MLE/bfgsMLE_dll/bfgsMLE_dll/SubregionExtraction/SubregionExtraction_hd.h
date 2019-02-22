#ifndef __SUSREGION_EXTRACTION_HD_H
#define __SUSREGION_EXTRACTION_HD_H


#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include <stdio.h>

#include "bfgs_CommonPara.h"



#define ThreadsPerBlock			32 // Threads Per Block, fixed 32, should not be changed




class HDROIExtractData_TypeDef
{
public:
	unsigned short * h_RawImg; // raw images

	unsigned short * d_RawImg; // raw images

	unsigned short * d_GaussImg; // gaussian smooth filtered images
	unsigned short * d_AnnuImg; // annular filtered images

	unsigned short * d_StdImg1; // raw std filtered images
	unsigned short * d_StdImg2; // annular std filtered images

	unsigned short * h_PointImg; // raw point images
	unsigned short * h_cPointImg; // clustered point images

	unsigned short * d_PointImg; // raw point images

	int *h_iPosMem; // raw pos
	int *h_iPosCount; // raw pos
	int *h_oPosMem; // clustered pos
	int *h_oPosCount; // clustered pos

	int * d_RegionMem; // raw pos
	int * d_RegionNum; // raw pos


	unsigned short * h_RegionMem; // memory to store sub-regions
	int *h_TotalRegionNum; // total region number extracted
	
	unsigned short * d_RegionMem; // memory to store sub-regions

};



// high density 2d subregion detection and extraction

BFGS_MLE_API void HDROIExtract_ExtractMolecules(HDROIExtractData_TypeDef *h_HDROIExtractData, int LargeImgWidth, int LargeImgHigh, int ImgHighAdj, int RegionSize, int FrameOffset, cudaStream_t cstream);

BFGS_MLE_API void HDROIExtract_Init(HDROIExtractData_TypeDef **ha_HDROIExtractData, int LargeImgWidth, int LargeImgHigh, int RegionSize);
BFGS_MLE_API void HDROIExtract_Deinit(HDROIExtractData_TypeDef *h_HDROIExtractData);

BFGS_MLE_API int HDROIExtract_GetRegionNum(HDROIExtractData_TypeDef * h_HDROIExtractData);
BFGS_MLE_API void HDROIExtract_ResetRegionNum(HDROIExtractData_TypeDef * h_HDROIExtractData);


//
void cpuPointCluster(unsigned short * h_PointImg, int *h_iPosCount,int *h_iPosMem,unsigned short * h_cPointImg,int *h_oPosCount,int *h_oPosMem, int ImageWidth, int ImageHigh);

__global__ void gpuImgFiller_gauss_std(unsigned short * d_RawImg, unsigned short * d_GaussImg, unsigned short * d_StdImg1, int ImageWidth, int ImageHigh);
__global__ void gpuImgFiller_annu_stdannu(unsigned short * d_GaussImg, unsigned short * d_StdImg1, unsigned short * d_AnnuImg, unsigned short * d_StdImg2, int ImageWidth, int ImageHigh);
__global__ void gpuSubregionFinding(unsigned short * d_AnnuImg, unsigned short * d_StdImg2, unsigned short * d_PointImg, int * d_RegionNum, int * d_RegionMem, int ImageWidth, int ImageHigh);
__global__ void gpuSubregionExtraction(unsigned short *d_RawImg, int * d_RegionNum, int * d_RegionMem, unsigned short * d_RegionMem, int RegionSize, int LargeImgWidth, int LargeImgHigh, int ImgHighAdj, int FrameOffset);




#endif // __SUSREGION_EXTRACTION_HD_H


#ifndef __BFGS_MLE_HD_H
#define __BFGS_MLE_HD_H

#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include "math_functions.h"

#include "bfgs_hd_1.h"
#include "bfgs_hd_2.h"
#include "bfgs_hd_3.h"
#include "bfgs_hd_seperatedata.h"

#include <stdio.h>

#include "bfgs_CommonPara.h"



#define	 		PointParaNum	4

#define 	 	OutParaNum1 	4 
#define 	 	OutParaNum2 	8 
#define 	 	OutParaNum3 	12

#define			MaxParaNumHD	12


class HDLocData_TypeDef
{
public:
	unsigned short * h_SubRegion;
	float * h_LocArry;
	int *h_PointNum;

	unsigned short * d_SubRegion;
	float * d_LocArry1;
	float * d_LocArry2;
	float * d_LocArry3;

	int *d_Point1Num;
	int *d_Point1PosArry;
	int *d_Point2Num;
	int *d_Point2PosArry;
	int *d_Point3Num;
	int *d_Point3PosArry;

};


BFGS_MLE_API	int HDLoc_BFGS_MLELocalization(unsigned short * h_SubRegion, HDLocData_TypeDef* h_HDLocData, int FluoNum, float Offset, float kadc, float PsfWidth, int isFPGAProc, cudaStream_t cstream);

BFGS_MLE_API	int HDLoc_GetSubRegionOffset(unsigned short * FluoPoint, int ROISize);
BFGS_MLE_API	void HDLoc_GetImgSizeFromRegion(unsigned short * FluoPoint, int ROISize, int *ImageWidth, int *ImageHigh);
BFGS_MLE_API	int HDLoc_GetFirstFrameFromRegion(unsigned short * FluoPoint, int ROISize);

BFGS_MLE_API	void HDLoc_Init(HDLocData_TypeDef** ha_HDLocData, int ROISize);
BFGS_MLE_API	void HDLoc_Deinit(HDLocData_TypeDef* h_HDLocData);


#endif //__BFGS_MLE_HD_H

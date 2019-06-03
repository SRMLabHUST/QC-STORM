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


#include "bfgs_CommonPara.h"


#define bfgs_debug	0


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


	// spatial resolution calculation
	int ImagesPerGroup;
	int IsHollowTube; // tube width is significantly larger than localization precision
	float StrucuteSize_2D;
	float RSCResolutionTh; // expected resolution threshold

	// calculation control
	int SpatialResolutionCalcEn;


	int ProcessingMode; // online or offline

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


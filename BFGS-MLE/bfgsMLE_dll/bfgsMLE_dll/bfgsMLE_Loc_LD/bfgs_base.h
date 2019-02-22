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


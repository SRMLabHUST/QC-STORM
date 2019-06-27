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

#include "bfgs_base.h"

LocalizationPara::LocalizationPara()
{
	// camera parameter
	KAdc = 0.45;
	Offset = 100;
	QE = 0.72;
	ReadNoise_e = 1.3; // e-
	//

	ROISize = 7;
	LocType = LocType_GS2D;

	MultiEmitterFitEn = 0;

	WLEEn = 1;

	ConsecFit_DistanceTh_nm = 80.0f;
	ConsecFitEn = 0;

	// localization and rendering
	PixelSize = 100;
	PixelZoom = 5;

	SNR_th = 5;

	ColorMode_3D = 0;
	//
	ImageWidth = 2048;
	ImageHigh = 2048;

	TotalFrameNum = 1; // only valid for ImageJ plugin


	// 3d calibration curve para
	MinZDepth = -500;
	MaxZDepth = 500;

	ZDepthCorrFactor = 1.0f;

	// calibration of sigma X >= sigma Y
	p4_XGY = 0;
	p3_XGY = 0;
	p2_XGY = 0;
	p1_XGY = 1;
	p0_XGY = 0;

	// calibration of sigma X < sigma Y
	p4_XLY = 0;
	p3_XLY = 0;
	p2_XLY = 0;
	p1_XLY = 1;
	p0_XLY = 0;


	// spatial resolution
	ImagesPerGroup = 50;
	IsHollowTube = 0;
	StrucuteSize_2D = 40; // tublin

	// calculation control
	SpatialResolutionCalcEn = 1;


	UpdateSRImageSize();
}


void LocalizationPara::UpdateSRImageSize()
{
	// update render image size
	SRImageWidth = (int)(ImageWidth*PixelZoom);
	SRImageHigh = (int)(ImageHigh*PixelZoom);

	SRImageWidth = (SRImageWidth + 3) / 4 * 4;
	SRImageHigh = (SRImageHigh + 3) / 4 * 4;

}

bool LocalizationPara::IsEqual(LocalizationPara & iPara)
{
	int SumNum = 0;

	// camera parameter
	SumNum += this->Offset != iPara.Offset; // DN
	SumNum += this->KAdc != iPara.KAdc; // e-/DN
	SumNum += this->QE != iPara.QE; // e-/pn
	SumNum += this->ReadNoise_e != iPara.ReadNoise_e; // e-

													  // localization and rendering
	SumNum += this->ROISize != iPara.ROISize;
	SumNum += this->LocType != iPara.LocType;

	SumNum += this->MultiEmitterFitEn != iPara.MultiEmitterFitEn;

	SumNum += this->WLEEn != iPara.WLEEn;

	// consecutive molecules filter or fit
	SumNum += this->ConsecFit_DistanceTh_nm != iPara.ConsecFit_DistanceTh_nm; // pixel
	SumNum += this->ConsecFitEn != iPara.ConsecFitEn;

	// rendering
	SumNum += this->PixelSize != iPara.PixelSize; // raw image pixel size
	SumNum += this->PixelZoom != iPara.PixelZoom; // the super-resolution image size is PixelSize/PixelZoom

	SumNum += this->SNR_th != iPara.SNR_th; // minimum SNR to render
	SumNum += this->ColorMode_3D != iPara.ColorMode_3D;


	// raw image size in pixel
	SumNum += this->ImageWidth != iPara.ImageWidth;
	SumNum += this->ImageHigh != iPara.ImageHigh;

	SumNum += this->TotalFrameNum != iPara.TotalFrameNum; // only valid for ImageJ plugin

	SumNum += this->SRImageWidth != iPara.SRImageWidth; // super-resolution image size in pixel
	SumNum += this->SRImageHigh != iPara.SRImageHigh;

	// 3d localization calibration curve para

	SumNum += this->MinZDepth != iPara.MinZDepth; // min z depth of 3d imaging
	SumNum += this->MaxZDepth != iPara.MaxZDepth; // max z depth of 3d imaging

	SumNum += this->ZDepthCorrFactor != iPara.ZDepthCorrFactor;

	// calibration of sigma X >= sigma Y
	SumNum += this->p4_XGY != iPara.p4_XGY;
	SumNum += this->p3_XGY != iPara.p3_XGY;
	SumNum += this->p2_XGY != iPara.p2_XGY;
	SumNum += this->p1_XGY != iPara.p1_XGY;
	SumNum += this->p0_XGY != iPara.p0_XGY;

	// calibration of sigma X < sigma Y
	SumNum += this->p4_XLY != iPara.p4_XLY;
	SumNum += this->p3_XLY != iPara.p3_XLY;
	SumNum += this->p2_XLY != iPara.p2_XLY;
	SumNum += this->p1_XLY != iPara.p1_XLY;
	SumNum += this->p0_XLY != iPara.p0_XLY;


	// spatial resolution calculation
	SumNum += this->ImagesPerGroup != iPara.ImagesPerGroup;
	SumNum += this->IsHollowTube != iPara.IsHollowTube; // tube width is significantly larger than localization precision
	SumNum += this->StrucuteSize_2D != iPara.StrucuteSize_2D;

															  // calculation control
	SumNum += this->SpatialResolutionCalcEn != iPara.SpatialResolutionCalcEn;


	return SumNum == 0;
}

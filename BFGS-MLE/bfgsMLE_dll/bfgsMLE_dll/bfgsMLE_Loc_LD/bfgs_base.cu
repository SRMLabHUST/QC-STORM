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

	BadFitFilterWithAutoThEn = 1;

	BackgroundNoiseFilterEn = 1;

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

	ZDepthCorrFactor = 0.75f;

	p4 = 0;
	p3 = 0;
	p2 = 0;
	p1 = 1;
	p0 = 0;

	// double-helix 3d loc para
	MeanDistance = 10.1;
	DistanceTh = 1.0;
	RotateType = 0;

	// spatial resolution
	ImagesPerGroup = 50;
	IsHollowTube = 0;
	StrucuteSize_2D = 40; // tublin
	RSCResolutionTh = 0;

	// calculation control
	OnTimeCalcEn = 1;
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

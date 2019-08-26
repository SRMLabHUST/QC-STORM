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

#include "stdafx.h"
#include "StatInfDisplay.h"


// image display resources

volatile int StatDispSel = 0;

unsigned char *ConvertInfImageBuf;


StatisticalInformationDisplay StatInfDisplay;


void DisplayStatInf_SingleItem(char *DispStr, float DispValue,  CImg <unsigned char> *CImg_Axis, int XPos, int YPos)
{
	char tbuf[120];

	sprintf(tbuf, "%s: %.2f", DispStr, DispValue);
	CImg_Axis->draw_text(XPos, YPos, tbuf, CImg_White, 0, 1, CImgDisp_FontSize);

}

void DispStatInfOnFigure_Single(CImg <unsigned char> *CImg_Axis, int XPos, int YPos)
{


	// display together

	// mean value of current data
	DisplayStatInf_SingleItem("Current Total Photon", FluoStatData.MeanTotalPhoton, CImg_Axis, XPos, YPos + 0);
	DisplayStatInf_SingleItem("Current LocPrecisionXY (nm)", FluoStatData.MeanLocPrecisionXY, CImg_Axis, XPos, YPos + 20);
	DisplayStatInf_SingleItem("Current SNR", FluoStatData.MeanPeakSNR, CImg_Axis, XPos, YPos + 40);
	DisplayStatInf_SingleItem("Current Background", FluoStatData.MeanBackground, CImg_Axis, XPos, YPos + 60);
	DisplayStatInf_SingleItem("Current PSF Width (pixel)", FluoStatData.MeanPSFWidth, CImg_Axis, XPos, YPos + 80);
	DisplayStatInf_SingleItem("Current LocDensity 2D (0 neighbor)", FluoStatData.MeanLocDensity2D, CImg_Axis, XPos, YPos + 100);

}


void DispSpatialResolutionInf_SingleImage(CImg <unsigned char> *CImg_Axis, int XPos, int YPos)
{
	float CurSpatialResolution = SpatialResolutionCalc.CurSpatialResolution;
	float CurNyquistResolution = SpatialResolutionCalc.CurNyquistResolution;


	// display together
	DisplayStatInf_SingleItem("Convolved Resolution XY (nm)", CurSpatialResolution, CImg_Axis, XPos, YPos + 0);
	DisplayStatInf_SingleItem("Nyquist Resolution (nm)", CurNyquistResolution, CImg_Axis, XPos, YPos + 20);

}


void DispOntimeInf()
{

	float XPos = 80;
	float YPos = 360;

	auto OnTimeImage = StatInfDisplay.InfDisp_Curve_Ontime->GetAxisImage();

	DisplayStatInf_SingleItem("Estimated on time", FluoStatData.MeanOntime, OnTimeImage, XPos, YPos + 0);

	DisplayStatInf_SingleItem("Emission in 1 frame  ", FluoStatData.h_OntimeRatio[0], OnTimeImage, XPos, YPos + 20);
	DisplayStatInf_SingleItem("Emission in 2 frames ", FluoStatData.h_OntimeRatio[1], OnTimeImage, XPos, YPos + 40);
	DisplayStatInf_SingleItem("Emission in 3 frames ", FluoStatData.h_OntimeRatio[2], OnTimeImage, XPos, YPos + 60);
	DisplayStatInf_SingleItem("Emission in >3 frames", FluoStatData.h_OntimeRatio[3], OnTimeImage, XPos, YPos + 80);

}

void DispStatInfOnFigure_All()
{
	if (IsLocResourceAllocated == false)return;

	float XPos_Hist = 580;
	float YPos_Hist = 40;

	float XPos_TimeVary = 80;
	float YPos_TimeVary = 340;


	DispStatInfOnFigure_Single(StatInfDisplay.InfDisp_Hist_TotalPhoton->GetAxisImage(), XPos_Hist, YPos_Hist);
	DispStatInfOnFigure_Single(StatInfDisplay.InfDisp_Hist_LocPrec->GetAxisImage(), XPos_Hist, YPos_Hist);
	DispStatInfOnFigure_Single(StatInfDisplay.InfDisp_Hist_SNR->GetAxisImage(), XPos_Hist, YPos_Hist);
	DispStatInfOnFigure_Single(StatInfDisplay.InfDisp_Hist_PSFWidth->GetAxisImage(), XPos_Hist, YPos_Hist);


	DispStatInfOnFigure_Single(StatInfDisplay.InfDisp_Curve_TotalPhoton->GetAxisImage(), XPos_TimeVary, YPos_TimeVary);
	DispStatInfOnFigure_Single(StatInfDisplay.InfDisp_Curve_LocPrec->GetAxisImage(), XPos_TimeVary, YPos_TimeVary);
	DispStatInfOnFigure_Single(StatInfDisplay.InfDisp_Curve_SNR->GetAxisImage(), XPos_TimeVary, YPos_TimeVary);
	DispStatInfOnFigure_Single(StatInfDisplay.InfDisp_Curve_Background->GetAxisImage(), XPos_TimeVary, YPos_TimeVary);
	DispStatInfOnFigure_Single(StatInfDisplay.InfDisp_Curve_PSFWidth->GetAxisImage(), XPos_TimeVary, YPos_TimeVary);
	DispStatInfOnFigure_Single(StatInfDisplay.InfDisp_Curve_LocDensity2D->GetAxisImage(), XPos_TimeVary, YPos_TimeVary);


	DispOntimeInf();

	DispSpatialResolutionInf_SingleImage(StatInfDisplay.InfDisp_Curve_SpatialResolution->GetAxisImage(), XPos_Hist, YPos_Hist);
	DispSpatialResolutionInf_SingleImage(StatInfDisplay.InfDisp_Curve_NyquistResolution->GetAxisImage(), XPos_Hist, YPos_Hist);
	DispSpatialResolutionInf_SingleImage(StatInfDisplay.InfDisp_Curve_DimensionFD->GetAxisImage(), XPos_TimeVary, YPos_TimeVary);
	DispSpatialResolutionInf_SingleImage(StatInfDisplay.InfDisp_Curve_LocDensityFD->GetAxisImage(), XPos_TimeVary, YPos_TimeVary);

}



void UpdateStatInfDisplay()
{


	if (IsLocResourceAllocated == false)return;


	// histogram of total photon
	StatInfDisplay.InfDisp_Hist_TotalPhoton->SetAllData(FluoStatData.h_Hist_TotalPhoton, StatInf_Hist_DatLen, 0);
	// histogram of loc prec
	StatInfDisplay.InfDisp_Hist_LocPrec->SetAllData(FluoStatData.h_Hist_LocPrecisionXY, StatInf_Hist_DatLen, 0);
	// histogram of SNR
	StatInfDisplay.InfDisp_Hist_SNR->SetAllData(FluoStatData.h_Hist_PeakSNR, StatInf_Hist_DatLen, 0);
	// histogram of PSF width
	StatInfDisplay.InfDisp_Hist_PSFWidth->SetAllData(FluoStatData.h_Hist_PSFWidth, StatInf_Hist_DatLen, 0);


	// time variation curve of total photon
	StatInfDisplay.InfDisp_Curve_TotalPhoton->SetAllData(FluoStatData.TimeVary_TotalPhoton.data(), FluoStatData.TimeVary_TotalPhoton.size(), 0);
	// time variation curve of loc prec
	StatInfDisplay.InfDisp_Curve_LocPrec->SetAllData(FluoStatData.TimeVary_LocPrecisionXY.data(), FluoStatData.TimeVary_LocPrecisionXY.size(), 0);
	// time variation curve of SNR
	StatInfDisplay.InfDisp_Curve_SNR->SetAllData(FluoStatData.TimeVary_PeakSNR.data(), FluoStatData.TimeVary_PeakSNR.size(), 0);
	// time variation curve of Background
	StatInfDisplay.InfDisp_Curve_Background->SetAllData(FluoStatData.TimeVary_Background.data(), FluoStatData.TimeVary_Background.size(), 0);
	// time variation curve of PSF width
	StatInfDisplay.InfDisp_Curve_PSFWidth->SetAllData(FluoStatData.TimeVary_PSFWidth.data(), FluoStatData.TimeVary_PSFWidth.size(), 0);
	StatInfDisplay.InfDisp_Curve_PSFWidth->SetAllData(FluoStatData.TimeVary_PSFWidth_Ctl.data(), FluoStatData.TimeVary_PSFWidth_Ctl.size(), 1);

	// time variation curve of Localization density 2D
	StatInfDisplay.InfDisp_Curve_LocDensity2D->SetAllData(FluoStatData.TimeVary_LocDensity2D.data(), FluoStatData.TimeVary_LocDensity2D.size(), 0);



	// time variation curve of ontime
	StatInfDisplay.InfDisp_Curve_Ontime->SetAllData(FluoStatData.TimeVary_OntimeF1Ratio.data(), FluoStatData.TimeVary_OntimeF1Ratio.size(), 0);
	StatInfDisplay.InfDisp_Curve_Ontime->SetAllData(FluoStatData.TimeVary_OntimeF2Ratio.data(), FluoStatData.TimeVary_OntimeF2Ratio.size(), 1);
	StatInfDisplay.InfDisp_Curve_Ontime->SetAllData(FluoStatData.TimeVary_OntimeF3Ratio.data(), FluoStatData.TimeVary_OntimeF3Ratio.size(), 2);
	StatInfDisplay.InfDisp_Curve_Ontime->SetAllData(FluoStatData.TimeVary_OntimeF4Ratio.data(), FluoStatData.TimeVary_OntimeF4Ratio.size(), 3);


	// time variation curve of Dimension FD
	StatInfDisplay.InfDisp_Curve_DimensionFD->SetAllData(SpatialResolutionCalc.Dimension_Vary_Group_FD.data(), SpatialResolutionCalc.Dimension_Vary_Group_FD.size(), 0);
	// time variation curve of Localization density FD
	StatInfDisplay.InfDisp_Curve_LocDensityFD->SetAllData(SpatialResolutionCalc.LocDensity_Vary_Group_FD.data(), SpatialResolutionCalc.LocDensity_Vary_Group_FD.size(), 0);

	// time variation curve of Spatial resolution
	StatInfDisplay.InfDisp_Curve_SpatialResolution->SetAllData(SpatialResolutionCalc.NyquistResolutionVary_10f, SpatialResolutionCalc.VaryDataLen_10f, 0);
	// time variation curve of Nyquist resolution
	StatInfDisplay.InfDisp_Curve_NyquistResolution->SetAllData(SpatialResolutionCalc.SpatialResolutionVary_10f, SpatialResolutionCalc.VaryDataLen_10f, 0);



	StatInfDisplay.CreateFigure_All();


	DispStatInfOnFigure_All();

}


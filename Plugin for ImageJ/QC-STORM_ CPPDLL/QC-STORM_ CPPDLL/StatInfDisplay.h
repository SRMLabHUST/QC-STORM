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

#include "OnlineLocalizationLD.h"

#include "CImgCurveDisplay.h"



class StatisticalInformationDisplay {

public:

	// histogram of total photon
	CImgCurveDisplay <int, DISP_TYPE_HIST, LOG_AXIS_OFF, FIXED_DATA_LEN_ON, StatInf_Hist_DatLen, LABEL_DISPLAY_ON, LABEL_DISPLAY_ON, LABEL_TYPE_INT, 1> *InfDisp_Hist_TotalPhoton;
	// histogram of loc prec
	CImgCurveDisplay <int, DISP_TYPE_HIST, LOG_AXIS_OFF, FIXED_DATA_LEN_ON, StatInf_Hist_DatLen, LABEL_DISPLAY_ON, LABEL_DISPLAY_ON, LABEL_TYPE_INT, 1> *InfDisp_Hist_LocPrec;
	// histogram of SNR
	CImgCurveDisplay <int, DISP_TYPE_HIST, LOG_AXIS_OFF, FIXED_DATA_LEN_ON, StatInf_Hist_DatLen, LABEL_DISPLAY_ON, LABEL_DISPLAY_ON, LABEL_TYPE_INT, 1> *InfDisp_Hist_SNR;
	// histogram of PSF width
	CImgCurveDisplay <int, DISP_TYPE_HIST, LOG_AXIS_OFF, FIXED_DATA_LEN_ON, StatInf_Hist_DatLen, LABEL_DISPLAY_ON, LABEL_DISPLAY_ON, LABEL_TYPE_INT, 1> *InfDisp_Hist_PSFWidth;



	// time variation curve of total photon
	CImgCurveDisplay <float, DISP_TYPE_CURVE, LOG_AXIS_OFF, FIXED_DATA_LEN_OFF, 0, LABEL_DISPLAY_OFF, LABEL_DISPLAY_ON, LABEL_TYPE_FLOAT, 1> *InfDisp_Curve_TotalPhoton;
	// time variation curve of loc prec
	CImgCurveDisplay <float, DISP_TYPE_CURVE, LOG_AXIS_OFF, FIXED_DATA_LEN_OFF, 0, LABEL_DISPLAY_OFF, LABEL_DISPLAY_ON, LABEL_TYPE_FLOAT, 1> *InfDisp_Curve_LocPrec;
	// time variation curve of SNR
	CImgCurveDisplay <float, DISP_TYPE_CURVE, LOG_AXIS_OFF, FIXED_DATA_LEN_OFF, 0, LABEL_DISPLAY_OFF, LABEL_DISPLAY_ON, LABEL_TYPE_FLOAT, 1> *InfDisp_Curve_SNR;
	// time variation curve of background
	CImgCurveDisplay <float, DISP_TYPE_CURVE, LOG_AXIS_OFF, FIXED_DATA_LEN_OFF, 0, LABEL_DISPLAY_OFF, LABEL_DISPLAY_ON, LABEL_TYPE_FLOAT, 1> *InfDisp_Curve_Background;
	// time variation curve of PSF width
	CImgCurveDisplay <float, DISP_TYPE_CURVE, LOG_AXIS_OFF, FIXED_DATA_LEN_OFF, 0, LABEL_DISPLAY_OFF, LABEL_DISPLAY_ON, LABEL_TYPE_FLOAT, 2> *InfDisp_Curve_PSFWidth;
	// time variation curve of Localization density 2D
	CImgCurveDisplay <float, DISP_TYPE_CURVE, LOG_AXIS_OFF, FIXED_DATA_LEN_OFF, 0, LABEL_DISPLAY_OFF, LABEL_DISPLAY_ON, LABEL_TYPE_FLOAT, 1> *InfDisp_Curve_LocDensity2D;
	// time variation curve of ontime
	CImgCurveDisplay <float, DISP_TYPE_CURVE, LOG_AXIS_OFF, FIXED_DATA_LEN_OFF, 0, LABEL_DISPLAY_OFF, LABEL_DISPLAY_ON, LABEL_TYPE_FLOAT, 4> *InfDisp_Curve_Ontime;


	// time variation curve of dimension fd
	CImgCurveDisplay <float, DISP_TYPE_CURVE, LOG_AXIS_OFF, FIXED_DATA_LEN_OFF, 0, LABEL_DISPLAY_OFF, LABEL_DISPLAY_ON, LABEL_TYPE_FLOAT, 1> *InfDisp_Curve_DimensionFD;
	// time variation curve of localization density fd
	CImgCurveDisplay <float, DISP_TYPE_CURVE, LOG_AXIS_OFF, FIXED_DATA_LEN_OFF, 0, LABEL_DISPLAY_OFF, LABEL_DISPLAY_ON, LABEL_TYPE_FLOAT, 1> *InfDisp_Curve_LocDensityFD;
	// time variation curve of spatial resolution
	CImgCurveDisplay <float, DISP_TYPE_CURVE, LOG_AXIS_OFF, FIXED_DATA_LEN_OFF, 0, LABEL_DISPLAY_OFF, LABEL_DISPLAY_ON, LABEL_TYPE_FLOAT, 1> *InfDisp_Curve_SpatialResolution;
	// time variation curve of nyquist resolution
	CImgCurveDisplay <float, DISP_TYPE_CURVE, LOG_AXIS_OFF, FIXED_DATA_LEN_OFF, 0, LABEL_DISPLAY_OFF, LABEL_DISPLAY_ON, LABEL_TYPE_FLOAT, 1> *InfDisp_Curve_NyquistResolution;




	StatisticalInformationDisplay()
	{

		// histogram of total photon
		InfDisp_Hist_TotalPhoton = new CImgCurveDisplay <int, DISP_TYPE_HIST, LOG_AXIS_OFF, FIXED_DATA_LEN_ON, StatInf_Hist_DatLen, LABEL_DISPLAY_ON, LABEL_DISPLAY_ON, LABEL_TYPE_INT, 1>("Total photon distribution", false, 0, 0, false);
		// histogram of loc prec
		InfDisp_Hist_LocPrec = new CImgCurveDisplay <int, DISP_TYPE_HIST, LOG_AXIS_OFF, FIXED_DATA_LEN_ON, StatInf_Hist_DatLen, LABEL_DISPLAY_ON, LABEL_DISPLAY_ON, LABEL_TYPE_INT, 1>("Localization precision distribution", false, 0, 0, false);
		// histogram of SNR
		InfDisp_Hist_SNR = new CImgCurveDisplay <int, DISP_TYPE_HIST, LOG_AXIS_OFF, FIXED_DATA_LEN_ON, StatInf_Hist_DatLen, LABEL_DISPLAY_ON, LABEL_DISPLAY_ON, LABEL_TYPE_INT, 1>("SNR distribution", false, 0, 0, false);
		// histogram of PSF width
		InfDisp_Hist_PSFWidth = new CImgCurveDisplay <int, DISP_TYPE_HIST, LOG_AXIS_OFF, FIXED_DATA_LEN_ON, StatInf_Hist_DatLen, LABEL_DISPLAY_ON, LABEL_DISPLAY_ON, LABEL_TYPE_INT, 1>("PSF width distribution", false, 0, 0, false);


		// time variation curve of total photon
		InfDisp_Curve_TotalPhoton = new CImgCurveDisplay <float, DISP_TYPE_CURVE, LOG_AXIS_OFF, FIXED_DATA_LEN_OFF, 0, LABEL_DISPLAY_OFF, LABEL_DISPLAY_ON, LABEL_TYPE_FLOAT, 1>("Total photon variation", true, 0, 0, true);
		// time variation curve of loc prec
		InfDisp_Curve_LocPrec = new CImgCurveDisplay <float, DISP_TYPE_CURVE, LOG_AXIS_OFF, FIXED_DATA_LEN_OFF, 0, LABEL_DISPLAY_OFF, LABEL_DISPLAY_ON, LABEL_TYPE_FLOAT, 1>("Localization precision variation", true, 0, 0, true);
		// time variation curve of SNR
		InfDisp_Curve_SNR = new CImgCurveDisplay <float, DISP_TYPE_CURVE, LOG_AXIS_OFF, FIXED_DATA_LEN_OFF, 0, LABEL_DISPLAY_OFF, LABEL_DISPLAY_ON, LABEL_TYPE_FLOAT, 1>("SNR variation", true, 0, 0, true);
		// time variation curve of Background
		InfDisp_Curve_Background = new CImgCurveDisplay <float, DISP_TYPE_CURVE, LOG_AXIS_OFF, FIXED_DATA_LEN_OFF, 0, LABEL_DISPLAY_OFF, LABEL_DISPLAY_ON, LABEL_TYPE_FLOAT, 1>("Background variation", true, 0, 0, true);
		// time variation curve of PSF width
		InfDisp_Curve_PSFWidth = new CImgCurveDisplay <float, DISP_TYPE_CURVE, LOG_AXIS_OFF, FIXED_DATA_LEN_OFF, 0, LABEL_DISPLAY_OFF, LABEL_DISPLAY_ON, LABEL_TYPE_FLOAT, 2>("PSF width variation", true, 0, 0, true);
		// time variation curve of Localization density 2D
		InfDisp_Curve_LocDensity2D = new CImgCurveDisplay <float, DISP_TYPE_CURVE, LOG_AXIS_OFF, FIXED_DATA_LEN_OFF, 0, LABEL_DISPLAY_OFF, LABEL_DISPLAY_ON, LABEL_TYPE_FLOAT, 1>("Localization density 2D variation", true, 0, 0, true);
		// time variation curve of ontime
		InfDisp_Curve_Ontime = new CImgCurveDisplay <float, DISP_TYPE_CURVE, LOG_AXIS_OFF, FIXED_DATA_LEN_OFF, 0, LABEL_DISPLAY_OFF, LABEL_DISPLAY_ON, LABEL_TYPE_FLOAT, 4>("Ontime variation", true, 1, 0, true);


		// time variation curve of Dimension FD
		InfDisp_Curve_DimensionFD = new CImgCurveDisplay <float, DISP_TYPE_CURVE, LOG_AXIS_OFF, FIXED_DATA_LEN_OFF, 0, LABEL_DISPLAY_OFF, LABEL_DISPLAY_ON, LABEL_TYPE_FLOAT, 1>("Dimension FD variation", true, 0, 0, true);
		// time variation curve of Localization density FD
		InfDisp_Curve_LocDensityFD = new CImgCurveDisplay <float, DISP_TYPE_CURVE, LOG_AXIS_OFF, FIXED_DATA_LEN_OFF, 0, LABEL_DISPLAY_OFF, LABEL_DISPLAY_ON, LABEL_TYPE_FLOAT, 1>("Localization density FD variation", true, 0, 0, true);
		// time variation curve of Spatial resolution
		InfDisp_Curve_SpatialResolution = new CImgCurveDisplay <float, DISP_TYPE_CURVE, LOG_AXIS_OFF, FIXED_DATA_LEN_OFF, 0, LABEL_DISPLAY_OFF, LABEL_DISPLAY_ON, LABEL_TYPE_FLOAT, 1>("Spatial resolution variation", true, 0, 0, true);
		// time variation curve of Nyquist resolution
		InfDisp_Curve_NyquistResolution = new CImgCurveDisplay <float, DISP_TYPE_CURVE, LOG_AXIS_OFF, FIXED_DATA_LEN_OFF, 0, LABEL_DISPLAY_OFF, LABEL_DISPLAY_ON, LABEL_TYPE_FLOAT, 1>("Nyquist resolution variation", true, 0, 0, true);


		SetDataLabel();
	}
	
	void CreateFigure_All()
	{
		InfDisp_Hist_TotalPhoton->CreateFigure();
		InfDisp_Hist_LocPrec->CreateFigure();
		InfDisp_Hist_SNR->CreateFigure();
		InfDisp_Hist_PSFWidth->CreateFigure();


		InfDisp_Curve_TotalPhoton->CreateFigure();
		InfDisp_Curve_LocPrec->CreateFigure();
		InfDisp_Curve_SNR->CreateFigure();
		InfDisp_Curve_Background->CreateFigure();
		InfDisp_Curve_PSFWidth->CreateFigure();
		InfDisp_Curve_LocDensity2D->CreateFigure();
		InfDisp_Curve_Ontime->CreateFigure();


		InfDisp_Curve_DimensionFD->CreateFigure();
		InfDisp_Curve_LocDensityFD->CreateFigure();
		InfDisp_Curve_SpatialResolution->CreateFigure();
		InfDisp_Curve_NyquistResolution->CreateFigure();


	}

	void SetDataLabel()
	{

		InfDisp_Hist_TotalPhoton->SetLabelValue(0, Hist_MaxTotalPhoton, 0); // x
		InfDisp_Hist_TotalPhoton->SetLabelValue(0, 100, 1); // y

		InfDisp_Hist_LocPrec->SetLabelValue(0, Hist_MaxLocPrecision, 0); // x
		InfDisp_Hist_LocPrec->SetLabelValue(0, 100, 1); // y

		InfDisp_Hist_SNR->SetLabelValue(0, Hist_MaxPeakSNR, 0); // x
		InfDisp_Hist_SNR->SetLabelValue(0, 100, 1); // y

		InfDisp_Hist_PSFWidth->SetLabelValue(0, Hist_MaxPSFWidth, 0); // x
		InfDisp_Hist_PSFWidth->SetLabelValue(0, 100, 1); // y

	}
};


// image display resources

extern volatile int StatDispSel;
extern unsigned char *ConvertInfImageBuf;

extern StatisticalInformationDisplay StatInfDisplay;


void UpdateStatInfDisplay();




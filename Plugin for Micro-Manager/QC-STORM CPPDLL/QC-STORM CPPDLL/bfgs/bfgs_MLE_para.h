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


// these are already defiled constant parameters in the library
// to change them, first change them in relevant library and keep these parameters the same with the library


// subregion extraction

#define MaxBatchedImageSize			(2048*2048*2*4)




// localziation type, for 2d or 3d selection
#define LocType_GS2D				0	// Gaussian symitrical 2d localization
#define LocType_AS3D				1	// astigmatism 3d


// use gaussian 2d or 3d localization
#define LocType_IsGS2D(x)				(x == LocType_GS2D)
#define LocType_IsAS3D(x)				(x == LocType_AS3D)

// use red hot 2d or color encoded depth 3d
#define RendType_Is2D(x)				(x == LocType_GS2D)
#define RendType_Is3D(x)				(x == LocType_AS3D)



// stored parameter number for each types of fitting for each fitted molecule

#define OutParaNumGS2D			12 // gaussian 2d
#define OutParaNumAS3D			12 // astigmatism 3d


// loc results information order for OutParaNumGS2D and OutParaNumAS3D and OutParaNumDH3D

#define Pos_PPho			0 // peak photon
#define Pos_XPos			1 // may have 0.5 or 1 pixel offset compared with other software
#define Pos_YPos			2 // may have 0.5 or 1 pixel offset compared with other software
#define Pos_ZPos			3 // may have 0.5 or 1 pixel offset compared with other software
#define Pos_SigX			4 // sigma x
#define Pos_SigY			5 // sigma y
#define Pos_TPho			6 // total photon
#define Pos_Bakg			7 // background
#define Pos_PSNR			8 // peak photon to background noise snr
#define Pos_CrbX			9 // crlb of x  
#define Pos_CrbY			10 // crlb of Y
#define Pos_Frme			11 // frame




// maximum molecules in allocated memory for each class
#define PointNumTh					20480
#define MaxPointNum					(20480*5)

// memory type
#define ImageSource_GPU					0
#define ImageSource_CPU_Pinned			1
#define ImageSource_CPU_Normal			2
#define ImageSource_GPU_NoCopy			3
#define ImageSource_ERR					1000


#define ProcessMode_Online				0
#define ProcessMode_Offline				1



// stastical information

#define StatInf_Hist_DatLen				400
#define StatInf_Hist_DatLenMax			(StatInf_Hist_DatLen*4)



#define Hist_MaxTotalPhoton				10000.0f


#define Hist_MaxLocPrecision			60.0f


#define Hist_MaxPeakSNR					50.0f


#define Hist_MaxPSFWidth				5.0f



#define LocDensityCalc_RadiusTh_pixel_2D		10
#define LocDensityCalc_RadiusTh_pixel_3D		20




/*
RenderingMode 0 : peak photon as weight for each molecule, rendered by localization precision calculated by CRLB
RenderingMode 1 : 1 as weight for each molecule, rendered by localization precision calculated by CRLB
RenderingMode 2 : 1 as weight for each molecule, rendered by fixed localization precision
*/

#define RenderingMode_FittedPhoton_CalculatedLocPrec		0
#define RenderingMode_1Photon_CalculatedLocPrec				1
#define RenderingMode_1Photon_FixedLocPrec					2


// regular RGB image
#define RGBImage_EncodeMode_3Bytes			0
// for ImageJ display
#define RGBImage_EncodeMode_4Bytes			1

// blue is low depth or reg is low depth
#define ImgRend_ColorMode_3D_BlueToRed		0
#define ImgRend_ColorMode_3D_RedToBlue		1





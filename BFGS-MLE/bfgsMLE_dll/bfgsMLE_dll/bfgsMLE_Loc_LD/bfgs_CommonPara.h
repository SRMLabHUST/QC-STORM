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


// this should only be defined in dll lib project
#define BFGS_EXPORT				1

// api declaration

#ifdef BFGS_EXPORT
// export api
#define BFGS_MLE_API		__declspec(dllexport)

#else
// import api
// extern: static lib,	__declspec(dllimport): DLL
#define BFGS_MLE_API		extern

#endif




// localziation type, for 2d or 3d selection
#define LocType_GS2D				0	// Gaussian symitrical 2d localization
#define LocType_AS3D				1	// astigmatism 3d
#define LocType_DH3D				2	// double-helix 3d


// use gaussian 2d or 3d localization
#define LocType_IsGS2D(x)				((x == LocType_GS2D)||(x == LocType_DH3D))
#define LocType_IsAS3D(x)				((x == LocType_AS3D))

// use red hot 2d or color encoded depth 3d
#define RendType_Is2D(x)				((x == LocType_GS2D))
#define RendType_Is3D(x)				((x == LocType_AS3D)||(x == LocType_DH3D))



// stored parameter number for each types of fitting for each fitted molecule

#define OutParaNumGS2D			12 // gaussian 2d
#define OutParaNumAS3D			12 // astigmatism 3d
#define OutParaNumDH3D			12 // double-helix 3d


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
#define PointNumTh				20480
#define MaxPointNum				(20480*4)

// memory type
#define ImageSource_GPU					0
#define ImageSource_CPU_Pinned			1
#define ImageSource_CPU_Normal			2
#define ImageSource_GPU_NoCopy			3
#define ImageSource_ERR					1000



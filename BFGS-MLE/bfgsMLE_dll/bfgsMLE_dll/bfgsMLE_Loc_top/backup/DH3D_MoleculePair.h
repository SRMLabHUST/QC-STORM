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

//#pragma once
#include "cuda_runtime.h"
#include "device_launch_parameters.h"




#include <stdio.h>

#include "bfgs_CommonPara.h"

#include "bfgs_base.h"

#include "statisticsInfo.h"


#define ThreadsPerBlock			32 //Threads Per Block


// for OutParaNumDH3D paired molecules information output

#define Pos_PPho0			0 // peak photon
#define Pos_XPos0			1 // may have 0.5 or 1 pixel offset compared with other software
#define Pos_YPos0			2 // may have 0.5 or 1 pixel offset compared with other software
#define Pos_PPho1			3 // may have 0.5 or 1 pixel offset compared with other software
#define Pos_XPos1			4 // molecule 1 x pos
#define Pos_YPos1			5 // molecule 1 y pos
#define Pos_ZPos0			6 // molecule 2 x pos
#define Pos_Dista			7 // molecule 2 y pos
#define Pos_Angle			8 // molecules pair angle
#define Pos_Frame			9 // frame



#define DistanceHistLen			300
#define MaxPixelDistance		30.0f

#define DistanceGap				(MaxPixelDistance/DistanceHistLen)








__global__ void gpuPairMolecules(float *d_LocArry, float *d_oLocArry, float *d_oPosArry, int *d_ValidFluoNum, float MeanDistance, float DistanceTh, int FluoNum, int PosArryEnable, float p4, float p3, float p2, float p1, float p0, int RotateType);
__global__ void gpuCalcDistanceDistribution(float *d_PosArry, int *d_DistanceDistrib, int FluoNum);



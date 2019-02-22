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

#include "cuda_runtime.h"
#include "device_launch_parameters.h"




#include <stdio.h>

#include "bfgs_CommonPara.h"



#define OverlapPairDistTh				3.0f

#define ThreadsPerBlock					32	//Threads Per Block


#define MaxOnTimeConsecutiveNum			10

// pair molecules in different frames and loc them as a single molecule


__global__ void gpuOntimeCalcFluoPair(float *d_LocArry, int *d_ForwardLinkID, int *d_BackwardLinkID, float PixelSize, float QE, int FluoNum);
__global__ void gpuOntimeCalc(float *d_LocArry, int *d_ForwardLinkID, int *d_BackwardLinkID, int *d_ConsecutiveNum, float PixelSize, float QE, int FluoNum);
__global__ void gpuGetOntimeDistribution(float *d_LocArry, int *d_ConsecutiveNum, int *d_OntimeDistrib, int *d_ValidFluoNum, int StartFrame, int EndFrame, int FluoNum);

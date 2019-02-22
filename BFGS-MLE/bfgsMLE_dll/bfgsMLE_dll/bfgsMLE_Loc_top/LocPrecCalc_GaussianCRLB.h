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


#define ThreadsPerBlock						32


#define GaussianCRLB_ParaNum				6

#define CRLB_Para_TPhn						0
#define CRLB_Para_Bkgn						1
#define CRLB_Para_XPos						2
#define CRLB_Para_YPos						3
#define CRLB_Para_SigX						4
#define CRLB_Para_SigY						5

#define GaussianCRLB_ROISize2D				9
#define GaussianCRLB_ROISize3D				15


#define TakeReadNoiseToBackground			1


__global__ void GaussianCRLB_Calc_top(float *d_LocArry, float ReadNoise_e, float QE, float PixelSize);


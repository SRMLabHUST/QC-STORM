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

#include "cuda_runtime.h"
#include "device_launch_parameters.h"

#include <stdio.h>


#include "bfgs_base.h"
#include "cudaWrapper.h"



#define ThreadsPerBlock				32 //Threads Per Block




#define MaxBatchedImageSize			(2048*2048*2*8)



// gaurantee pixels in a thread block are in the same line
#define GetImageWidth_Thread(x)			((x + ThreadsPerBlock - 1) / ThreadsPerBlock * ThreadsPerBlock)



// mark a possible molecule
// x, y, molecule type, Nearest neighbor distance valid, Nearest neighbor distance possible, null
#define ROIMarkInfNum					4


#define ROIMarkInf_XPos					0
#define ROIMarkInf_YPos					1
#define ROIMarkInf_TypeFit				2


#define ROIType_Fit_Single				0	
#define ROIType_Fit_Multi				1



// image filtering

// calculate background intensity
#define BlockFilter_BatchLineNum		8


#define BkgFilterSize					9			
#define BkgFilterSize_Half				((int)(BkgFilterSize/2))


// seperable convolution
#define LineFilterSize					5


#define SharpFilterSize					3		
#define SharpFilterSize_Half			((int)(SharpFilterSize/2))


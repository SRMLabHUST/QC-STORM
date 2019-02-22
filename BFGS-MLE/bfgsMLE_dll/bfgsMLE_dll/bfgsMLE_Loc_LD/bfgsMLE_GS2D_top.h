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

#include "bfgs_base.h"


#include "bfgsMLE_GS2D_core.h"


// bfgs 2d loc function
void LDLoc_BFGS_MLELocalizationGS2D(unsigned short * d_SubRegion, float * d_LocArry, LocalizationPara& LocPara, int FluoNum, cudaStream_t cstream);







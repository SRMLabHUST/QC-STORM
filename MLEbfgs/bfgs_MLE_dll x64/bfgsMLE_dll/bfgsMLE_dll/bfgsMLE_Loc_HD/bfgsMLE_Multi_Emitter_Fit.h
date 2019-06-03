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

#include "bfgs_base.h"





void HDLoc_BFGS_MLELocalization_2D_2Emitter(float * d_LocArry, unsigned short * d_ImageROI, float *d_WLEPara, int MultiFitFluoNum_2E, FitPosInf_TypeDef* d_FitPosInf, CoreFittingPara *d_FitPara, int ROISize, cudaStream_t cstream);


void HDLoc_BFGS_MLELocalization_2D_3Emitter(float * d_LocArry, unsigned short * d_ImageROI, float *d_WLEPara, int MultiFitFluoNum_3E, FitPosInf_TypeDef* d_FitPosInf, CoreFittingPara *d_FitPara, int ROISize, cudaStream_t cstream);



void HDLoc_BFGS_MLELocalization_AS3D_2Emitter(float * d_LocArry, unsigned short * d_ImageROI, float *d_WLEPara, int MultiFitFluoNum_2E, FitPosInf_TypeDef* d_FitPosInf, CoreFittingPara *d_FitPara, int ROISize, cudaStream_t cstream);

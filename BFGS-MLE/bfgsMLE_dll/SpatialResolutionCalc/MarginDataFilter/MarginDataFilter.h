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


#include "bfgs_CommonPara.h"


#define ThreadsPerBlock								32		//Threads Per Block

// use only center part of image, for multi-ROI stitching
class MarginDataFilter_TypeDef
{
public:
	float * h_LocArry;
	float * d_LocArry;

	float ValidFluoNum;


	void FilterMarginData(float* ih_LocArry, int FluoNum, float MarginPercentage, int ImageWidth, int ImageHigh, cudaStream_t cstream);

	void Init();
	void DeInit();

};


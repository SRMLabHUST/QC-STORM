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



#include "DimensionLocDensityCalc_Para.h"


#define ThreadsPerBlock				32	//Threads Per Block



// filter molecules in consecutive frames by a radius threshold, only keep the molecule in the first frame
class NyqConsecutiveFilter_TypeDef
{
public:

	//
	int *d_ForwardLinkID;
	int *d_BackwardLinkID;
	int *d_MaxSNRID;

public:
	void Init(int TotalFluoNum); // create CPU&GPU memory
	void DeInit(); // release CPU&GPU memory

	void FilterConsecutiveFluo(float * d_LocArry, int FluoNum, float Distance_th_pixel, cudaStream_t cstream);

};


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

#include <iostream>
using namespace std;


#include "cuda_runtime.h"
#include "device_launch_parameters.h"


#include "bfgs_CommonPara.h"

#define ThreadsPerBlock					32 //Threads Per Block




// remove invalid localizations with zero value
class ZeroLocalizationsRemovel_TypeDef {
public:
	int ValidFluoNum;

	float *h_LocArry;

public:
	// sort frame is necessary for multi emitter fitting that extra localizations are added on existing LocArry
	void RemoveZeroLocalizations(float *ih_LocArry,int iFluoNum, int SortFrameEn, int FirstFrame, int EndFrame, cudaStream_t cstream);

	void Init();

	void Deinit();


private:
	void RemoveZero_SortFrame(float *ih_LocArry, int iFluoNum, int FirstFrame, int EndFrame, cudaStream_t cstream);

	void RemoveZeroWithoutSortFrame(float *ih_LocArry, int iFluoNum, cudaStream_t cstream);
	void FindCopyID(float *ih_LocArry, int iFluoNum);

private:
	float *d_LocArry_Raw;
	float *d_LocArry_Valid;


	int *h_FluoID_Valid;
	int *d_FluoID_Valid;

	int *h_FluoNum_Valid;
	int *d_FluoNum_Valid;

};



__global__ void gpuRemoveZeorLocalizations(float * d_LocArry_Valid, float * d_LocArry_Raw, int *d_FluoID_Valid, int ValidFluoNum);

__global__ void gpuRemoveZeorLoc_SortFrame(float * d_LocArry_Valid, float * d_LocArry_Raw, int *d_FluoNum_Valid, int CurFrame, int FluoNum);


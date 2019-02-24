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
#include "bfgs_base.h"


#define ThreadsPerBlock				32	//Threads Per Block


#define MAX_CONSEC_NUM				8

// consecutive fit to group molecules scattered in consecutive frames into one
class ConsecutiveFit_TypeDef
{
public:

	int OutFluoNum;
	float * h_OutLocArry;


public:
	void Init(); // create CPU&GPU memory
	void Deinit(); // release CPU&GPU memory

	void FitConsecutiveFluo(float * d_iLocArry, LocalizationPara & LocPara, int FluoNum, cudaStream_t cstream, int IsLastProc); // d_iLocArry come from localization data 

	void GetAvailableData(cudaStream_t cstream); // get last data
	void GetResidualData(cudaStream_t cstream); // get last data+cur data,unsed at the last processing


	void ResetData();

private:

	int CurInputID;
	int BufFluoNum[2];

	float * d_BufLocArry[2];

	//
	int *d_ForwardLinkID;
	int *d_BackwardLinkID;


};



__global__ void gpuFindConsecutiveMolecules(float *d_LastLocArry, int LastFluoNum, float *d_CurLocArry, int CurFluoNum, int *d_ForwardLinkID, int *d_BackwardLinkID, float DistanceTh_nm, float PixelSize, float QE);

__global__ void gpuConsecutiveFit(float *d_LastLocArry, int LastFluoNum, float *d_CurLocArry, int CurFluoNum, int *d_ForwardLinkID, int *d_BackwardLinkID, int ConsecFitFluoNum, float QE);

__global__ void gpuRemoveConsecutiveFluo(float *d_LastLocArry, int LastFluoNum, float *d_CurLocArry, int CurFluoNum, int *d_ForwardLinkID, int *d_BackwardLinkID, int ConsecFitFluoNum);


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

#include "statisticsInfo.h"



#define ThreadsPerBlock				32	//Threads Per Block




// consecutive fit to merge molecules in consecutive frames
class ConsecutiveFit_TypeDef
{
public:

	int OutFluoNum;
	float * h_OutLocArry;
	
	// on time calculate

	float *h_OntimeRatio;
private:

	int FluoNum_LastGroup;

	float * d_LocArry_LastGroup;
	float * d_LocArry_ConsecFit;


	// consecutive molecule find
	int *d_ForwardLinkID;
	int *d_BackwardLinkID;


	// on time calculate
	int *h_OntimeDistrib;
	int *d_OntimeDistrib;

	int *h_ValidFluoNum; // for ontime distribution 
	int *d_ValidFluoNum; // for ontime distribution 

public:
	void Init(); 
	void Deinit();

	void ConsecutiveFit_WeightedAvg(float * h_iLocArry, int FluoNum_CurGrop, int IsEndProc, LocalizationPara & LocPara, cudaStream_t cstream); // d_iLocArry come from localization data 

};



__global__ void gpuFindConsecutiveMolecules(float * d_LocArry_ConsecFit, int TotalFluoNum, int *d_ForwardLinkID, int *d_BackwardLinkID, float DistanceTh_Pixel);

__global__ void gpuConsecutiveFit(float * d_LocArry_ConsecFit, int TotalFluoNum, int FluoNum_LastGroup, int *d_ForwardLinkID, int *d_BackwardLinkID, int *d_OntimeDistrib, int *d_ValidFluoNum, float QE);
__global__ void gpuRemoveConsecutiveFluo(float * d_LocArry_ConsecFit, int TotalFluoNum, int *d_ForwardLinkID, int *d_BackwardLinkID);

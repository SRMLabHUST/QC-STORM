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


#include "bfgsMLE_GS2D_core.h"

#include "bfgsMLE_AS3D_core.h"





void LDLoc_BFGS_MLELocalizationGS2D(float * d_LocArry, unsigned short * d_ImageROI, float *d_WLEPara, int SingleFitNum, FitPosInf_TypeDef* d_FitPosInf, CoreFittingPara *d_FitPara, int ROISize, cudaStream_t cstream)
{
	cudaError_t err;

	int BlockDim = ThreadsPerBlock;
	int BlockNum = ((SingleFitNum + ThreadsPerBlock - 1) / ThreadsPerBlock);
	

	// MLEROILocTop execute
	switch (ROISize)
	{	
	case 5:
		bfgsMLELoc_Gauss2D<5, FitParaNum_2D> << <BlockNum, BlockDim, 0, cstream >> > (d_LocArry, d_ImageROI, d_WLEPara, SingleFitNum, d_FitPosInf, d_FitPara);
		break;
	case 7:
		bfgsMLELoc_Gauss2D<7, FitParaNum_2D> << <BlockNum, BlockDim, 0, cstream >> > (d_LocArry, d_ImageROI, d_WLEPara, SingleFitNum, d_FitPosInf, d_FitPara);
		break;
	case 9:
		bfgsMLELoc_Gauss2D<9, FitParaNum_2D> << <BlockNum, BlockDim, 0, cstream >> > (d_LocArry, d_ImageROI, d_WLEPara, SingleFitNum, d_FitPosInf, d_FitPara);
		break;
	case 11:
		bfgsMLELoc_Gauss2D<11, FitParaNum_2D> << <BlockNum, BlockDim, 0, cstream >> > (d_LocArry, d_ImageROI, d_WLEPara, SingleFitNum, d_FitPosInf, d_FitPara);
		break;
	case 13:
		bfgsMLELoc_Gauss2D<13, FitParaNum_2D> << <BlockNum, BlockDim, 0, cstream >> > (d_LocArry, d_ImageROI, d_WLEPara, SingleFitNum, d_FitPosInf, d_FitPara);
		break;
	case 15:
		bfgsMLELoc_Gauss2D<15, FitParaNum_2D> << <BlockNum, BlockDim, 0, cstream >> > (d_LocArry, d_ImageROI, d_WLEPara, SingleFitNum, d_FitPosInf, d_FitPara);
		break;
	case 17:
		bfgsMLELoc_Gauss2D<17, FitParaNum_2D> << <BlockNum, BlockDim, 0, cstream >> > (d_LocArry, d_ImageROI, d_WLEPara, SingleFitNum, d_FitPosInf, d_FitPara);
		break;

	default:

		break;
	}

	cudaStreamQuery(cstream);

}




void LDLoc_BFGS_MLELocalizationAS3D(float * d_LocArry, unsigned short * d_ImageROI, float *d_WLEPara, int SingleFitNum, FitPosInf_TypeDef* d_FitPosInf, CoreFittingPara *d_FitPara, int ROISize, cudaStream_t cstream)
{
	cudaError_t err;


	int BlockDim = ThreadsPerBlock;
	int BlockNum = ((SingleFitNum + ThreadsPerBlock - 1) / ThreadsPerBlock);


	// MLEROILocTop_AS3D execute
	switch (ROISize)
	{

	case 5:
		bfgsMLELoc_AS3D<5, FitParaNum_AS3D> << <BlockNum, BlockDim, 0, cstream >> >(d_LocArry, d_ImageROI, d_WLEPara, SingleFitNum, d_FitPosInf, d_FitPara);
		break;
	case 7:
		bfgsMLELoc_AS3D<7, FitParaNum_AS3D> << <BlockNum, BlockDim, 0, cstream >> >(d_LocArry, d_ImageROI, d_WLEPara, SingleFitNum, d_FitPosInf, d_FitPara);
		break;
	case 9:
		bfgsMLELoc_AS3D<8, FitParaNum_AS3D> << <BlockNum, BlockDim, 0, cstream >> >(d_LocArry, d_ImageROI, d_WLEPara, SingleFitNum, d_FitPosInf, d_FitPara);
		break;
	case 11:
		bfgsMLELoc_AS3D<11, FitParaNum_AS3D> << <BlockNum, BlockDim, 0, cstream >> >(d_LocArry, d_ImageROI, d_WLEPara, SingleFitNum, d_FitPosInf, d_FitPara);
		break;
	case 13:
		bfgsMLELoc_AS3D<13, FitParaNum_AS3D> << <BlockNum, BlockDim, 0, cstream >> >(d_LocArry, d_ImageROI, d_WLEPara, SingleFitNum, d_FitPosInf, d_FitPara);
		break;
	case 15:
		bfgsMLELoc_AS3D<15, FitParaNum_AS3D> << <BlockNum, BlockDim, 0, cstream >> >(d_LocArry, d_ImageROI, d_WLEPara, SingleFitNum, d_FitPosInf, d_FitPara);
		break;
	case 17:
		bfgsMLELoc_AS3D<17, FitParaNum_AS3D> << <BlockNum, BlockDim, 0, cstream >> >(d_LocArry, d_ImageROI, d_WLEPara, SingleFitNum, d_FitPosInf, d_FitPara);
		break;


	default:

		break;
	}

	cudaStreamQuery(cstream);

}

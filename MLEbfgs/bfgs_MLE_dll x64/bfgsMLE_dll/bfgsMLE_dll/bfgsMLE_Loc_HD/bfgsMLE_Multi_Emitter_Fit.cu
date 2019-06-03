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


#include "bfgsMLE_GS2D_2Emitter_core.h"

#include "bfgsMLE_GS2D_3Emitter_core.h"

#include "bfgsMLE_AS3D_2Emitter_core.h"


#include "bfgsMLE_Multi_Emitter_Fit.h"



void HDLoc_BFGS_MLELocalization_2D_2Emitter(float * d_LocArry, unsigned short * d_ImageROI, float *d_WLEPara, int MultiFitFluoNum_2E, FitPosInf_TypeDef* d_FitPosInf, CoreFittingPara *d_FitPara, int ROISize, cudaStream_t cstream)
{
	cudaError_t err;


	int BlockDim = ThreadsPerBlock;
	int BlockNum = ((MultiFitFluoNum_2E + ThreadsPerBlock - 1) / ThreadsPerBlock);
	
	
	// MLEROILocTop execute
	switch (ROISize)
	{	
	case 5:
		bfgsMLELoc_Gauss2D_2E<5, FitParaNum_2D_2E> << <BlockNum, BlockDim, 0, cstream >> > (d_LocArry, d_ImageROI, d_WLEPara, MultiFitFluoNum_2E, d_FitPosInf, d_FitPara);
		break;
	case 7:
		bfgsMLELoc_Gauss2D_2E<7, FitParaNum_2D_2E> << <BlockNum, BlockDim, 0, cstream >> > (d_LocArry, d_ImageROI, d_WLEPara, MultiFitFluoNum_2E, d_FitPosInf, d_FitPara);
		break;
	case 9:
		bfgsMLELoc_Gauss2D_2E<9, FitParaNum_2D_2E> << <BlockNum, BlockDim, 0, cstream >> > (d_LocArry, d_ImageROI, d_WLEPara, MultiFitFluoNum_2E, d_FitPosInf, d_FitPara);
		break;
	case 11:
		bfgsMLELoc_Gauss2D_2E<11, FitParaNum_2D_2E> << <BlockNum, BlockDim, 0, cstream >> > (d_LocArry, d_ImageROI, d_WLEPara, MultiFitFluoNum_2E, d_FitPosInf, d_FitPara);
		break;
	case 13:
		bfgsMLELoc_Gauss2D_2E<13, FitParaNum_2D_2E> << <BlockNum, BlockDim, 0, cstream >> > (d_LocArry, d_ImageROI, d_WLEPara, MultiFitFluoNum_2E, d_FitPosInf, d_FitPara);
		break;
	
		/*
	// shared memory is not enough
	case 15:

	break;
	case 17:

	break;
		*/
	default:

		break;
	}

	cudaStreamQuery(cstream);

}



void HDLoc_BFGS_MLELocalization_2D_3Emitter(float * d_LocArry, unsigned short * d_ImageROI, float *d_WLEPara, int MultiFitFluoNum_3E, FitPosInf_TypeDef* d_FitPosInf, CoreFittingPara *d_FitPara, int ROISize, cudaStream_t cstream)
{
	cudaError_t err;


	int BlockDim = ThreadsPerBlock;
	int BlockNum = ((MultiFitFluoNum_3E + ThreadsPerBlock - 1) / ThreadsPerBlock);


	// MLEROILocTop execute
	switch (ROISize)
	{
	case 5:
		bfgsMLELoc_Gauss2D_3E<5, FitParaNum_2D_3E> << <BlockNum, BlockDim, 0, cstream >> > (d_LocArry, d_ImageROI, d_WLEPara, MultiFitFluoNum_3E, d_FitPosInf, d_FitPara);
		break;
	case 7:
		bfgsMLELoc_Gauss2D_3E<7, FitParaNum_2D_3E> << <BlockNum, BlockDim, 0, cstream >> > (d_LocArry, d_ImageROI, d_WLEPara, MultiFitFluoNum_3E, d_FitPosInf, d_FitPara);
		break;
	case 9:
		bfgsMLELoc_Gauss2D_3E<9, FitParaNum_2D_3E> << <BlockNum, BlockDim, 0, cstream >> > (d_LocArry, d_ImageROI, d_WLEPara, MultiFitFluoNum_3E, d_FitPosInf, d_FitPara);
		break;
	case 11:
		bfgsMLELoc_Gauss2D_3E<11, FitParaNum_2D_3E> << <BlockNum, BlockDim, 0, cstream >> > (d_LocArry, d_ImageROI, d_WLEPara, MultiFitFluoNum_3E, d_FitPosInf, d_FitPara);
		break;
	case 13:
		bfgsMLELoc_Gauss2D_3E<13, FitParaNum_2D_3E> << <BlockNum, BlockDim, 0, cstream >> > (d_LocArry, d_ImageROI, d_WLEPara, MultiFitFluoNum_3E, d_FitPosInf, d_FitPara);
		break;

		/*
		// shared memory is not enough
		case 15:

		break;
		case 17:

		break;
		*/
	default:

		break;
	}

	cudaStreamQuery(cstream);

}



void HDLoc_BFGS_MLELocalization_AS3D_2Emitter(float * d_LocArry, unsigned short * d_ImageROI, float *d_WLEPara, int MultiFitFluoNum_2E, FitPosInf_TypeDef* d_FitPosInf, CoreFittingPara *d_FitPara, int ROISize, cudaStream_t cstream)
{
	cudaError_t err;


	int BlockDim = ThreadsPerBlock;
	int BlockNum = ((MultiFitFluoNum_2E + ThreadsPerBlock - 1) / ThreadsPerBlock);


	// MLEROILocTop_AS3D execute
	switch (ROISize)
	{

	case 5:
		bfgsMLELoc_AS3D_2E<5, FitParaNum_AS3D_2E> << <BlockNum, BlockDim, 0, cstream >> >(d_LocArry, d_ImageROI, d_WLEPara, MultiFitFluoNum_2E, d_FitPosInf, d_FitPara);
		break;
	case 7:
		bfgsMLELoc_AS3D_2E<7, FitParaNum_AS3D_2E> << <BlockNum, BlockDim, 0, cstream >> >(d_LocArry, d_ImageROI, d_WLEPara, MultiFitFluoNum_2E, d_FitPosInf, d_FitPara);
		break;
	case 9:
		bfgsMLELoc_AS3D_2E<8, FitParaNum_AS3D_2E> << <BlockNum, BlockDim, 0, cstream >> >(d_LocArry, d_ImageROI, d_WLEPara, MultiFitFluoNum_2E, d_FitPosInf, d_FitPara);
		break;
	case 11:
		bfgsMLELoc_AS3D_2E<11, FitParaNum_AS3D_2E> << <BlockNum, BlockDim, 0, cstream >> >(d_LocArry, d_ImageROI, d_WLEPara, MultiFitFluoNum_2E, d_FitPosInf, d_FitPara);
		break;
	case 13:
		bfgsMLELoc_AS3D_2E<13, FitParaNum_AS3D_2E> << <BlockNum, BlockDim, 0, cstream >> >(d_LocArry, d_ImageROI, d_WLEPara, MultiFitFluoNum_2E, d_FitPosInf, d_FitPara);
		break;

		/*
		// shared memory is not enough
		case 15:

		break;
		case 17:

		break;

		*/

	default:

		break;
	}

	cudaStreamQuery(cstream);

}

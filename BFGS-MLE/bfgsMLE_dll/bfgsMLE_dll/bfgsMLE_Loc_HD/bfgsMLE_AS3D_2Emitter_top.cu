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



#include "bfgsMLE_AS3D_2Emitter_core.h"


#include "bfgsMLE_Multi_Emitter_Fit.h"


void HDLoc_BFGS_MLELocalization_AS3D_2Emitter(float * d_LocArry, unsigned short * d_ImageROI, float *d_WLEPara, int MultiFitFluoNum, int * d_MultiFitFluoPos, LocalizationPara& LocPara, cudaStream_t cstream)
{
	cudaError_t err;
	// MaxFluoNum must be the integer multiples of 32
	int ROISize = LocPara.ROISize;


	int BlockDim = ThreadsPerBlock;
	int BlockNum = ((MultiFitFluoNum + ThreadsPerBlock - 1) / ThreadsPerBlock);


	// MLEROILocTop_AS3D execute
	switch (ROISize)
	{

	case 5:
		bfgsMLELoc_AS3D_2E<5, FitParaNum_AS3D_2E> << <BlockNum, BlockDim, 0, cstream >> >(d_LocArry, d_ImageROI, d_WLEPara, MultiFitFluoNum, d_MultiFitFluoPos, LocPara.Offset, LocPara.KAdc, LocPara.QE, LocPara.ZDepthCorrFactor, LocPara.p4, LocPara.p3, LocPara.p2, LocPara.p1, LocPara.p0);
		break;
	case 7:
		bfgsMLELoc_AS3D_2E<7, FitParaNum_AS3D_2E> << <BlockNum, BlockDim, 0, cstream >> >(d_LocArry, d_ImageROI, d_WLEPara, MultiFitFluoNum, d_MultiFitFluoPos, LocPara.Offset, LocPara.KAdc, LocPara.QE, LocPara.ZDepthCorrFactor, LocPara.p4, LocPara.p3, LocPara.p2, LocPara.p1, LocPara.p0);
		break;
	case 9:
		bfgsMLELoc_AS3D_2E<8, FitParaNum_AS3D_2E> << <BlockNum, BlockDim, 0, cstream >> >(d_LocArry, d_ImageROI, d_WLEPara, MultiFitFluoNum, d_MultiFitFluoPos, LocPara.Offset, LocPara.KAdc, LocPara.QE, LocPara.ZDepthCorrFactor, LocPara.p4, LocPara.p3, LocPara.p2, LocPara.p1, LocPara.p0);
		break;
	case 11:
		bfgsMLELoc_AS3D_2E<11, FitParaNum_AS3D_2E> << <BlockNum, BlockDim, 0, cstream >> >(d_LocArry, d_ImageROI, d_WLEPara, MultiFitFluoNum, d_MultiFitFluoPos, LocPara.Offset, LocPara.KAdc, LocPara.QE, LocPara.ZDepthCorrFactor, LocPara.p4, LocPara.p3, LocPara.p2, LocPara.p1, LocPara.p0);
		break;
	case 13:
		bfgsMLELoc_AS3D_2E<13, FitParaNum_AS3D_2E> << <BlockNum, BlockDim, 0, cstream >> >(d_LocArry, d_ImageROI, d_WLEPara, MultiFitFluoNum, d_MultiFitFluoPos, LocPara.Offset, LocPara.KAdc, LocPara.QE, LocPara.ZDepthCorrFactor, LocPara.p4, LocPara.p3, LocPara.p2, LocPara.p1, LocPara.p0);
		break;
		
		/*
		// shared memory is not enough
		case 15:
		bfgsMLELoc_AS3D_2E<15, FitParaNum_AS3D_2E> << <BlockNum, BlockDim, 0, cstream >> >(d_LocArry, d_ImageROI, d_WLEPara, MultiFitFluoNum, d_MultiFitFluoPos, LocPara.Offset, LocPara.KAdc, LocPara.QE, LocPara.ZDepthCorrFactor, LocPara.p4, LocPara.p3, LocPara.p2, LocPara.p1, LocPara.p0);
		break;
	case 17:
		bfgsMLELoc_AS3D_2E<17, FitParaNum_AS3D_2E> << <BlockNum, BlockDim, 0, cstream >> >(d_LocArry, d_ImageROI, d_WLEPara, MultiFitFluoNum, d_MultiFitFluoPos, LocPara.Offset, LocPara.KAdc, LocPara.QE, LocPara.ZDepthCorrFactor, LocPara.p4, LocPara.p3, LocPara.p2, LocPara.p1, LocPara.p0);
		break;

		*/

	default:

		break;
	}

	cudaStreamQuery(cstream);

}

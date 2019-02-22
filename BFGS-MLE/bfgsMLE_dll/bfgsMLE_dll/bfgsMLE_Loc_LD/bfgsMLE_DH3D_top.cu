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

#include "bfgsMLE_DH3D_top.h"



// top application use

void LDLoc_BFGS_MLELocalizationDH3D(unsigned short * d_SubRegion, float * d_LocArry, LocalizationPara& LocPara, int FluoNum, cudaStream_t cstream)
{
	cudaError_t err;
	// MaxFluoNum must be the integer multiples of 32
	int ROISize = LocPara.ROISize;

	int BlockDim = ThreadsPerBlock;
	int BlockNum = ((FluoNum + ThreadsPerBlock - 1) / ThreadsPerBlock);

	
	// MLEROILocTop execute
	switch (ROISize)
	{	
	case 5:
		bfgsMLELoc_DH3D<5, FitParaNum_DH3D> << <BlockNum, BlockDim, 0, cstream >> >(d_SubRegion, d_LocArry, LocPara.Offset, LocPara.KAdc, LocPara.QE, FluoNum);
		break;
	case 7:
		bfgsMLELoc_DH3D<7, FitParaNum_DH3D> << <BlockNum, BlockDim, 0, cstream >> >(d_SubRegion, d_LocArry, LocPara.Offset, LocPara.KAdc, LocPara.QE, FluoNum);
		break;
	case 9:
		bfgsMLELoc_DH3D<9, FitParaNum_DH3D> << <BlockNum, BlockDim, 0, cstream >> >(d_SubRegion, d_LocArry, LocPara.Offset, LocPara.KAdc, LocPara.QE, FluoNum);
		break;
	case 11:
		bfgsMLELoc_DH3D<11, FitParaNum_DH3D> << <BlockNum, BlockDim, 0, cstream >> >(d_SubRegion, d_LocArry, LocPara.Offset, LocPara.KAdc, LocPara.QE, FluoNum);
		break;
	case 13:
		bfgsMLELoc_DH3D<13, FitParaNum_DH3D> << <BlockNum, BlockDim, 0, cstream >> >(d_SubRegion, d_LocArry, LocPara.Offset, LocPara.KAdc, LocPara.QE, FluoNum);
		break;
	case 15:
		bfgsMLELoc_DH3D<15, FitParaNum_DH3D> << <BlockNum, BlockDim, 0, cstream >> >(d_SubRegion, d_LocArry, LocPara.Offset, LocPara.KAdc, LocPara.QE, FluoNum);
		break;
	case 17:
		bfgsMLELoc_DH3D<17, FitParaNum_DH3D> << <BlockNum, BlockDim, 0, cstream >> >(d_SubRegion, d_LocArry, LocPara.Offset, LocPara.KAdc, LocPara.QE, FluoNum);
		break;
	default:

		break;
	}

	cudaStreamQuery(cstream);

}

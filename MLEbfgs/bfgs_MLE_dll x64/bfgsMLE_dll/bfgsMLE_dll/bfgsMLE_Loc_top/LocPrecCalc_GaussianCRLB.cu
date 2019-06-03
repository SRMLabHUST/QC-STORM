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


#include "LocPrecCalc_GaussianCRLB.h"

#include "bfgs_top.h"



void LDLocData_TypeDef::LocPrecCalc_GaussianCRLB(float* d_LocArry, LocalizationPara & LocPara, int FluoNum, cudaStream_t cstream)
{
	// only suitable for 2d localization with symmetric Gaussian PSF

	int BlockDim = ThreadsPerBlock;
	int BlockNum = ((FluoNum + ThreadsPerBlock - 1) / ThreadsPerBlock);

	if (LocPara.LocType != LocType_AS3D)
	{
		GaussianCRLB_Calc_top<GaussianCRLB_ROISize2D> << <BlockNum, BlockDim, 0, cstream >> >(d_LocArry, LocPara.ReadNoise_e, LocPara.QE, LocPara.PixelSize);
	}
	else
	{
		GaussianCRLB_Calc_top<GaussianCRLB_ROISize3D> << <BlockNum, BlockDim, 0, cstream >> >(d_LocArry, LocPara.ReadNoise_e, LocPara.QE, LocPara.PixelSize);
	}

	cudaStreamSynchronize(cstream);
}



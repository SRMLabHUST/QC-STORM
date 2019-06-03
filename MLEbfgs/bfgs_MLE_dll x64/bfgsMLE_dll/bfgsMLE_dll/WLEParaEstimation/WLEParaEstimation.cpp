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

#include "stdafx.h"
#include "WLEParaEstimation.h"


void WLEParameterEstimation_TypeDef::WLEParameterEstimate(unsigned short * d_ImageROI, int LocType, int MultiEmitterFitEn, int ROISize, int FluoNum, cudaStream_t cstream)
{
	const int ROIWholeSize = ROISize*(ROISize + 1);

	// copy from CPU
//	cudaMemcpyAsync(d_ImageROI, h_ImageROI, FluoNum * ROIWholeSize * sizeof(unsigned short), cudaMemcpyHostToDevice, cstream);


	CalculatePSFWidth(d_ImageROI, d_WLEPara, FluoNum, ROISize, cstream);

	CalculateNearestNeighborDistance(d_ImageROI, ROISize, d_WLEPara, FluoNum, cstream);


	MoleculeTypeClasify(d_ImageROI, LocType, MultiEmitterFitEn, ROISize, d_WLEPara, FluoNum, cstream);

//	cudaMemcpyAsync(h_WLEPara, d_WLEPara, FluoNum * WLE_ParaNumber * sizeof(float), cudaMemcpyDeviceToHost, cstream);


	cudaStreamSynchronize(cstream);

}


void WLEParameterEstimation_TypeDef::Init(LocalizationPara & LocPara)
{
	cudaError_t err;
	const int ROIWholeSize = LocPara.ROISize*(LocPara.ROISize + 1);

//	err = cudaMalloc((void **)&d_ImageROI, MaxPointNum * ROIWholeSize * sizeof(unsigned short));

//	HandleErr(err, "cudaMalloc d_ImageROI");

	err = cudaMallocHost((void **)&h_WLEPara, MaxPointNum * WLE_ParaNumber * sizeof(float));
	err = cudaMalloc((void **)&d_WLEPara, MaxPointNum * WLE_ParaNumber * sizeof(float));

//	HandleErr(err, "cudaMalloc d_WLEPara");


}

void WLEParameterEstimation_TypeDef::Deinit()
{
	cudaError_t err;
//	err = cudaFree(d_ImageROI);

	err = cudaFreeHost(h_WLEPara);
	err = cudaFree(d_WLEPara);


}


#include "stdafx.h"
#include "WLEParaEstimation.h"


void WLEParameterEstimation_TypeDef::WLEParameterEstimate(unsigned short * h_ROIMem, int ROISize, int FluoNum, cudaStream_t cstream)
{
	const int ROIWholeSize = ROISize*(ROISize + 1);

	// copy from CPU
	cudaMemcpyAsync(d_ROIMem, h_ROIMem, FluoNum * ROIWholeSize * sizeof(unsigned short), cudaMemcpyHostToDevice, cstream);


	CalculatePSFWidth(d_ROIMem, d_WLEPara, FluoNum, ROISize, cstream);

	CalculateNearestNeighborDistance(d_ROIMem, ROISize, d_WLEPara, FluoNum, cstream);

	MoleculeTypeClasify(ROISize, d_WLEPara, FluoNum, cstream);


//	cudaMemcpyAsync(h_WLEPara, d_WLEPara, FluoNum * WLE_ParaNumber * sizeof(float), cudaMemcpyDeviceToHost, cstream);


	cudaStreamSynchronize(cstream);

}


void WLEParameterEstimation_TypeDef::Init(LocalizationPara & LocPara)
{
	cudaError_t err;
	const int ROIWholeSize = LocPara.ROISize*(LocPara.ROISize + 1);

	err = cudaMalloc((void **)&d_ROIMem, MaxPointNum * ROIWholeSize * sizeof(unsigned short));

	HandleErr(err, "cudaMalloc d_ROIMem");

	err = cudaMallocHost((void **)&h_WLEPara, MaxPointNum * WLE_ParaNumber * sizeof(float));
	err = cudaMalloc((void **)&d_WLEPara, MaxPointNum * WLE_ParaNumber * sizeof(float));

	HandleErr(err, "cudaMalloc d_WLEPara");


}

void WLEParameterEstimation_TypeDef::Deinit()
{
	cudaError_t err;
	err = cudaFree(d_ROIMem);

	err = cudaFreeHost(h_WLEPara);
	err = cudaFree(d_WLEPara);


}


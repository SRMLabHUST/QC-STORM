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

#include "DH3D_MoleculePair.h"


void DH3DPairData_TypeDef::PairMolecules(float *d_LocArry, LocalizationPara & LocPara, int FluoNum, cudaStream_t cstream)
{

	int PosArryEnable = 1;

	cudaMemsetAsync(d_ValidFluoNum, 0, sizeof(int), cstream);

//	cudaMemcpyAsync(d_LocArry, h_LocArry, FluoNum*OutParaNumGS2D*sizeof(float), cudaMemcpyHostToDevice, cstream);


	int BlockDim = ThreadsPerBlock;
	int BlockNum = ((FluoNum + ThreadsPerBlock - 1) / ThreadsPerBlock);


	gpuPairMolecules << <BlockNum, BlockDim, 0, cstream >> >(d_LocArry, d_PairedLocArry, d_PairedPosArry, d_ValidFluoNum, LocPara.MeanDistance, LocPara.DistanceTh, FluoNum, PosArryEnable, LocPara.p4, LocPara.p3, LocPara.p2, LocPara.p1, LocPara.p0, LocPara.RotateType);

	cudaMemcpyAsync(h_ValidFluoNum, d_ValidFluoNum, sizeof(int), cudaMemcpyDeviceToHost, cstream);
	cudaStreamSynchronize(cstream);

	oValidFluoNum = *(h_ValidFluoNum);

	// get paired localization results
//	cudaMemcpyAsync(h_PairedLocArry, d_PairedLocArry, oValidFluoNum*OutParaNumGS2D*sizeof(float), cudaMemcpyDeviceToHost, cstream);
//	cudaStreamSynchronize(cstream);


	if (PosArryEnable)
	{
		// get paired molecules group
		BlockNum = ((oValidFluoNum + ThreadsPerBlock - 1) / ThreadsPerBlock);
		gpuCalcDistanceDistribution << <BlockNum, BlockDim, 0, cstream >> >(d_PairedPosArry, d_DistanceDistrib, oValidFluoNum);
		
		cudaMemcpyAsync(h_PairedPosArry, d_PairedPosArry, oValidFluoNum*OutParaNumGS2D*sizeof(float), cudaMemcpyDeviceToHost, cstream);
	
		cudaMemcpyAsync(h_DistanceDistrib, d_DistanceDistrib, DistanceHistLen*sizeof(int), cudaMemcpyDeviceToHost, cstream);
		cudaStreamSynchronize(cstream);
		
		
		MeanDistance = DistanceGap*FluoStatisticData_TypeDef::GetHistogramMaxDataPos(h_DistanceDistrib, DistanceHistLen);
		DistanceHistWidth = 2.355f* DistanceGap*FluoStatisticData_TypeDef::GetHistogramWidth(h_DistanceDistrib, FluoStatisticData_TypeDef::GetHistogramMaxDataPos(h_DistanceDistrib, DistanceHistLen), DistanceHistLen);

		
		printf("Double-Helix 3d, Mean and Width th for pairs distance th:%f %f\n", MeanDistance, DistanceHistWidth);
	}
	
	printf("Double-Helix oval fluonum: %d - %d\n", FluoNum, oValidFluoNum);
}



void DH3DPairData_TypeDef::Init()
{
	cudaError_t err;


	// host and gpu

	cudaMallocHost((void **)&h_PairedLocArry, MaxPointNum*OutParaNumAS3D*sizeof(float));
	cudaMalloc((void **)&d_PairedLocArry, MaxPointNum*OutParaNumAS3D*sizeof(float));

	cudaMallocHost((void **)&h_PairedPosArry, MaxPointNum*OutParaNumAS3D*sizeof(float));
	cudaMalloc((void **)&d_PairedPosArry, MaxPointNum*OutParaNumAS3D*sizeof(float));

	cudaMallocHost((void **)&h_ValidFluoNum, sizeof(int));
	cudaMalloc((void **)&d_ValidFluoNum, sizeof(int));

	cudaMallocHost((void **)&h_DistanceDistrib, DistanceHistLen*sizeof(int));
	cudaMalloc((void **)&d_DistanceDistrib, DistanceHistLen*sizeof(int));


}

void DH3DPairData_TypeDef::Deinit()
{

	cudaFreeHost(h_PairedLocArry);
	cudaFree(d_PairedLocArry);

	cudaFreeHost(h_PairedPosArry);
	cudaFree(d_PairedPosArry);

	cudaFreeHost(h_ValidFluoNum);
	cudaFree(d_ValidFluoNum);

	cudaFreeHost(h_DistanceDistrib);
	cudaFree(d_DistanceDistrib);

}


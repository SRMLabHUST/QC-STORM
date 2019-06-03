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

#include "ConsecutiveFit.h"
#include "bfgs_top.h"



// d_iLocArry come from localization data 
void ConsecutiveFit_TypeDef::ConsecutiveFit_WeightedAvg(float * h_iLocArry, int FluoNum_CurGrop, int IsEndProc, LocalizationPara & LocPara, cudaStream_t cstream)
{
	int TotalFluoNum = FluoNum_LastGroup + FluoNum_CurGrop;

	// copy data, find consecutive molecules should process the acjecent molecules between two independent image groups
	cudaMemcpyAsync(d_LocArry_ConsecFit, d_LocArry_LastGroup, FluoNum_LastGroup*OutParaNumGS2D*sizeof(float), cudaMemcpyDeviceToDevice, cstream);

	cudaMemcpyAsync(&d_LocArry_ConsecFit[FluoNum_LastGroup*OutParaNumGS2D], h_iLocArry, FluoNum_CurGrop*OutParaNumGS2D * sizeof(float), cudaMemcpyHostToDevice, cstream);


	// find consecutive molecules
	float DistanceTh_Pixel = LocPara.ConsecFit_DistanceTh_nm / LocPara.PixelSize;

	int BlockDim = ThreadsPerBlock;
	int BlockNum = ((TotalFluoNum + ThreadsPerBlock - 1) / ThreadsPerBlock);
	
	cudaMemsetAsync(d_ForwardLinkID, 0, TotalFluoNum * sizeof(int), cstream);
	cudaMemsetAsync(d_BackwardLinkID, 0, TotalFluoNum * sizeof(int), cstream);

	// ontime calculation
	cudaMemsetAsync(d_OntimeDistrib, 0, MaxOnTimeConsecutiveNum * sizeof(int), cstream);
	cudaMemsetAsync(d_ValidFluoNum, 0, sizeof(int), cstream);


	gpuFindConsecutiveMolecules << <BlockNum, BlockDim, 0, cstream >> > (d_LocArry_ConsecFit, TotalFluoNum, d_ForwardLinkID, d_BackwardLinkID, DistanceTh_Pixel);
	
	gpuConsecutiveFit << <BlockNum, BlockDim, 0, cstream >> > (d_LocArry_ConsecFit, TotalFluoNum, FluoNum_LastGroup, d_ForwardLinkID, d_BackwardLinkID, d_OntimeDistrib, d_ValidFluoNum, LocPara.QE);


	gpuRemoveConsecutiveFluo << <BlockNum, BlockDim, 0, cstream >> > (d_LocArry_ConsecFit, TotalFluoNum, d_ForwardLinkID, d_BackwardLinkID);

	cudaMemcpyAsync(h_ValidFluoNum, d_ValidFluoNum, sizeof(int), cudaMemcpyDeviceToHost, cstream);
	cudaMemcpyAsync(h_OntimeDistrib, d_OntimeDistrib, MaxOnTimeConsecutiveNum * sizeof(int), cudaMemcpyDeviceToHost, cstream);


	cudaStreamSynchronize(cstream);

	// ontime distribution ratio
	if (*h_ValidFluoNum <= 0)*h_ValidFluoNum = 1;
	for (int cnt = 0; cnt < MaxOnTimeConsecutiveNum; cnt++)
	{
		h_OntimeRatio[cnt] = h_OntimeDistrib[cnt] * 1.0f / (*(h_ValidFluoNum));
	}


	// localization precision update for consecutive fitted molecules
	LDLocData_TypeDef::LocPrecCalc_GaussianCRLB(d_LocArry_ConsecFit, LocPara, FluoNum_LastGroup, cstream);


	if (IsEndProc == 0)
	{
		// output last group
		OutFluoNum = FluoNum_LastGroup;

		cudaMemcpyAsync(d_LocArry_LastGroup, &d_LocArry_ConsecFit[FluoNum_LastGroup*OutParaNumGS2D], FluoNum_CurGrop*OutParaNumGS2D * sizeof(float), cudaMemcpyDeviceToDevice, cstream);
		
		FluoNum_LastGroup = FluoNum_CurGrop;
	}
	else
	{
		// output all
		OutFluoNum = TotalFluoNum;
		FluoNum_LastGroup = 0;
	}

	cudaMemcpyAsync(h_OutLocArry, d_LocArry_ConsecFit, OutFluoNum*OutParaNumGS2D * sizeof(float), cudaMemcpyDeviceToHost, cstream);
	cudaStreamSynchronize(cstream);

}



void ConsecutiveFit_TypeDef::Init()
{
	OutFluoNum = 0;
	FluoNum_LastGroup = 0;

	cudaMallocHost((void **)&h_OutLocArry, 2 * MaxPointNum * OutParaNumGS2D * sizeof(float));

	cudaMalloc((void **)&d_LocArry_LastGroup, MaxPointNum * OutParaNumGS2D * sizeof(float));
	
	cudaMalloc((void **)&d_LocArry_ConsecFit, 2 * MaxPointNum * OutParaNumGS2D * sizeof(float));


	// consecutive molecule find

	cudaMalloc((void **)&d_ForwardLinkID, 2 * MaxPointNum * sizeof(int));
	cudaMalloc((void **)&d_BackwardLinkID, 2 * MaxPointNum * sizeof(int));


	// on time calculate

	cudaMallocHost((void **)&h_OntimeDistrib, MaxOnTimeConsecutiveNum * sizeof(int));
	cudaMalloc((void **)&d_OntimeDistrib, MaxOnTimeConsecutiveNum * sizeof(int));

	cudaMallocHost((void **)&h_ValidFluoNum, sizeof(int));
	cudaMalloc((void **)&d_ValidFluoNum, sizeof(int));

	h_OntimeRatio = new float[MaxOnTimeConsecutiveNum];

}

void ConsecutiveFit_TypeDef::Deinit()
{
	cudaFreeHost(h_OutLocArry);

	cudaFree(d_LocArry_LastGroup);

	cudaFree(d_LocArry_ConsecFit);

	// consecutive molecule find
	cudaFree(d_ForwardLinkID);
	cudaFree(d_BackwardLinkID);


	// on time calculate

	cudaFreeHost(h_OntimeDistrib);
	cudaFree(d_OntimeDistrib);

	cudaFreeHost(h_ValidFluoNum);
	cudaFree(d_ValidFluoNum);

	delete[]h_OntimeRatio;
}




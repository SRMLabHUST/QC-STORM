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
void ConsecutiveFit_TypeDef::FitConsecutiveFluo(float * d_iLocArry, LocalizationPara & LocPara, int FluoNum, cudaStream_t cstream, int IsLastProc)
{
	int TotalFluoNum;


	int LastId = (CurInputID + 1) % 2;

	BufFluoNum[CurInputID] = FluoNum;
	float *d_CurLocArry = d_BufLocArry[CurInputID];
	int CurFluoNum = BufFluoNum[CurInputID];

	float *d_LastLocArry = d_BufLocArry[LastId];
	int LastFluoNum = BufFluoNum[LastId];
	int ConsecFitFluoNum = 0;


	//
	CurInputID = (CurInputID + 1) % 2;
	TotalFluoNum = LastFluoNum + CurFluoNum;

	//
	cudaMemcpyAsync(d_CurLocArry, d_iLocArry, FluoNum*OutParaNumGS2D*sizeof(float), cudaMemcpyDeviceToDevice, cstream);


	cudaMemsetAsync(d_ForwardLinkID, 0, TotalFluoNum * sizeof(int), cstream);
	cudaMemsetAsync(d_BackwardLinkID, 0, TotalFluoNum * sizeof(int), cstream);


	int BlockDim = ThreadsPerBlock;
	int BlockNum = ((TotalFluoNum + ThreadsPerBlock - 1) / ThreadsPerBlock);
	
	gpuFindConsecutiveMolecules << <BlockNum, BlockDim, 0, cstream >> >(d_LastLocArry, LastFluoNum, d_CurLocArry, CurFluoNum, d_ForwardLinkID, d_BackwardLinkID, LocPara.ConsecFit_DistanceTh_nm, LocPara.PixelSize, LocPara.QE);
	

	// only the bigining part of 2th part are fited with the first part, except the last time 
	if (IsLastProc == 0)
	{
		ConsecFitFluoNum = LastFluoNum;
	}
	else
	{
		ConsecFitFluoNum = TotalFluoNum;
	}


	BlockNum = ((ConsecFitFluoNum + ThreadsPerBlock - 1) / ThreadsPerBlock);

	gpuConsecutiveFit << <BlockNum, BlockDim, 0, cstream >> >(d_LastLocArry, LastFluoNum, d_CurLocArry, CurFluoNum, d_ForwardLinkID, d_BackwardLinkID, ConsecFitFluoNum, LocPara.QE);

	gpuRemoveConsecutiveFluo << <BlockNum, BlockDim, 0, cstream >> >(d_LastLocArry, LastFluoNum, d_CurLocArry, CurFluoNum, d_ForwardLinkID, d_BackwardLinkID, ConsecFitFluoNum);



	// localization precision calculated by CRLB
	LDLocData_TypeDef::LocPrecCalc_GaussianCRLB(d_LastLocArry, LocPara, LastFluoNum, cstream);

	cudaStreamSynchronize(cstream);


}

void ConsecutiveFit_TypeDef::GetAvailableData(cudaStream_t cstream)
{
	int LastID = CurInputID;
	float *d_CurLocArry = d_BufLocArry[LastID];
	OutFluoNum = BufFluoNum[LastID];

	cudaMemcpyAsync(h_OutLocArry, d_CurLocArry, OutFluoNum*OutParaNumGS2D*sizeof(float), cudaMemcpyDeviceToHost, cstream);
	cudaStreamSynchronize(cstream);

}


void ConsecutiveFit_TypeDef::GetResidualData(cudaStream_t cstream)
{
	int LastID = CurInputID;
	int LastFluoNum = BufFluoNum[LastID];

	int CurId = (CurInputID + 1) % 2;
	int CurFluoNum = BufFluoNum[CurId];

	OutFluoNum = LastFluoNum + CurFluoNum;


	int AddrOffset = 0;
	// 1 part
	float *d_CurLocArry;

	d_CurLocArry = d_BufLocArry[LastID];
	AddrOffset = 0;
	cudaMemcpyAsync(&h_OutLocArry[AddrOffset], d_CurLocArry, LastFluoNum*OutParaNumGS2D*sizeof(float), cudaMemcpyDeviceToHost, cstream);

	// 2 part
	d_CurLocArry = d_BufLocArry[CurId];
	AddrOffset = LastFluoNum*OutParaNumGS2D;
	cudaMemcpyAsync(&h_OutLocArry[AddrOffset], d_CurLocArry, CurFluoNum*OutParaNumGS2D*sizeof(float), cudaMemcpyDeviceToHost, cstream);

	cudaStreamSynchronize(cstream);


}

void ConsecutiveFit_TypeDef::Init()
{

	cudaMallocHost((void **)&h_OutLocArry, PointNumTh * 3 * OutParaNumGS2D*sizeof(float));

	cudaMalloc((void **)&d_BufLocArry[0], PointNumTh * 3 * OutParaNumGS2D*sizeof(float));
	cudaMalloc((void **)&d_BufLocArry[1], PointNumTh * 3 * OutParaNumGS2D*sizeof(float));

	cudaMalloc((void **)&d_ForwardLinkID, 2 * PointNumTh * 3 * sizeof(int));
	cudaMalloc((void **)&d_BackwardLinkID, 2 * PointNumTh * 3 * sizeof(int));


	ResetData();

}

void ConsecutiveFit_TypeDef::Deinit()
{
	cudaFreeHost(h_OutLocArry);

	cudaFree(d_BufLocArry[0]);
	cudaFree(d_BufLocArry[1]);

	cudaFree(d_ForwardLinkID);
	cudaFree(d_BackwardLinkID);

}

void ConsecutiveFit_TypeDef::ResetData()
{
	CurInputID = 0;
	BufFluoNum[0] = 0;
	BufFluoNum[1] = 0;
	OutFluoNum = 0;

	cudaStream_t cstream;
	cudaStreamCreate(&cstream);

	cudaMemsetAsync(d_BufLocArry[0], 0, MaxPointNum*OutParaNumGS2D*sizeof(float), cstream);
	cudaMemsetAsync(d_BufLocArry[1], 0, MaxPointNum*OutParaNumGS2D*sizeof(float), cstream);

	cudaStreamDestroy(cstream);

}


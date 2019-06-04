#include "ConsecutiveFluoMerge.h"




void ConsecutiveFluoMerger_TypeDef::MergeConsecutiveFluo(float * h_LocArry, int FluoNum, LocalizationPara & LocPara, int FilterMode, float Distance_th_pixel, cudaStream_t cstream)
{

	int BlockDim = ThreadsPerBlock;
	int BlockNum = (FluoNum + ThreadsPerBlock - 1) / ThreadsPerBlock;


	cudaMemcpyAsync(d_LocArry, h_LocArry, FluoNum*OutParaNumGS2D*sizeof(float), cudaMemcpyHostToDevice, cstream);


	cudaMemsetAsync(d_ForwardLinkID, 0, FluoNum * sizeof(int), cstream);
	cudaMemsetAsync(d_BackwardLinkID, 0, FluoNum * sizeof(int), cstream);

	
	gpuFindConsecutiveFilterPair << <BlockNum, BlockDim, 0, cstream >> >(d_LocArry, d_ForwardLinkID, d_BackwardLinkID, Distance_th_pixel, FluoNum);

	gpuConsecutiveFit << <BlockNum, BlockDim, 0, cstream >> > (d_LocArry, d_ForwardLinkID, d_BackwardLinkID, LocPara.QE, FluoNum);

	gpuRemoveConsecutiveFluo_KeepFirst << <BlockNum, BlockDim, 0, cstream >> >(d_LocArry, d_ForwardLinkID, d_BackwardLinkID, FluoNum);

	// localization precision calculated by CRLB
	LDLocData_TypeDef::LocPrecCalc_GaussianCRLB(d_LocArry, LocPara, FluoNum, cstream);

	cudaMemcpyAsync(h_LocArry, d_LocArry, FluoNum*OutParaNumGS2D*sizeof(float), cudaMemcpyDeviceToHost, cstream);


	cudaStreamSynchronize(cstream);

}



void ConsecutiveFluoMerger_TypeDef::Init(unsigned int TotalFluoNum)
{
	cudaError_t err;

	printf("TotalFluoNum:%d\n", TotalFluoNum);

	err = cudaMallocHost((void **)&h_LocArry, TotalFluoNum * OutParaNumGS2D * sizeof(float));
	HandleErr(err,"cudaMallocHost h_LocArry");

	err = cudaMalloc((void **)&d_LocArry, TotalFluoNum * OutParaNumGS2D*sizeof(float));
	HandleErr(err, "cudaMalloc h_LocArry");


	cudaMalloc((void **)&d_ForwardLinkID, TotalFluoNum * sizeof(int));
	cudaMalloc((void **)&d_BackwardLinkID, TotalFluoNum * sizeof(int));

}


void ConsecutiveFluoMerger_TypeDef::DeInit()
{
	cudaFreeHost(h_LocArry);
	cudaFree(d_LocArry);


	cudaFree(d_ForwardLinkID);
	cudaFree(d_BackwardLinkID);

}



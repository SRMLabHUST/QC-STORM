#include "DH3D_MoleculePair.h"


void DH3D_MoleculePair_TypeDef::MoleculePair(float *h_iLocArry, int FluoNum, LocalizationPara & LocPara, cudaStream_t cstream)
{
	cudaMemsetAsync(d_ValidoFluoNum, 0, sizeof(int), cstream);

	cudaMemcpyAsync(d_LocArry, h_iLocArry, FluoNum*OutParaNumGS2D * sizeof(float), cudaMemcpyHostToDevice, cstream);


	// pair molecules
	int BlockDim = ThreadsPerBlock;
	int BlockNum = ((FluoNum + ThreadsPerBlock - 1) / ThreadsPerBlock);

	gpu_MoleculePair << <BlockNum, BlockDim, 0, cstream >> >(d_LocArry, FluoNum, d_PairID, d_ValidoFluoNum, LocPara.DH_MeanDistance, LocPara.DH_DistanceTh);

	cudaMemcpyAsync(h_ValidoFluoNum, d_ValidoFluoNum, sizeof(int), cudaMemcpyDeviceToHost, cstream);
	cudaStreamSynchronize(cstream);

	oValidFluoNum = *h_ValidoFluoNum;

	// convert paired molecules into single localization


	BlockDim = ThreadsPerBlock;
	BlockNum = ((oValidFluoNum + ThreadsPerBlock - 1) / ThreadsPerBlock);

	gpu_MoleculeMerge << <BlockNum, BlockDim, 0, cstream >> > (d_LocArry, d_oLocArry, d_PairID, oValidFluoNum, LocPara.DH_RotateType, LocPara.p4_XGY, LocPara.p3_XGY, LocPara.p2_XGY, LocPara.p1_XGY, LocPara.p0_XGY);

	cudaMemcpyAsync(h_oLocArry, d_oLocArry, oValidFluoNum*OutParaNumGS2D * sizeof(float), cudaMemcpyDeviceToHost, cstream);

	cudaStreamSynchronize(cstream);

}


void DH3D_MoleculePair_TypeDef::Init()
{
	cudaMallocHost((void **)&h_LocArry, MaxPointNum * OutParaNumGS2D * sizeof(float));
	cudaMalloc((void **)&d_LocArry, MaxPointNum * OutParaNumGS2D * sizeof(float));

	cudaMallocHost((void **)&h_oLocArry, (MaxPointNum / 2) * OutParaNumGS2D * sizeof(float));
	cudaMalloc((void **)&d_oLocArry, (MaxPointNum / 2) * OutParaNumGS2D * sizeof(float));


	cudaMallocHost((void **)&h_PairID, (MaxPointNum / 2) * PAIR_ID_LEN * sizeof(int));
	cudaMalloc((void **)&d_PairID, (MaxPointNum / 2)  *PAIR_ID_LEN * sizeof(int));

	cudaMallocHost((void **)&h_ValidoFluoNum, sizeof(int));
	cudaMalloc((void **)&d_ValidoFluoNum, sizeof(int));

	oValidFluoNum = 0;

}

void DH3D_MoleculePair_TypeDef::Deinit()
{
	cudaError_t err;

	err = cudaFreeHost(h_LocArry);
	err = cudaFree(d_LocArry);

	err = cudaFreeHost(h_oLocArry);
	err = cudaFree(d_oLocArry);


	err = cudaFreeHost(h_PairID);
	err = cudaFree(d_PairID);

	err = cudaFreeHost(h_ValidoFluoNum);
	err = cudaFree(d_ValidoFluoNum);

}


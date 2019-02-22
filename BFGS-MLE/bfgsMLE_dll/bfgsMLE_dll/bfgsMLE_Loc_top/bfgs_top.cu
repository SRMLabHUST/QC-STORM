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

#include "bfgs_top.h"

#include "bfgsMLE_GS2D_top.h"

#include "bfgsMLE_AS3D_top.h"

#include "bfgsMLE_DH3D_top.h"


#include "bfgs_LocalizationFilter.h"
#include "bfgs_OntimeCalc.h"



void LDLocData_TypeDef::BFGS_MLELocalization(unsigned short * h_SubRegion, LocalizationPara & LocPara, int FluoNum, cudaStream_t cstream)
{
	// note there are zeros points in the localization results, there are filtered failed results,
	// they should not be classified into false positive or negative in algorithm evaluation

	cudaError_t err;

	int LocType = LocPara.LocType;
	int ROISize = LocPara.ROISize;

	int RegionDataSize = (ROISize*(ROISize + 1));

	oValidFluoNum = FluoNum;

	err = cudaMemcpyAsync(d_SubRegion, h_SubRegion, FluoNum*RegionDataSize*sizeof(short), cudaMemcpyHostToDevice, cstream);
//	HandleErr(err, "loc memcpy region");


	if (LocType_IsGS2D(LocPara.LocType))
	{
		// for 2d round Gaussian localization
		LDLoc_BFGS_MLELocalizationGS2D(d_SubRegion, d_LocArry, LocPara, FluoNum, cstream);

	}
	else if (LocType_IsAS3D(LocPara.LocType))
	{
		// for 3d astigmatism elliptical Gaussian localization
		LDLoc_BFGS_MLELocalizationAS3D(d_SubRegion, d_LocArry, LocPara, FluoNum, cstream);

	}


	// localization precision calculated by CRLB
	LocPrecCalc_GaussianCRLB(d_LocArry, LocPara, FluoNum, cstream);


	// filter bad molecules
	FilterBadFit(LocPara, FluoNum, cstream);


	// double helix localziation or commom 2d as3d
	if (LocPara.LocType == LocType_DH3D)
	{
		DH3DPairData.PairMolecules(d_LocArry, LocPara, oValidFluoNum, cstream);

		oValidFluoNum = DH3DPairData.oValidFluoNum;

		// get paired localization results
		cudaMemcpyAsync(h_LocArry, DH3DPairData.d_PairedLocArry, oValidFluoNum*OutParaNumGS2D*sizeof(float), cudaMemcpyDeviceToHost, cstream);
		cudaStreamSynchronize(cstream);

	}
	else
	{
		// 2d and AS3d
		cudaMemcpyAsync(h_LocArry, d_LocArry, oValidFluoNum*OutParaNumGS2D*sizeof(float), cudaMemcpyDeviceToHost, cstream);
		cudaStreamSynchronize(cstream);

	}

	// then call filterbadfit and OntimeCalc

}


void LDLocData_TypeDef::Init(LocalizationPara & LocPara)
{

	cudaError_t err;

	int ROIWholeLen = LocPara.ROISize*(LocPara.ROISize + 1);


	// host and gpu
	err = cudaMallocHost((void **)&h_SubRegion, MaxPointNum*ROIWholeLen*sizeof(short));
	HandleErr(err, "cudaMallocHost LDLoc h_SubRegion");

	err = cudaMalloc((void **)&d_SubRegion, MaxPointNum*ROIWholeLen*sizeof(short));
	HandleErr(err, "cudaMalloc LDLoc d_SubRegion");

	cudaMallocHost((void **)&h_LocArry, MaxPointNum*OutParaNumAS3D*sizeof(float));
	cudaMalloc((void **)&d_LocArry, MaxPointNum*OutParaNumAS3D*sizeof(float));


	// Consecutive fitting from adjecent frames

	cudaMalloc((void **)&d_ForwardLinkID, MaxPointNum*sizeof(int));
	cudaMalloc((void **)&d_BackwardLinkID, MaxPointNum*sizeof(int));
	cudaMalloc((void **)&d_ConsecutiveNum, MaxPointNum*sizeof(int));


	cudaMallocHost((void **)&h_OntimeDistrib, MaxOnTimeConsecutiveNum*sizeof(int));
	cudaMalloc((void **)&d_OntimeDistrib, MaxOnTimeConsecutiveNum*sizeof(int));

	cudaMallocHost((void **)&h_ValidFluoNum, sizeof(int));
	cudaMalloc((void **)&d_ValidFluoNum, sizeof(int));

	h_OntimeRatio = new float[MaxOnTimeConsecutiveNum];

	// double-helix 3d localization, just pair molecules after 2d localization

	if (LocPara.LocType == LocType_DH3D)
	{
		DH3DPairData.Init();
	}

	// for loc filter
	err = cudaMallocHost((void **)&h_SNRSumUp, sizeof(float));
	err = cudaMallocHost((void **)&h_ValidNum, sizeof(int));

	cudaMalloc((void **)&d_SNRSumUp, sizeof(float));
	cudaMalloc((void **)&d_ValidNum, sizeof(int));

}

void LDLocData_TypeDef::Deinit( LocalizationPara & LocPara)
{
	cudaError_t err;

	err = cudaFreeHost(h_SubRegion);
	HandleErr(err, "cudaFreeHost LDLoc h_SubRegion");
	err = cudaFree(d_SubRegion);
	HandleErr(err, "cudaFree LDLoc d_SubRegion");

	err = cudaFreeHost(h_LocArry);
	err = cudaFree(d_LocArry);

	// Consecutive fitting from adjecent frames

	cudaFree(d_ForwardLinkID);
	cudaFree(d_BackwardLinkID);
	cudaFree(d_ConsecutiveNum);

	cudaFreeHost(h_OntimeDistrib);
	cudaFree(d_OntimeDistrib);

	cudaFreeHost(h_ValidFluoNum);
	cudaFree(d_ValidFluoNum);

	delete[](h_OntimeRatio);

	// double-helix 3d localization, just pair molecules after 2d localization
	if (LocPara.LocType == LocType_DH3D)
	{
		DH3DPairData.Deinit();

	}

	// for loc filter
	cudaFreeHost(h_SNRSumUp);
	cudaFreeHost(h_ValidNum);

	cudaFree(d_SNRSumUp);
	cudaFree(d_ValidNum);

}


int LDLocData_TypeDef::GetFirstFrame(float * h_LocArry, int FluoNum)
{
	float(*pLocArry)[OutParaNumGS2D]; // for parameter array

	pLocArry = (float(*)[OutParaNumGS2D])h_LocArry;

	int cnt = 0;
	int curFrame = 0;
	for (cnt = 0; cnt < FluoNum; cnt++)
	{
		curFrame = pLocArry[cnt][Pos_Frme];

		if (curFrame != 0)break;
	}


	return curFrame;
}

int LDLocData_TypeDef::GetLastFrame(float * h_LocArry, int FluoNum)
{
	float(*pLocArry)[OutParaNumGS2D]; // for parameter array

	pLocArry = (float(*)[OutParaNumGS2D])h_LocArry;

	int cnt = 0;
	int curFrame = 0;
	for (cnt = FluoNum - 1; cnt > 0; cnt--)
	{
		curFrame = pLocArry[cnt][Pos_Frme];

		if (curFrame != 0)break;
	}


	return curFrame;
}
int LDLocData_TypeDef::GetFirstFrameFromROI(unsigned short * h_SubRegion, int ROISize, int FluoNum)
{

	int WholeROISize = ROISize*(ROISize + 1);
	int ROIPixelSize = ROISize*ROISize;

	int FrameAddrBias = ROIPixelSize + 2;
	int curFrame = h_SubRegion[FrameAddrBias];

	return curFrame;
}

int LDLocData_TypeDef::GetLastFrameFromROI(unsigned short * h_SubRegion, int ROISize, int FluoNum)
{
	int WholeROISize = ROISize*(ROISize + 1);
	int ROIPixelSize = ROISize*ROISize;

	int FrameAddrBias = WholeROISize*(FluoNum - 2) + ROIPixelSize + 2;
	int curFrame = h_SubRegion[FrameAddrBias];

	return curFrame;
}


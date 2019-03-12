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

#include "bfgsMLE_core.h"


#include "bfgs_LocalizationFilter.h"
#include "bfgs_OntimeCalc.h"

#include "LDROIExtraction.h"


void LDLocData_TypeDef::BFGS_MLELocalization(unsigned short * h_ImageROI, float *h_WLEPara, LocalizationPara & LocPara, int FluoNum, cudaStream_t cstream)
{
	// note there are zeros points in the localization results, there are filtered failed results,
	// they should not be classified into false positive or negative in algorithm evaluation

	cudaError_t err;

	int LocType = LocPara.LocType;
	int ROISize = LocPara.ROISize;

	int ROIWholeSize = (ROISize*(ROISize + 1));

	oValidFluoNum = FluoNum;

	err = cudaMemcpyAsync(d_ImageROI, h_ImageROI, FluoNum * ROIWholeSize *sizeof(short), cudaMemcpyHostToDevice, cstream);

	if (WLE_ENABLE)
	{
		cudaMemcpyAsync(d_WLEPara, h_WLEPara, FluoNum * WLE_ParaNumber * sizeof(float), cudaMemcpyHostToDevice, cstream);
	}

	// for multi emitter fitting
	cudaMemsetAsync(d_MultiFitFluoNum, 0, sizeof(int), cstream);


	// Low density fitting
	if (LocPara.LocType == LocType_GS2D)
	{
		// for 2d round Gaussian localization
		LDLoc_BFGS_MLELocalizationGS2D(d_LocArry, d_ImageROI, d_WLEPara, d_MultiFitFluoNum, d_MultiFitFluoPos, LocPara, FluoNum, cstream);
	
	}
	else if (LocPara.LocType == LocType_AS3D)
	{
		// for 3d astigmatism elliptical Gaussian localization
		LDLoc_BFGS_MLELocalizationAS3D(d_LocArry, d_ImageROI, d_WLEPara, d_MultiFitFluoNum, d_MultiFitFluoPos, LocPara, FluoNum, cstream);
	
	}

	// high density fitting
	if (LocPara.MultiEmitterFitEn)
	{
		cudaMemcpyAsync(h_MultiFitFluoNum, d_MultiFitFluoNum, sizeof(int), cudaMemcpyDeviceToHost, cstream);
		cudaStreamSynchronize(cstream);

		int MultiFitFluoNum = *h_MultiFitFluoNum;

		MultiFitRatio = (float)MultiFitFluoNum / FluoNum;

		if (LocPara.LocType == LocType_GS2D)
		{
			HDLoc_BFGS_MLELocalization_2D_2Emitter(d_LocArry, d_ImageROI, d_WLEPara, MultiFitFluoNum, d_MultiFitFluoPos, LocPara, cstream);
		
		}
		else if (LocPara.LocType == LocType_AS3D)
		{
			HDLoc_BFGS_MLELocalization_AS3D_2Emitter(d_LocArry, d_ImageROI, d_WLEPara, MultiFitFluoNum, d_MultiFitFluoPos, LocPara, cstream);
		
		}

		cudaStreamSynchronize(cstream);
	}




	// localization precision calculated by CRLB
	LocPrecCalc_GaussianCRLB(d_LocArry, LocPara, FluoNum, cstream);


	// filter bad molecules
	FilterBadFit(LocPara, FluoNum, cstream);


	// 2d and AS3d
	cudaMemcpyAsync(h_LocArry, d_LocArry, oValidFluoNum*OutParaNumGS2D * sizeof(float), cudaMemcpyDeviceToHost, cstream);
	cudaStreamSynchronize(cstream);


}


void LDLocData_TypeDef::Init(LocalizationPara & LocPara)
{

	cudaError_t err;

	int ROIWholeSize = LocPara.ROISize*(LocPara.ROISize + 1);


	// host and gpu
	err = cudaMallocHost((void **)&h_ImageROI, MaxPointNum*ROIWholeSize*sizeof(short));
	HandleErr(err, "cudaMallocHost LDLoc h_ImageROI");

	err = cudaMalloc((void **)&d_ImageROI, MaxPointNum*ROIWholeSize*sizeof(short));
	HandleErr(err, "cudaMalloc LDLoc d_ImageROI");

	cudaMallocHost((void **)&h_LocArry, MaxPointNum*OutParaNumAS3D*sizeof(float));
	cudaMalloc((void **)&d_LocArry, MaxPointNum*OutParaNumAS3D*sizeof(float));


	err = cudaMalloc((void **)&d_WLEPara, MaxPointNum * WLE_ParaNumber * sizeof(float));


	// Consecutive fitting from adjecent frames

	cudaMalloc((void **)&d_ForwardLinkID, MaxPointNum*sizeof(int));
	cudaMalloc((void **)&d_BackwardLinkID, MaxPointNum*sizeof(int));
	cudaMalloc((void **)&d_ConsecutiveNum, MaxPointNum*sizeof(int));


	cudaMallocHost((void **)&h_OntimeDistrib, MaxOnTimeConsecutiveNum*sizeof(int));
	cudaMalloc((void **)&d_OntimeDistrib, MaxOnTimeConsecutiveNum*sizeof(int));

	cudaMallocHost((void **)&h_ValidFluoNum, sizeof(int));
	cudaMalloc((void **)&d_ValidFluoNum, sizeof(int));

	h_OntimeRatio = new float[MaxOnTimeConsecutiveNum];


	// for loc filter
	err = cudaMallocHost((void **)&h_SNRSumUp, sizeof(float));
	err = cudaMallocHost((void **)&h_ValidNum, sizeof(int));

	cudaMalloc((void **)&d_SNRSumUp, sizeof(float));
	cudaMalloc((void **)&d_ValidNum, sizeof(int));


	// multi emitter fitting
	cudaMallocHost((void **)&h_MultiFitFluoNum, sizeof(int));
	cudaMalloc((void **)&d_MultiFitFluoNum, sizeof(int));

	cudaMalloc((void **)&d_MultiFitFluoPos, MaxPointNum * sizeof(int));

	MultiFitRatio = 0;

}

void LDLocData_TypeDef::Deinit( LocalizationPara & LocPara)
{
	cudaError_t err;

	err = cudaFreeHost(h_ImageROI);
	HandleErr(err, "cudaFreeHost LDLoc h_ImageROI");
	err = cudaFree(d_ImageROI);
	HandleErr(err, "cudaFree LDLoc d_ImageROI");

	err = cudaFreeHost(h_LocArry);
	err = cudaFree(d_LocArry);

	err = cudaFree(d_WLEPara);
	// Consecutive fitting from adjecent frames

	cudaFree(d_ForwardLinkID);
	cudaFree(d_BackwardLinkID);
	cudaFree(d_ConsecutiveNum);

	cudaFreeHost(h_OntimeDistrib);
	cudaFree(d_OntimeDistrib);

	cudaFreeHost(h_ValidFluoNum);
	cudaFree(d_ValidFluoNum);

	delete[](h_OntimeRatio);


	// for loc filter
	cudaFreeHost(h_SNRSumUp);
	cudaFreeHost(h_ValidNum);

	cudaFree(d_SNRSumUp);
	cudaFree(d_ValidNum);



	// multi emitter fitting
	cudaFreeHost(h_MultiFitFluoNum);

	cudaFree(d_MultiFitFluoNum);
	cudaFree(d_MultiFitFluoPos);


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
int LDLocData_TypeDef::GetFirstFrameFromROI(unsigned short * h_ImageROI, int ROISize, int FluoNum)
{

	int WholeROISize = ROISize*(ROISize + 1);
	int ROIPixelSize = ROISize*ROISize;

	int FrameAddrBias = ROIPixelSize + 2;
	int curFrame = h_ImageROI[FrameAddrBias];

	return curFrame;
}

int LDLocData_TypeDef::GetLastFrameFromROI(unsigned short * h_ImageROI, int ROISize, int FluoNum)
{
	int WholeROISize = ROISize*(ROISize + 1);
	int ROIPixelSize = ROISize*ROISize;

	int FrameAddrBias = WholeROISize*(FluoNum - 2) + ROIPixelSize + 2;
	int curFrame = h_ImageROI[FrameAddrBias];

	return curFrame;
}


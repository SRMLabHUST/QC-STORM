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


#include "LDROIExtraction.h"

#include <time.h>


void LDLocData_TypeDef::BFGS_MLELocalization(unsigned short * h_ImageROI, float *h_WLEPara, LocalizationPara & LocPara, int FluoNum, cudaStream_t cstream)
{

	CopyFittingPara(LocPara, cstream);


	// initial for multi emitter fitting, filled by single molecule localization
	ResetNumbers(cstream);

	if (FluoNum <= 0)return;
	if (FluoNum > MaxPointNum)FluoNum = MaxPointNum;


	int ROIWholeSize = (LocPara.ROISize*(LocPara.ROISize + 1));


	cudaMemcpyAsync(d_ImageROI, h_ImageROI, FluoNum * ROIWholeSize *sizeof(short), cudaMemcpyHostToDevice, cstream);

	if (WLE_ENABLE)
	{
		cudaMemcpyAsync(d_WLEPara, h_WLEPara, FluoNum * WLE_ParaNumber * sizeof(float), cudaMemcpyHostToDevice, cstream);
	}



	// initial molecule add position for multi emitter fitting
	*(h_FitPosInf->h_MultiFit_AddedFluoNum) = FluoNum;
	cudaMemcpyAsync(h_FitPosInf->d_MultiFit_AddedFluoNum, h_FitPosInf->h_MultiFit_AddedFluoNum, sizeof(int), cudaMemcpyHostToDevice, cstream);

	
	// clasify molecule by preestimated data
	MoleculePreFitClasify(LocPara.ROISize, LocPara.MultiEmitterFitEn, FluoNum, cstream);


	int SingleFitNum = *(h_FitPosInf->h_SingleFitFluoNum);

	int E2FitNum_o = *(h_FitPosInf->h_MultiFitFluoNum_2E);

	int FluoNum_div = FluoNum;
	if (FluoNum_div < 1)FluoNum_div = 1;


	FitRatio_1E = (float)SingleFitNum / FluoNum_div;

	// will be updated later for multi-emitter fitting
	FitRatio_Final_1E = (float)(SingleFitNum) / FluoNum_div;


	// Low density fitting
	if ((LocPara.LocType == LocType_GS2D) ||(LocPara.LocType == LocType_DH3D))
	{
		// for 2d round Gaussian localization
		LDLoc_BFGS_MLELocalizationGS2D(d_LocArry, d_ImageROI, d_WLEPara, SingleFitNum, d_FitPosInf, d_FitPara, LocPara.ROISize, cstream);

	}
	else if (LocPara.LocType == LocType_AS3D)
	{
		// for 3d astigmatism elliptical Gaussian localization
		LDLoc_BFGS_MLELocalizationAS3D(d_LocArry, d_ImageROI, d_WLEPara, SingleFitNum, d_FitPosInf, d_FitPara, LocPara.ROISize, cstream);

	}

	oValidFluoNum = FluoNum;


	// high density fitting
	if (LocPara.MultiEmitterFitEn)
	{

#if(SPEED_TEST) // only for speed testing
		int time1 = clock();

#endif // SPEED_TEST


		// two emitter fitting
		cudaMemcpyAsync(h_FitPosInf->h_MultiFitFluoNum_2E, h_FitPosInf->d_MultiFitFluoNum_2E, sizeof(int), cudaMemcpyDeviceToHost, cstream);
		cudaStreamSynchronize(cstream);

		int MultiFitFluoNum_2E = *(h_FitPosInf->h_MultiFitFluoNum_2E);



		if (LocPara.LocType == LocType_GS2D)
		{
			HDLoc_BFGS_MLELocalization_2D_2Emitter(d_LocArry, d_ImageROI, d_WLEPara, MultiFitFluoNum_2E, d_FitPosInf, d_FitPara, LocPara.ROISize, cstream);
	
		}
		else if (LocPara.LocType == LocType_AS3D)
		{
			HDLoc_BFGS_MLELocalization_AS3D_2Emitter(d_LocArry, d_ImageROI, d_WLEPara, MultiFitFluoNum_2E, d_FitPosInf, d_FitPara, LocPara.ROISize, cstream);
	
		}

		cudaStreamSynchronize(cstream);


#if(SPEED_TEST) // only for speed testing
		int time2 = clock();

#endif // SPEED_TEST


		// three emitter fitting
		cudaMemcpyAsync(h_FitPosInf->h_MultiFitFluoNum_3E, h_FitPosInf->d_MultiFitFluoNum_3E, sizeof(int), cudaMemcpyDeviceToHost, cstream);
		cudaStreamSynchronize(cstream);

		int MultiFitFluoNum_3E = *(h_FitPosInf->h_MultiFitFluoNum_3E);


		if (LocPara.LocType == LocType_GS2D)
		{
			HDLoc_BFGS_MLELocalization_2D_3Emitter(d_LocArry, d_ImageROI, d_WLEPara, MultiFitFluoNum_3E, d_FitPosInf, d_FitPara, LocPara.ROISize, cstream);
		}


		// added by multi emitter fitting
		cudaMemcpyAsync(h_FitPosInf->h_MultiFit_AddedFluoNum, h_FitPosInf->d_MultiFit_AddedFluoNum, sizeof(int), cudaMemcpyDeviceToHost, cstream);

		cudaMemcpyAsync(h_FitPosInf->h_RejectedFluoNum, h_FitPosInf->d_RejectedFluoNum, sizeof(int), cudaMemcpyDeviceToHost, cstream);

		cudaStreamSynchronize(cstream);


		oValidFluoNum = *(h_FitPosInf->h_MultiFit_AddedFluoNum);

		int RejectedFluoNum = *(h_FitPosInf->h_RejectedFluoNum);


		FitRatio_2E = (float)MultiFitFluoNum_2E / FluoNum_div;
		FitRatio_3E = (float)MultiFitFluoNum_3E / FluoNum_div;


		//
		FitRatio_Final_1E = (float)(SingleFitNum - (MultiFitFluoNum_2E - E2FitNum_o)) / FluoNum_div;

		if (LocPara.LocType == LocType_GS2D)
		{
			FitRatio_Final_2E = (float)(MultiFitFluoNum_2E - MultiFitFluoNum_3E) / FluoNum_div;
			FitRatio_Final_3E = (float)(MultiFitFluoNum_3E - RejectedFluoNum) / FluoNum_div;

		}
		else
		{
			FitRatio_Final_2E = (float)(MultiFitFluoNum_2E - RejectedFluoNum) / FluoNum_div;
		}

		FitRatio_Final_4E = (float)RejectedFluoNum / FluoNum_div;

//		printf("FitRatio_1E-4E: %.4f %.4f %.4f %.4f\n", FitRatio_Final_1E, FitRatio_Final_2E, FitRatio_Final_3E, FitRatio_Final_4E);


#if(SPEED_TEST) // only for speed testing
		int time3 = clock();

		int time_2f = time2 - time1;
		int time_3f = time3 - time2;
		if (time_2f < 1)time_2f = 1;
		if (time_3f < 1)time_3f = 1;
		int LocSpeed_2f = MultiFitFluoNum_2E * 1000 / time_2f;
		int LocSpeed_3f = MultiFitFluoNum_3E * 1000 / time_3f;

		printf("MultiFitFluoNum_2E:%d - speed %d /s\n", MultiFitFluoNum_2E, LocSpeed_2f);
		printf("MultiFitFluoNum_3E:%d - speed %d /s\n", MultiFitFluoNum_3E, LocSpeed_3f);

#endif // SPEED_TEST

	}


	// localization precision calculated by CRLB
	LocPrecCalc_GaussianCRLB(d_LocArry, LocPara, oValidFluoNum, cstream);


	// 2d and AS3d
	cudaMemcpyAsync(h_LocArry, d_LocArry, oValidFluoNum * OutParaNumGS2D * sizeof(float), cudaMemcpyDeviceToHost, cstream);
	cudaStreamSynchronize(cstream);

}

void LDLocData_TypeDef::ResetNumbers(cudaStream_t cstream)
{

	oValidFluoNum = 0;


	FitRatio_1E = 0; // single molecule fit ratio
	FitRatio_2E = 0; // two emitter fit ratio
	FitRatio_3E = 0; // three emitter fit ratio

	FitRatio_Final_1E = 0; // single molecule fitting ratio
	FitRatio_Final_2E = 0; // two emitter fitting ratio
	FitRatio_Final_3E = 0; // three emitter fitting ratio
	FitRatio_Final_4E = 0; // four or more emitter fitting ratio


	//

	cudaMemsetAsync(h_FitPosInf->d_SingleFitFluoNum, 0, sizeof(int), cstream);

	cudaMemsetAsync(h_FitPosInf->d_MultiFitFluoNum_2E, 0, sizeof(int), cstream);
	cudaMemsetAsync(h_FitPosInf->d_MultiFitFluoNum_3E, 0, sizeof(int), cstream);
	cudaMemsetAsync(h_FitPosInf->d_RejectedFluoNum, 0, sizeof(int), cstream);

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


	err = cudaMalloc((void **)&d_WLEPara, MaxPointNum * WLE_ParaNumber * sizeof(float));


	cudaMallocHost((void **)&h_LocArry, MaxPointNum*OutParaNumAS3D*sizeof(float));
	cudaMalloc((void **)&d_LocArry, MaxPointNum*OutParaNumAS3D*sizeof(float));


	// for convinient parameter transimission
	cudaMallocHost((void **)&h_FitPara, sizeof(CoreFittingPara));
	cudaMalloc((void **)&d_FitPara, sizeof(CoreFittingPara));

	cudaMallocHost((void **)&h_FitPosInf, sizeof(FitPosInf_TypeDef));
	cudaMalloc((void **)&d_FitPosInf, sizeof(FitPosInf_TypeDef));


	// single molecule fitting

	cudaMallocHost((void **)&h_FitPosInf->h_SingleFitFluoNum, sizeof(int));
	cudaMalloc((void **)&h_FitPosInf->d_SingleFitFluoNum, sizeof(int));
	cudaMalloc((void **)&h_FitPosInf->d_SingleFitFluoPos, MaxPointNum * sizeof(int));


	// multi emitter fitting
	cudaMallocHost((void **)&h_FitPosInf->h_MultiFitFluoNum_2E, sizeof(int));
	cudaMalloc((void **)&h_FitPosInf->d_MultiFitFluoNum_2E, sizeof(int));
	cudaMalloc((void **)&h_FitPosInf->d_MultiFitFluoPos_2E, MaxPointNum * sizeof(int));

	cudaMallocHost((void **)&h_FitPosInf->h_MultiFitFluoNum_3E, sizeof(int));
	cudaMalloc((void **)&h_FitPosInf->d_MultiFitFluoNum_3E, sizeof(int));
	cudaMalloc((void **)&h_FitPosInf->d_MultiFitFluoPos_3E, MaxPointNum * sizeof(int));

	cudaMallocHost((void **)&h_FitPosInf->h_RejectedFluoNum, sizeof(int));
	cudaMalloc((void **)&h_FitPosInf->d_RejectedFluoNum, sizeof(int));


	cudaMallocHost((void **)&h_FitPosInf->h_MultiFit_AddedFluoNum, sizeof(int));
	cudaMalloc((void **)&h_FitPosInf->d_MultiFit_AddedFluoNum, sizeof(int));


	// copy assigned memory for GPU use
	cudaMemcpy(d_FitPosInf, h_FitPosInf, sizeof(FitPosInf_TypeDef), cudaMemcpyHostToDevice);


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



	// single molecule fitting
	cudaFreeHost(h_FitPosInf->h_SingleFitFluoNum);
	cudaFree(h_FitPosInf->d_SingleFitFluoNum);
	cudaFree(h_FitPosInf->d_SingleFitFluoPos);


	// multi emitter fitting
	cudaFreeHost(h_FitPosInf->h_MultiFitFluoNum_2E);
	cudaFree(h_FitPosInf->d_MultiFitFluoNum_2E);
	cudaFree(h_FitPosInf->d_MultiFitFluoPos_2E);

	cudaFreeHost(h_FitPosInf->h_MultiFitFluoNum_3E);
	cudaFree(h_FitPosInf->d_MultiFitFluoNum_3E);
	cudaFree(h_FitPosInf->d_MultiFitFluoPos_3E);

	cudaFreeHost(h_FitPosInf->h_RejectedFluoNum);
	cudaFree(h_FitPosInf->d_RejectedFluoNum);


	cudaFreeHost(h_FitPosInf->h_MultiFit_AddedFluoNum);
	cudaFree(h_FitPosInf->d_MultiFit_AddedFluoNum);


	// must last release
	err = cudaFreeHost(h_FitPara);
	err = cudaFree(d_FitPara);

	err = cudaFreeHost(h_FitPosInf);
	err = cudaFree(d_FitPosInf);

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

void LDLocData_TypeDef::CopyFittingPara(LocalizationPara & LocPara, cudaStream_t cstream)
{
	h_FitPara->MultiEmitterFitEn = LocPara.MultiEmitterFitEn;
	h_FitPara->WLEEn = LocPara.WLEEn;

	// camera parameter
	h_FitPara->Offset = LocPara.Offset;
	h_FitPara->KAdc = LocPara.KAdc;
	h_FitPara->QE = LocPara.QE;
	h_FitPara->ReadNoise_e = LocPara.ReadNoise_e;

	// 3D imaging
	h_FitPara->ZDepthCorrFactor = LocPara.ZDepthCorrFactor;

	// calibration of sigma X >= sigma Y
	h_FitPara->p4_XGY = LocPara.p4_XGY;
	h_FitPara->p3_XGY = LocPara.p3_XGY;
	h_FitPara->p2_XGY = LocPara.p2_XGY;
	h_FitPara->p1_XGY = LocPara.p1_XGY;
	h_FitPara->p0_XGY = LocPara.p0_XGY;

	// calibration of sigma X < sigma Y
	h_FitPara->p4_XLY = LocPara.p4_XLY;
	h_FitPara->p3_XLY = LocPara.p3_XLY;
	h_FitPara->p2_XLY = LocPara.p2_XLY;
	h_FitPara->p1_XLY = LocPara.p1_XLY;
	h_FitPara->p0_XLY = LocPara.p0_XLY;


	cudaMemcpyAsync(d_FitPara, h_FitPara, sizeof(CoreFittingPara), cudaMemcpyHostToDevice, cstream);

}

/*
// estimage activation density by ratio of single-molecule fit ROI
// is replaced by molecule ratio in FluoStatisticData_TypeDef, combining the two methods shold be optimal
float LDLocData_TypeDef::GetDensityEstimation(float FitRatio_Final_1E, int ImagingMode)
{
	// LocType_GS2D
	float p1 =   0.0006503f;
	float p2 =     -0.1201f;
	float p3 =       5.771f;


	if(ImagingMode == LocType_AS3D)
	{
		p1 =   0.0006733f;
		p2 =     -0.1229f;
		p3 =       5.762f;
	}

	float ActivationDensity = p1*FitRatio_Final_1E*FitRatio_Final_1E + p2*FitRatio_Final_1E + p3;
	return ActivationDensity;
}

*/

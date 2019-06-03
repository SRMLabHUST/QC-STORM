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

#pragma once

#include "bfgsMLE_core.h"

#include "WLEParaEstimation_Parameters.h"


// MLE fitting iteration number 

#define IterateNum_2D			8  // total iteration number, fixed iteration number
#define IterateNum_2D_bs		11 // bisection iteration to find best walk length


// number of fitting parameter
#define FitParaNum_2D			5

// fitting parameters ID
#define Fit2D_Peak				0
#define Fit2D_XPos				1
#define Fit2D_YPos				2
#define Fit2D_Sigm				3
#define Fit2D_Bakg				4


template <int ROISize, int ROIPixelNum, int FitParaNum, int IterateNum, int IterateNum_bs>
__device__ void MLELocalization_GS2D(float ImageROI[][ROIPixelNum], float Ininf[][ThreadsPerBlock], float grad[][ThreadsPerBlock], float d0[][ThreadsPerBlock], float D0[], float* WLE_Weight, int tid);


template <int ROISize, int ROIPixelNum, int FitParaNum>
__device__ void poissonfGradient_GS2D(float ImageROI[][ROIPixelNum], float Ininf[][ThreadsPerBlock], float grad[][ThreadsPerBlock], float* WLE_Weight, int tid);


template <int ROISize, int ROIPixelNum, int FitParaNum>
__device__ float MLEFit_TargetF_GS2D(float ImageROI[][ROIPixelNum], float Ininf[], float* WLE_Weight, int tid);


template <int ROISize, int ROIPixelNum, int FitParaNum>
__device__ void PreEstimation_GS2D(float ImageROI[][ROIPixelNum], float Ininf[][ThreadsPerBlock], float WLE_SigmaX, int tid);


template <int ROISize, int FitParaNum>
__global__ void bfgsMLELoc_Gauss2D(float *d_LocArry, unsigned short *d_ImageROI, float *d_WLEPara, int SingleFitNum, FitPosInf_TypeDef* d_FitPosInf, CoreFittingPara *d_FitPara)
{

	enum {
		ROIPixelNum = (ROISize*ROISize),
		ROIWholeSize = (ROISize*(ROISize + 1)),

		ROISize_Half = (int)(ROISize / 2),

		LoadLoopNum = ((ROIPixelNum + ThreadsPerBlock - 1) / ThreadsPerBlock)

	};


#define XYPosVaryBound		1.5f

	float XYLBound = ROISize / 2.0f - XYPosVaryBound; // 1.0f
	float XYUBound = ROISize / 2.0f + XYPosVaryBound;

	float QE = d_FitPara->QE;
	int MultiEmitterFitEn = d_FitPara->MultiEmitterFitEn;

	
	int *d_SingleFitFluoPos = d_FitPosInf->d_SingleFitFluoPos;
	int * d_MultiFitFluoNum_2E = d_FitPosInf->d_MultiFitFluoNum_2E;
	int * d_MultiFitFluoPos_2E = d_FitPosInf->d_MultiFitFluoPos_2E;


	float Simga_Max_Th = ROISize / 1.5f / 2.35f;

	if (MultiEmitterFitEn)
	{
		Simga_Max_Th = ROISize / 2.0f / 2.35f;
	}


	//
	__shared__ float ImageROI[ThreadsPerBlock][ROIPixelNum];
	__shared__ float D0[ThreadsPerBlock][FitParaNum*FitParaNum];	// inv of matrix hessian 

	// avoid bank conflict
	__shared__ float Ininf[FitParaNum][ThreadsPerBlock]; // cur position
	__shared__ float grad[FitParaNum][ThreadsPerBlock];  // gradient
	__shared__ float d0[FitParaNum][ThreadsPerBlock];	// direction


	int gid = threadIdx.x + blockDim.x*blockIdx.x;
	int gid_0 = blockDim.x*blockIdx.x;

	int tid = threadIdx.x;


	float(*pLocArry)[OutParaNumGS2D]; // for parameter array
	pLocArry = (float(*)[OutParaNumGS2D])d_LocArry;


	float(*pD0)[FitParaNum] = (float(*)[FitParaNum])&D0[tid][0];

	// WLE parameter array
	float(*pWLEPara)[WLE_ParaNumber] = (float(*)[WLE_ParaNumber])d_WLEPara;

	float WLE_Weight[ROIPixelNum];



	// load image region from global memory
	for (int fcnt = 0; fcnt < ThreadsPerBlock; fcnt++)
	{
		int CurId = gid_0 + fcnt;

		CurId = min(CurId, SingleFitNum - 1);

		// assign correct position for detected molecules
		int CurMapID = CurId;
		if (MultiEmitterFitEn) CurMapID = d_SingleFitFluoPos[CurId];


		int AddrOffset = CurMapID*ROIWholeSize;

#pragma unroll
		for (int cnt = 0; cnt < LoadLoopNum; cnt++)
		{
			int CurLoadId = cnt*ThreadsPerBlock + tid;

			if (CurLoadId < ROIPixelNum)
			{
				ImageROI[fcnt][CurLoadId] = fmaxf((d_ImageROI[AddrOffset + CurLoadId] - d_FitPara->Offset)*d_FitPara->KAdc, 1.0f);
			}
		}
	}


	int rcnt, ccnt;

	//initial D0
#pragma unroll
	for (rcnt = 0; rcnt < FitParaNum; rcnt++)
	{
#pragma unroll
		for (ccnt = 0; ccnt < FitParaNum; ccnt++)
		{
			if (rcnt == ccnt)pD0[rcnt][ccnt] = 1.0f;
			else pD0[rcnt][ccnt] = 0.0f;
		}
	}
	

	// why this function impact the time so much? it's should be placed here
	__syncthreads();



	if (gid < SingleFitNum)
	{
		// assign correct position for detected molecules
		int CurFluoID = gid;
		if (MultiEmitterFitEn) CurFluoID = d_SingleFitFluoPos[gid];
		

		int XOffset = d_ImageROI[CurFluoID*ROIWholeSize + ROIPixelNum + 0];
		int YOffset = d_ImageROI[CurFluoID*ROIWholeSize + ROIPixelNum + 1];
		int CurFrame = d_ImageROI[CurFluoID*ROIWholeSize + ROIPixelNum + 2] + 65536 * d_ImageROI[CurFluoID*ROIWholeSize + ROIPixelNum + 3];


		// estimated PSF width for mask of WLE, weighted MLE
		int MoleculeType = 0;

		float WLE_SigmaX = 0;


#if(WLE_ENABLE == 1)

		float SigmaL = pWLEPara[CurFluoID][WLE_Para_SigmaL];
		float SigmaR = pWLEPara[CurFluoID][WLE_Para_SigmaR];
		float SigmaU = pWLEPara[CurFluoID][WLE_Para_SigmaU];
		float SigmaD = pWLEPara[CurFluoID][WLE_Para_SigmaD];

		MoleculeType = pWLEPara[CurFluoID][WLE_Para_FluoType];

		SigmaL = min(SigmaL, SigmaR);
		SigmaU = min(SigmaU, SigmaD);

		WLE_SigmaX = min(SigmaL, SigmaU);

		WLE_SigmaX = min(WLE_SigmaX, ROISize / 2.0f / 2.35f);

#endif // WLE_ENABLE


		// pre-estimation
		PreEstimation_GS2D<ROISize, ROIPixelNum, FitParaNum>(ImageROI, Ininf, WLE_SigmaX, tid);


#if(WLE_ENABLE == 1)

		if (MoleculeType != MoleculeType_MLEFit)
		{
			// smaller Width for contaminated molecules
			WLE_SigmaX = WLE_SigmaX / 1.2f;
			WLE_SigmaX = min(WLE_SigmaX, ROISize / 2.0f / 2.35f);
		}

		if (d_FitPara->WLEEn == 0)
		{
			if (MoleculeType != MoleculeType_MultiFit)
			{
				WLE_SigmaX = ROISize;
			}
		}

		WLEWeightsCalc<ROISize>(WLE_Weight, WLE_SigmaX, WLE_SigmaX, MoleculeType);

#endif // WLE_ENABLE


		MLELocalization_GS2D<ROISize, ROIPixelNum, FitParaNum_2D, IterateNum_2D, IterateNum_2D_bs>(ImageROI, Ininf, grad, d0, &D0[tid][0], WLE_Weight, tid);


		// remove scaling factor
		float PeakPhoton = Ininf[Fit2D_Peak][tid] * AScalFactor / QE; // PeakPhoton
		float Background = Ininf[Fit2D_Bakg][tid] * BScalFactor / QE; // bkg

		float FittedXPos = Ininf[Fit2D_XPos][tid];
		float FittedYPos = Ininf[Fit2D_YPos][tid];

		float SimgaX = Ininf[Fit2D_Sigm][tid] * SScalFactor;

		SimgaX = sqrtf(0.5f / SimgaX); // pixel

		// remove nan
		FittedXPos = isnan(FittedXPos) ? -1 : FittedXPos;
		FittedYPos = isnan(FittedYPos) ? -1 : FittedYPos;
		SimgaX = isnan(SimgaX) ? -1 : SimgaX;

		float CurSNR = sqrtf(QE)*PeakPhoton / sqrtf(PeakPhoton + Background);

		float TotalPhoton = PeakPhoton * 2 * 3.141593f*SimgaX*SimgaX;

		float XPos = FittedXPos - ROISize_Half + XOffset; // X
		float YPos = FittedYPos - ROISize_Half + YOffset; // Y

		// remove nan
		int IsFittingValid = isnan(Background) ? 0 : 1;


		if (IsFittingValid && (CurSNR > 4.0f) && (SimgaX <= Simga_Max_Th) && (FittedXPos > XYLBound) && (FittedXPos < XYUBound) && (FittedYPos > XYLBound) && (FittedYPos < XYUBound)) //
		{
			pLocArry[CurFluoID][Pos_PPho] = PeakPhoton; // peak photon
			pLocArry[CurFluoID][Pos_XPos] = XPos; // may have 0.5 or 1 pixel offset compared with other software
			pLocArry[CurFluoID][Pos_YPos] = YPos; // may have 0.5 or 1 pixel offset compared with other software
			pLocArry[CurFluoID][Pos_ZPos] = 0.0f; // 
			pLocArry[CurFluoID][Pos_SigX] = SimgaX; // sigma x
			pLocArry[CurFluoID][Pos_SigY] = SimgaX; // sigma y
			pLocArry[CurFluoID][Pos_TPho] = TotalPhoton; // total photon
			pLocArry[CurFluoID][Pos_Bakg] = Background; // background
			pLocArry[CurFluoID][Pos_PSNR] = CurSNR; // peak snr
			pLocArry[CurFluoID][Pos_Frme] = CurFrame; // frame

		}
		else 
		{
			if (MultiEmitterFitEn)
			{
				// mark position for multi emitter fitting
				if (IsFittingValid && (CurSNR > 6.0f))
				{

#if(!SPEED_TEST)

					// mark position for multi emitter fitting
					int StorePos = atomicAdd(d_MultiFitFluoNum_2E, 1);

					if (StorePos < MaxPointNum)
					{
						d_MultiFitFluoPos_2E[StorePos] = CurFluoID;
					}

#endif // SPEED_TEST

				}
			}


			// remove invalid point
#pragma unroll
			for (rcnt = 0; rcnt < OutParaNumGS2D; rcnt++)
			{
				pLocArry[CurFluoID][rcnt] = 0;

			}
		}
	}
}


template <int ROISize, int ROIPixelNum, int FitParaNum, int IterateNum, int IterateNum_bs>
__device__ void MLELocalization_GS2D(float ImageROI[][ROIPixelNum], float Ininf[][ThreadsPerBlock], float grad[][ThreadsPerBlock], float d0[][ThreadsPerBlock], float D0[], float* WLE_Weight, int tid)
{
	// adjust d0
	float td0[FitParaNum];
	float td0_total;

	float sk[FitParaNum]; // bfgs quasi-newton method
	float yk[FitParaNum];
	float tgrad[FitParaNum];


	float xd[FitParaNum * 2];

	float ddat[2];
	float dpos[2];


	// initialize
	poissonfGradient_GS2D<ROISize, ROIPixelNum, FitParaNum>(ImageROI, Ininf, grad, WLE_Weight, tid);

#pragma unroll
	for (int rcnt = 0; rcnt < FitParaNum; rcnt++)
	{
		d0[rcnt][tid] = -grad[rcnt][tid];
	}


	// fitting iteration
	for (int itcnt = 0; itcnt<IterateNum; itcnt++)
	{
		td0_total = 0;

		// adjust d0
#pragma unroll
		for (int rcnt = 0; rcnt < FitParaNum; rcnt++)
		{
			td0[rcnt] = abs(d0[rcnt][tid]);
			td0_total += td0[rcnt];
		}
		// reciprocal of mean
		td0_total = ((float)FitParaNum) / td0_total;


#pragma unroll
		for (int rcnt = 0; rcnt < FitParaNum; rcnt++)
		{
			d0[rcnt][tid] = d0[rcnt][tid] * td0_total;

		}


		dpos[0] = 0.00001f; // scale factor left limit, should not equal to 0 and smaller
		dpos[1] = 1.0f; // scale factor right limit, should not lager than 2

		VectorAddMul1<FitParaNum>(&xd[0], Ininf, d0, 0.0001f, tid);
		VectorAddMul1<FitParaNum>(&xd[FitParaNum], Ininf, d0, 1.0f, tid);

		ddat[0] = MLEFit_TargetF_GS2D<ROISize, ROIPixelNum, FitParaNum>(ImageROI, &xd[0], WLE_Weight, tid);
		ddat[1] = MLEFit_TargetF_GS2D<ROISize, ROIPixelNum, FitParaNum>(ImageROI, &xd[FitParaNum], WLE_Weight, tid);

		// bisection method to find optimal walk length
		for (int cnt = 0; cnt < IterateNum_bs; cnt++)
		{
			// which part shrink
			int xdsel = (ddat[0] < ddat[1]);


			dpos[xdsel] = (dpos[0] + dpos[1])*0.5f; //  /2.0f which one shift

			if (cnt < IterateNum_bs - 1)
			{
				VectorAddMul1<FitParaNum>(&xd[xdsel*FitParaNum], Ininf, d0, dpos[xdsel], tid);// xd=ininf+d0*scale
				ddat[xdsel] = MLEFit_TargetF_GS2D<ROISize, ROIPixelNum, FitParaNum>(ImageROI, &xd[xdsel*FitParaNum], WLE_Weight, tid);
			}

		}
		float scale = (dpos[0] + dpos[1])*0.5f; //  /2.0f calculated direction

		// sk 		= d0*base;
		// inInf 	= inInf+d0*base;
#pragma unroll
		for (int rcnt = 0; rcnt < FitParaNum; rcnt++)
		{
			sk[rcnt] = d0[rcnt][tid] * scale;
			Ininf[rcnt][tid] = Ininf[rcnt][tid] + sk[rcnt];

		}

		IninfConstrain<FitParaNum>(Ininf, tid); // 


		if (itcnt < IterateNum - 1)
		{
			// tgrad = grad;
#pragma unroll
			for (int rcnt = 0; rcnt < FitParaNum; rcnt++)
			{
				tgrad[rcnt] = grad[rcnt][tid];
			}

			poissonfGradient_GS2D<ROISize, ROIPixelNum, FitParaNum>(ImageROI, Ininf, grad, WLE_Weight, tid);

#pragma unroll
			for (int rcnt = 0; rcnt < FitParaNum; rcnt++)
			{
				yk[rcnt] = grad[rcnt][tid] - tgrad[rcnt];
			}

			ConstructD0<FitParaNum>(D0, sk, yk, tid);
			MatMulVector<FitParaNum>(D0, grad, d0, tid);
		}
	}
}


template <int ROISize, int ROIPixelNum, int FitParaNum>
__device__ void poissonfGradient_GS2D(float ImageROI[][ROIPixelNum], float Ininf[][ThreadsPerBlock], float grad[][ThreadsPerBlock], float* WLE_Weight, int tid)
{
	float tIninf[FitParaNum];
	float tgradn;
	float tgradp;
	int cnt;

#pragma unroll
	for (cnt = 0; cnt < FitParaNum; cnt++)
	{
		tIninf[cnt] = Ininf[cnt][tid];
	}

	for (cnt = 0; cnt<FitParaNum; cnt++)
	{
		tIninf[cnt] = tIninf[cnt] - 0.002f;
		tgradn = MLEFit_TargetF_GS2D<ROISize, ROIPixelNum, FitParaNum>(ImageROI, tIninf, WLE_Weight, tid);

		tIninf[cnt] = tIninf[cnt] + 0.004f;
		tgradp = MLEFit_TargetF_GS2D<ROISize, ROIPixelNum, FitParaNum>(ImageROI, tIninf, WLE_Weight, tid);

		grad[cnt][tid] = (tgradp - tgradn)*250.0f; // /0.004f;
		tIninf[cnt] = tIninf[cnt] - 0.002f;
	}
}



template <int ROISize, int ROIPixelNum, int FitParaNum>
__device__ float MLEFit_TargetF_GS2D(float ImageROI[][ROIPixelNum], float Ininf[], float* WLE_Weight, int tid)
{
	int row, col;

	float rowpos, colpos;
	float TotalLoss = 0;
	float IVal;
	float pixval;
	float tdat;

	float(*pSubImg)[ROISize] = (float(*)[ROISize])&ImageROI[tid][0];

	float(*pWLE_Weight)[ROISize] = (float(*)[ROISize])WLE_Weight;


#pragma unroll
	for (row = 0; row < FitParaNum; row++)
	{
		if (Ininf[row]< 0.0f)Ininf[row] = 0.0f;
	}

	// slightly improve performance than explicitly define a variable
#define Peak_Scaled   Ininf[0] * AScalFactor
#define XPos_Scaled   Ininf[1]
#define YPos_Scaled   Ininf[2]
#define Sigm_Scaled   Ininf[3] * SScalFactor
#define Bakg_Scaled   Ininf[4] * BScalFactor


	for (row = 0; row<ROISize; row++)
	{
		for (col = 0; col<ROISize; col++)
		{
			// pf goal function
			pixval = pSubImg[row][col]; //ImageROI[tid][row*ROISize + col]; // coloffset+

			rowpos = row + 0.5f; // pixel center position
			colpos = col + 0.5f;

			// get model intensity
			tdat = -((colpos - XPos_Scaled)*(colpos - XPos_Scaled) + (rowpos - YPos_Scaled)*(rowpos - YPos_Scaled))*Sigm_Scaled;

			IVal = Peak_Scaled * __expf(tdat) + Bakg_Scaled;

			tdat = IVal - pixval - pixval*(__logf(IVal) - __logf(pixval));


#if(WLE_ENABLE == 1)
			// weighted MLE, improve localization precision
			TotalLoss = TotalLoss + pWLE_Weight[row][col] *tdat;
#else
			TotalLoss = TotalLoss + tdat;

#endif // WLE_ENABLE

		}
	}

	return TotalLoss;
}


template <int ROISize, int ROIPixelNum, int FitParaNum>
__device__ void PreEstimation_GS2D(float ImageROI[][ROIPixelNum], float Ininf[][ThreadsPerBlock], float WLE_SigmaX, int tid)
{
	float(*pSubImg)[ROISize] = (float(*)[ROISize])&ImageROI[tid][0];

	int CenterPos = (int)(ROISize / 2);

	// get bkg
	float Bkg = 0;
	float Mean1 = 0;
	float Mean2 = 0;
	float Mean3 = 0;
	float Mean4 = 0;

#pragma unroll
	for (int cnt = 0; cnt < ROISize; cnt++)
	{
		Mean1 += pSubImg[cnt][0];
		Mean2 += pSubImg[cnt][ROISize - 1];
		Mean3 += pSubImg[0][cnt];
		Mean4 += pSubImg[ROISize - 1][cnt];
	}
	Mean1 = min(Mean1, Mean2);
	Mean3 = min(Mean3, Mean4);

	Bkg = (Mean1 + Mean3) / 2 / ROISize;

	float Amp = pSubImg[CenterPos][CenterPos] - Bkg;


	// estimate sigma width

#if(WLE_ENABLE == 1)

	float SigmaX = WLE_SigmaX;

#else
	float SigmaX = ROISize / 2.2f / 2.35f;

#endif //WLE_ENABLE


	//
	float ParaSigma1 = (1.0f / (2.0f*SigmaX*SigmaX));

	// scaling
	Ininf[Fit2D_Peak][tid] = Amp*rAScalFactor; // /256.0f;	// Amplitude

	Ininf[Fit2D_XPos][tid] = (float)CenterPos + 0.5f; // sumDatx / sumdat - 0.5f;	// x
	Ininf[Fit2D_YPos][tid] = (float)CenterPos + 0.5f; // sumDaty / sumdat - 0.5f;	// y
	Ininf[Fit2D_Sigm][tid] = ParaSigma1 * rSScalFactor;	// sigma

	Ininf[Fit2D_Bakg][tid] = Bkg*rBScalFactor; // /128.0f;	// bg

}

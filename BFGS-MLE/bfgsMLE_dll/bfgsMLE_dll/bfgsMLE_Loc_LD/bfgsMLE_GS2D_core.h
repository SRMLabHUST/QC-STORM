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
#include "WLEParaEstimation_core.h"


template <int ROISize, int SharedImageWidth1, int FitParaNum, int IterateNum, int IterateNum_bs>
__device__ void MLELocalization_GS2D(float subimg[][SharedImageWidth1], float Ininf[][ThreadsPerBlock], float grad[][ThreadsPerBlock], float d0[][ThreadsPerBlock], float D0[], float* WLE_ParaArray, int tid);


template <int ROISize, int SharedImageWidth1, int FitParaNum>
__device__ void poissonfGradient_GS2D(float subimg[][SharedImageWidth1], float Ininf[][ThreadsPerBlock], float grad[][ThreadsPerBlock], float* WLE_ParaArray, int tid);


template <int ROISize, int SharedImageWidth1, int FitParaNum>
__device__ float MLEFit_TargetF_GS2D(float subimg[][SharedImageWidth1], float Ininf[], float* WLE_ParaArray, int tid);


template <int ROISize, int SharedImageWidth1, int FitParaNum>
__device__ void PreEstimation_GS2D(float subimg[][SharedImageWidth1], float Ininf[][ThreadsPerBlock], float WLE_SigmaX, int tid);


// algorithm core codes
template <int ROISize, int FitParaNum>
__global__ void bfgsMLELoc_Gauss2D(float *d_LocArry, unsigned short *d_SubRegion, float *d_WLEPara, float Offset, float kadc, float QE, int FluoNum)
{
	enum {
		ROIDataLen = (ROISize*ROISize),
		ROIWholeLen = (ROISize*(ROISize + 1)),

		LoadNum = (ROIDataLen / ThreadsPerBlock),
		LoadResidual = (ROIDataLen%ThreadsPerBlock),

		// avoid bank conflict, if no conflict, +1 can be avoid
		SharedImageWidth1 = (ROIWholeLen + 1),

		HalfRegionLen = (int)(ROISize / 2),
	};


#define XYPosVaryBound		1.5f

	float XYLBound = ROISize / 2.0f - XYPosVaryBound; // 1.0f
	float XYUBound = ROISize / 2.0f + XYPosVaryBound;


	//
	__shared__ float ImageRegion[ThreadsPerBlock][SharedImageWidth1];
	__shared__ float D0[ThreadsPerBlock][FitParaNum*FitParaNum];	// inv of matrix hessian 

	// avoid bank conflict
	__shared__ float Ininf[FitParaNum][ThreadsPerBlock]; // cur position
	__shared__ float grad[FitParaNum][ThreadsPerBlock];  // gradient
	__shared__ float d0[FitParaNum][ThreadsPerBlock];	// direction


	int gid = threadIdx.x + blockDim.x*blockIdx.x;
	int tid = threadIdx.x;

	//	if (gid >= FluoNum)return; // conflict with the __syncthreads below


	int BlockOffset = blockDim.x*blockIdx.x;
	int gMemPos;

	float(*pLocArry)[OutParaNumGS2D]; // for parameter array
	pLocArry = (float(*)[OutParaNumGS2D])d_LocArry;

	// read sub image into shared memory
	int cnt = 0;

	float(*pD0)[FitParaNum] = (float(*)[FitParaNum])&D0[tid][0];

	// WLE parameter array
	float(*pWLEPara)[WLE_ParaNumber] = (float(*)[WLE_ParaNumber])d_WLEPara;

	float WLE_ParaArray[WLE_ParaArrayLen];



	// load image region from global memory, to achieve Aligned and Coalesced Access
	int loadcnt = 0;
#pragma unroll
	for (cnt = 0; cnt < ThreadsPerBlock; cnt++)
	{
		gMemPos = (BlockOffset + cnt)*ROIWholeLen;
#pragma unroll
		for (loadcnt = 0; loadcnt < LoadNum; loadcnt++)
		{
			//to avoid negative value in logf in target function, that will be a disaster
			ImageRegion[cnt][loadcnt*ThreadsPerBlock + tid] = fmaxf((d_SubRegion[gMemPos + loadcnt*ThreadsPerBlock + tid] - Offset)*kadc, 1.0f);
		}
	}
	if (tid < LoadResidual)
	{
#pragma unroll
		for (cnt = 0; cnt < ThreadsPerBlock; cnt++)
		{
			gMemPos = (BlockOffset + cnt)*ROIWholeLen;
			ImageRegion[cnt][LoadNum*ThreadsPerBlock + tid] = fmaxf((d_SubRegion[gMemPos + LoadNum*ThreadsPerBlock + tid] - Offset)*kadc, 1.0f);
		}
	}
	if (tid < ROISize)
	{
#pragma unroll
		for (cnt = 0; cnt < ThreadsPerBlock; cnt++)
		{
			gMemPos = (BlockOffset + cnt)*ROIWholeLen;
			ImageRegion[cnt][ROIDataLen + tid] = d_SubRegion[gMemPos + ROIDataLen + tid];
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
	

	// why this function impact the time so much?
	__syncthreads();


	if (gid < FluoNum)
	{

		// estimated PSF width for mask of WLE, weighted MLE
		int MoleculeType = 0;

		float WLE_SigmaX = 0;


#if(WLE_ENABLE == 1)
		float SigmaL = pWLEPara[gid][WLE_Para_SigmaL];
		float SigmaR = pWLEPara[gid][WLE_Para_SigmaR];
		float SigmaU = pWLEPara[gid][WLE_Para_SigmaU];
		float SigmaD = pWLEPara[gid][WLE_Para_SigmaD];

		MoleculeType = pWLEPara[gid][WLE_Para_FluoType];


		SigmaL = min(SigmaL, SigmaR);
		SigmaU = min(SigmaU, SigmaD);

		WLE_SigmaX = min(SigmaL, SigmaU);

		
#endif // WLE_ENABLE


		// pre-estimation
		PreEstimation_GS2D<ROISize, SharedImageWidth1, FitParaNum>(ImageRegion, Ininf, WLE_SigmaX, tid);


		if (MoleculeType == MoleculeType_MLEFit)
		{
			WLE_SigmaX = WLE_SigmaX*1.2f;

			WLE_SigmaX = max(WLE_SigmaX, WLE_SigmaX / 2.0f / 2.35f);
		}
		else
		{
			WLE_SigmaX = WLE_SigmaX /1.2f;
		}

		WLE_ParaArray[WLE_Fit_SigmaX] = WLE_SigmaX;
		WLE_ParaArray[WLE_Fit_SigmaY] = WLE_SigmaX; // not necessary for 2d


		poissonfGradient_GS2D<ROISize, SharedImageWidth1, FitParaNum>(ImageRegion, Ininf, grad, WLE_ParaArray, tid);

#pragma unroll
		for (rcnt = 0; rcnt < FitParaNum; rcnt++)
		{
			d0[rcnt][tid] = -grad[rcnt][tid];
		}


		MLELocalization_GS2D<ROISize, SharedImageWidth1, FitParaNum_2D, IterateNum_2D, IterateNum_bs_2D>(ImageRegion, Ininf, grad, d0, &D0[tid][0], WLE_ParaArray, tid);

		// remove scaling factor
		float PeakPhoton = Ininf[Fit2D_Peak][tid] * AScalFactor; // PeakPhoton
		float Background = Ininf[Fit2D_Bakg][tid] * BScalFactor; // bkg

		float FittedXPos = Ininf[Fit2D_XPos][tid];
		float FittedYPos = Ininf[Fit2D_YPos][tid];

		float SimgaX = Ininf[Fit2D_Sigm][tid] * SScalFactor;

		SimgaX = sqrtf(0.5f / SimgaX); // pixel

									   // remove nan
		FittedXPos = isnan(FittedXPos) ? -1 : FittedXPos;
		FittedYPos = isnan(FittedYPos) ? -1 : FittedYPos;
		SimgaX = isnan(SimgaX) ? -1 : SimgaX;

		float CurSNR = PeakPhoton / sqrtf(PeakPhoton + Background); // e-

																	// convert e- to photon
		PeakPhoton = PeakPhoton / QE;
		Background = Background / QE;

		float TotalPhoton = PeakPhoton * 2 * 3.141593f*SimgaX*SimgaX;


		float XPos = FittedXPos - HalfRegionLen + ImageRegion[tid][ROIDataLen + 0]; // X
		float YPos = FittedYPos - HalfRegionLen + ImageRegion[tid][ROIDataLen + 1]; // Y

		int CurFrame = ImageRegion[tid][ROIDataLen + 3] * 65536 + ImageRegion[tid][ROIDataLen + 2];

		//	printf("dat:%.2f %.2f %.2f %.2f %.2f %d\n", PeakPhoton, XPos, YPos, SimgaX, Background, CurFrame);

#if(WLE_TEST)
		if (1) //
#else
		if ((MoleculeType <= 1) && (FittedXPos > XYLBound) && (FittedXPos < XYUBound) && (FittedYPos > XYLBound) && (FittedYPos < XYUBound)) //
#endif //WLE_TEST
		{

			pLocArry[gid][Pos_PPho] = PeakPhoton; // peak photon
			pLocArry[gid][Pos_XPos] = XPos; // may have 0.5 or 1 pixel offset compared with other software
			pLocArry[gid][Pos_YPos] = YPos; // may have 0.5 or 1 pixel offset compared with other software
			pLocArry[gid][Pos_ZPos] = 0.0f; // 
			pLocArry[gid][Pos_SigX] = SimgaX; // sigma x
			pLocArry[gid][Pos_SigY] = SimgaX; // sigma y
			pLocArry[gid][Pos_TPho] = TotalPhoton;   // total photon
			pLocArry[gid][Pos_Bakg] = Background; // background
			pLocArry[gid][Pos_PSNR] = CurSNR;        // peal snr
			pLocArry[gid][Pos_Frme] = CurFrame; // frame

		}
		else
		{

			// invalid point
#pragma unroll
			for (rcnt = 0; rcnt < OutParaNumGS2D; rcnt++)
			{
				pLocArry[gid][rcnt] = 0;

			}
		}
	}
}


template <int ROISize, int SharedImageWidth1, int FitParaNum, int IterateNum, int IterateNum_bs>
__device__ void MLELocalization_GS2D(float subimg[][SharedImageWidth1], float Ininf[][ThreadsPerBlock], float grad[][ThreadsPerBlock], float d0[][ThreadsPerBlock], float D0[], float* WLE_ParaArray, int tid)
{
	// adjust d0
	float td0[FitParaNum];
	float td0_total;

	float sk[FitParaNum]; // bfgs quasi-newton method
	float yk[FitParaNum];
	float tgrad[FitParaNum];

	int cnt = 0;
	int rcnt;

	int itcnt = 0; // iteration number

	// find work length, divide 2 method


	float scale;

	float xd[FitParaNum * 2];

	float ddat[2];
	float dpos[2];

	int xdsel = 0;

	for (itcnt = 0; itcnt<IterateNum; itcnt++)
	{
		td0_total = 0;

		// adjust d0
#pragma unroll
		for (rcnt = 0; rcnt < FitParaNum; rcnt++)
		{
			td0[rcnt] = abs(d0[rcnt][tid]);
			td0_total += td0[rcnt];
		}
		// reciprocal of mean
		td0_total = ((float)FitParaNum) / td0_total;


#pragma unroll
		for (rcnt = 0; rcnt < FitParaNum; rcnt++)
		{
			d0[rcnt][tid] = d0[rcnt][tid] * td0_total;

		}


		dpos[0] = 0.00001f; // scale factor left limit, should not equal to 0 and smaller
		dpos[1] = 1.0f; // scale factor right limit, should not lager than 2

		VectorAddMul1<FitParaNum>(&xd[0], Ininf, d0, 0.0001f, tid);
		VectorAddMul1<FitParaNum>(&xd[FitParaNum], Ininf, d0, 1.0f, tid);

		ddat[0] = MLEFit_TargetF_GS2D<ROISize, SharedImageWidth1, FitParaNum>(subimg, &xd[0], WLE_ParaArray, tid);
		ddat[1] = MLEFit_TargetF_GS2D<ROISize, SharedImageWidth1, FitParaNum>(subimg, &xd[FitParaNum], WLE_ParaArray, tid);
		// 
		for (cnt = 0; cnt < IterateNum_bs; cnt++)
		{
			// which part shrink
			xdsel = (ddat[0] < ddat[1]);


			dpos[xdsel] = (dpos[0] + dpos[1])*0.5f; //  /2.0f which one shift

			if (cnt < IterateNum_bs - 1)
			{
				VectorAddMul1<FitParaNum>(&xd[xdsel*FitParaNum], Ininf, d0, dpos[xdsel], tid);// xd=ininf+d0*scale
				ddat[xdsel] = MLEFit_TargetF_GS2D<ROISize, SharedImageWidth1, FitParaNum>(subimg, &xd[xdsel*FitParaNum], WLE_ParaArray, tid);
			}

		}
		scale = (dpos[0] + dpos[1])*0.5f; //  /2.0f calculated direction

		// sk 		= d0*base;
		// inInf 	= inInf+d0*base;
#pragma unroll
		for (rcnt = 0; rcnt < FitParaNum; rcnt++)
		{
			sk[rcnt] = d0[rcnt][tid] * scale;
			Ininf[rcnt][tid] = Ininf[rcnt][tid] + sk[rcnt];

		}

		IninfConstrain<FitParaNum>(Ininf, tid); // 


		if (itcnt < IterateNum - 1)
		{
			// tgrad = grad;
#pragma unroll
			for (rcnt = 0; rcnt < FitParaNum; rcnt++)
			{
				tgrad[rcnt] = grad[rcnt][tid];
			}

			poissonfGradient_GS2D<ROISize, SharedImageWidth1, FitParaNum>(subimg, Ininf, grad, WLE_ParaArray, tid);

#pragma unroll
			for (rcnt = 0; rcnt < FitParaNum; rcnt++)
			{
				yk[rcnt] = grad[rcnt][tid] - tgrad[rcnt];
			}

			ConstructD0<FitParaNum>(D0, sk, yk, tid);
			MatMulVector<FitParaNum>(D0, grad, d0, tid);
		}
	}
}


template <int ROISize, int SharedImageWidth1, int FitParaNum>
__device__ void poissonfGradient_GS2D(float subimg[][SharedImageWidth1], float Ininf[][ThreadsPerBlock], float grad[][ThreadsPerBlock], float* WLE_ParaArray, int tid)
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
		tgradn = MLEFit_TargetF_GS2D<ROISize, SharedImageWidth1, FitParaNum>(subimg, tIninf, WLE_ParaArray, tid);

		tIninf[cnt] = tIninf[cnt] + 0.004f;
		tgradp = MLEFit_TargetF_GS2D<ROISize, SharedImageWidth1, FitParaNum>(subimg, tIninf, WLE_ParaArray, tid);

		grad[cnt][tid] = (tgradp - tgradn)*250.0f; // /0.004f;
		tIninf[cnt] = tIninf[cnt] - 0.002f;
	}
}




template <int ROISize, int SharedImageWidth1, int FitParaNum>
__device__ float MLEFit_TargetF_GS2D(float subimg[][SharedImageWidth1], float Ininf[], float* WLE_ParaArray, int tid)
{
	int row, col;

	float rowpos, colpos;
	float TotalLoss = 0;
	float IVal;
	float pixval;
	float tdat;

	float(*pSubImg)[ROISize] = (float(*)[ROISize])&subimg[tid][0];

	// weighted MLE
	float WLE_Weight = 0;

	float WLE_SigmaX1 = WLE_ParaArray[WLE_Fit_SigmaX];

	WLE_SigmaX1 = 1.0f / (2.0f * WLE_SigmaX1 * WLE_SigmaX1);

	float ROICenter = ROISize / 2.0f; // use pre-estimated center pos don't help improve the localization precision



#pragma unroll
	for (row = 0; row < FitParaNum; row++)
	{
		if (Ininf[row]< 0.0f)Ininf[row] = 0.0f;
	}

	// slightly improve performance than explicitly define a variable
#define Peak_Scale   Ininf[0] * AScalFactor
#define XPos_Scale   Ininf[1]
#define YPos_Scale   Ininf[2]
#define Sigm_Scale   Ininf[3] * SScalFactor
#define Bakg_Scale   Ininf[4] * BScalFactor


	for (row = 0; row<ROISize; row++)
	{
		for (col = 0; col<ROISize; col++)
		{
			// pf goal function
			pixval = pSubImg[row][col]; //subimg[tid][row*ROISize + col]; // coloffset+

			rowpos = row + 0.5f; // pixel center position
			colpos = col + 0.5f;

			// get model intensity
			tdat = -((colpos - XPos_Scale)*(colpos - XPos_Scale) + (rowpos - YPos_Scale)*(rowpos - YPos_Scale))*Sigm_Scale;

			IVal = Peak_Scale * __expf(tdat) + Bakg_Scale;

			tdat = IVal - pixval - pixval*(__logf(IVal) - __logf(pixval));


#if(WLE_ENABLE == 1)
			// weighted MLE, improve localization precision
			WLE_Weight = __expf(-((colpos - ROICenter)*(colpos - ROICenter) + (rowpos - ROICenter)*(rowpos - ROICenter))*WLE_SigmaX1);

			TotalLoss = TotalLoss + WLE_Weight*tdat;
#else
			TotalLoss = TotalLoss + tdat;

#endif // WLE_ENABLE

		}
	}

	return TotalLoss;
}


template <int ROISize, int SharedImageWidth1, int FitParaNum>
__device__ void PreEstimation_GS2D(float subimg[][SharedImageWidth1], float Ininf[][ThreadsPerBlock], float WLE_SigmaX, int tid)
{
	float(*pSubImg)[ROISize] = (float(*)[ROISize])&subimg[tid][0];

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
	float SigmaX = ROISize / 2.0f / 2.35f;

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

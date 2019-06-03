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
#define FittedFluoNum_AS3D_2E			2

#define IterateNum_AS3D_2E				16  // total iteration number, fixed iteration number
#define IterateNum_AS3D_2E_bs			11 // bisection iteration to find best walk length


// number of fitting parameter
#define FitParaNum_AS3D_2E				11

#define FreeParaPerFluo					5


// fitting parameters ID
#define FitAS3D_2E_Peak0				0
#define FitAS3D_2E_XPos0				1
#define FitAS3D_2E_YPos0				2
#define FitAS3D_2E_SigX0				3
#define FitAS3D_2E_SigY0				4

#define FitAS3D_2E_Peak1				5
#define FitAS3D_2E_XPos1				6
#define FitAS3D_2E_YPos1				7
#define FitAS3D_2E_SigX1				8
#define FitAS3D_2E_SigY1				9

#define FitAS3D_2E_Bakg					10



// better to disable
#define WLE_ENABLE_MultiEmitterFit_3D				(0 && WLE_ENABLE)





template <int ROISize, int ROIPixelNum, int FitParaNum, int IterateNum, int IterateNum_bs>
__device__ void MLELocalization_AS3D_2E(float ImageROI[][ROIPixelNum], float Ininf[][ThreadsPerBlock], float grad[][ThreadsPerBlock], float d0[][ThreadsPerBlock], float D0[], float* WLE_Weight, int tid);


template <int ROISize, int ROIPixelNum, int FitParaNum>
__device__ void poissonfGradient_AS3D_2E(float ImageROI[][ROIPixelNum], float Ininf[][ThreadsPerBlock], float grad[][ThreadsPerBlock], float* WLE_Weight, int tid);


template <int ROISize, int ROIPixelNum, int FitParaNum>
__device__ float MLEFit_TargetF_AS3D_2E(float ImageROI[][ROIPixelNum], float Ininf[], float* WLE_Weight, int tid);


template <int ROISize, int ROIPixelNum, int FitParaNum>
__device__ void PreEstimation_AS3D_2E(float ImageROI[][ROIPixelNum], float Ininf[][ThreadsPerBlock], int tid);



template <int ROISize, int FitParaNum>
__global__ void bfgsMLELoc_AS3D_2E(float *d_LocArry, unsigned short *d_ImageROI, float *d_WLEPara, int MultiFitFluoNum_2E, FitPosInf_TypeDef* d_FitPosInf, CoreFittingPara *d_FitPara)
{
	
	enum {
		ROIPixelNum = (ROISize*ROISize),
		ROIWholeSize = (ROISize*(ROISize + 1)),


		ROISize_Half = (int)(ROISize / 2),

		LoadLoopNum = ((ROIPixelNum + ThreadsPerBlock - 1) / ThreadsPerBlock)

	};


#define XYPosVaryBound			2.5f

	float XYLBound = ROISize / 2.0f - XYPosVaryBound; // 1.5f
	float XYUBound = ROISize / 2.0f + XYPosVaryBound; // 1.5f
	
	float Simga_Max_Th = ROISize / 1.9f / 2.35f;
	float Simga_Min_Th = 0.8f;


	float QE = d_FitPara->QE;

	int * d_MultiFitFluoPos_2E = d_FitPosInf->d_MultiFitFluoPos_2E;
	int * d_MultiFit_AddedFluoNum = d_FitPosInf->d_MultiFit_AddedFluoNum;
	int * d_RejectedFluoNum = d_FitPosInf->d_RejectedFluoNum;


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

	const float ROICenter = ROISize / 2.0f;


	// load image region from global memory for multi emitter fitting
	for (int fcnt = 0; fcnt < ThreadsPerBlock; fcnt++)
	{
		int CurdId = gid_0 + fcnt;

		CurdId = min(CurdId, MultiFitFluoNum_2E - 1);

		int CurMapID = d_MultiFitFluoPos_2E[CurdId];
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


	//	if (gid >= FluoNum)return; // conflict with the __syncthreads below
	// why this function impact the time so much?
	__syncthreads();



	float Distance[FittedFluoNum_AS3D_2E];


	if (gid < MultiFitFluoNum_2E)
	{
		int CurFluoID = d_MultiFitFluoPos_2E[gid];


		int XOffset = d_ImageROI[CurFluoID*ROIWholeSize + ROIPixelNum + 0];
		int YOffset = d_ImageROI[CurFluoID*ROIWholeSize + ROIPixelNum + 1];
		int CurFrame = d_ImageROI[CurFluoID*ROIWholeSize + ROIPixelNum + 2] + 65536 * d_ImageROI[CurFluoID*ROIWholeSize + ROIPixelNum + 3];

		// fixed weighted MLE PSF width
		float WLE_SigmaX = ROISize / 1.9f / 2.35f;
		float WLE_SigmaY = ROISize / 1.9f / 2.35f;

		// nearest neighbor molecule in ROI detected
		float NearNeighborDistance = 1000;

#if(WLE_ENABLE_MultiEmitterFit_3D == 1)
		float SigmaL = pWLEPara[CurFluoID][WLE_Para_SigmaL];
		float SigmaR = pWLEPara[CurFluoID][WLE_Para_SigmaR];
		float SigmaU = pWLEPara[CurFluoID][WLE_Para_SigmaU];
		float SigmaD = pWLEPara[CurFluoID][WLE_Para_SigmaD];


		SigmaL = min(SigmaL, SigmaR);
		SigmaU = min(SigmaU, SigmaD);

		WLE_SigmaX = SigmaL * 1.2f;
		WLE_SigmaY = SigmaU * 1.2f;

		WLE_SigmaX = min(WLE_SigmaX, ROISize / 1.6f / 2.35f);
		WLE_SigmaY = min(WLE_SigmaY, ROISize / 1.6f / 2.35f);

		NearNeighborDistance = pWLEPara[CurFluoID][WLE_Para_NearDistance];

#endif // WLE_ENABLE_MultiEmitterFit_3D


		// pre-estimation
		PreEstimation_AS3D_2E<ROISize, ROIPixelNum, FitParaNum>(ImageROI, Ininf, tid);

#if(WLE_ENABLE_MultiEmitterFit_3D == 1)

		if (d_FitPara->WLEEn == 0)
		{
			WLE_SigmaX = ROISize;
			WLE_SigmaY = ROISize;
		}

		WLEWeightsCalc<ROISize>(WLE_Weight, WLE_SigmaX, WLE_SigmaY, MoleculeType_MultiFit);

#endif // WLE_ENABLE_MultiEmitterFit_3D


		MLELocalization_AS3D_2E<ROISize, ROIPixelNum, FitParaNum, IterateNum_AS3D_2E, IterateNum_AS3D_2E_bs>(ImageROI, Ininf, grad, d0, &D0[tid][0], WLE_Weight, tid);


		// remove scaling factor
		Ininf[FitAS3D_2E_Peak0][tid] = Ininf[FitAS3D_2E_Peak0][tid] * AScalFactor / QE;
		Ininf[FitAS3D_2E_Peak1][tid] = Ininf[FitAS3D_2E_Peak1][tid] * AScalFactor / QE;


		Ininf[FitAS3D_2E_SigX0][tid] = sqrtf(0.5f / (Ininf[FitAS3D_2E_SigX0][tid] * SScalFactor));
		Ininf[FitAS3D_2E_SigY0][tid] = sqrtf(0.5f / (Ininf[FitAS3D_2E_SigY0][tid] * SScalFactor));
		Ininf[FitAS3D_2E_SigX1][tid] = sqrtf(0.5f / (Ininf[FitAS3D_2E_SigX1][tid] * SScalFactor));
		Ininf[FitAS3D_2E_SigY1][tid] = sqrtf(0.5f / (Ininf[FitAS3D_2E_SigY1][tid] * SScalFactor));


		float Background = Ininf[FitAS3D_2E_Bakg][tid] * BScalFactor / QE;


		// remove nan
		int IsFittingValid = isnan(Background) ? 0 : 1;


		// molecule selection
		for (int cnt = 0; cnt <FittedFluoNum_AS3D_2E; cnt++)
		{
			float curPeak = Ininf[FreeParaPerFluo * cnt + 0][tid];
			float curx = Ininf[FreeParaPerFluo * cnt + 1][tid];
			float cury = Ininf[FreeParaPerFluo * cnt + 2][tid];

			if (abs(curx - ROICenter) >= ROISize_Half)curx = -100;
			if (abs(cury - ROICenter) >= ROISize_Half)cury = -100;


			Distance[cnt] = sqrtf((curx - ROICenter)*(curx - ROICenter) + (cury - ROICenter)*(cury - ROICenter));
			Distance[cnt] = curPeak*(ROISize - Distance[cnt]);
		}


		int MoleculeSel0 = 0;
		int MoleculeSel1 = 0;

		for (int cnt = 0; cnt < FittedFluoNum_AS3D_2E; cnt++)
		{
			if (Distance[cnt] > Distance[MoleculeSel0])
			{
				MoleculeSel0 = cnt;
			}
		}

		MoleculeSel1 = 1 - MoleculeSel0;


		// the first molecule
		float PeakPhoton0 = Ininf[FreeParaPerFluo * MoleculeSel0 + 0][tid];

		float XPos0_0 = Ininf[FreeParaPerFluo * MoleculeSel0 + 1][tid];
		float YPos0_0 = Ininf[FreeParaPerFluo * MoleculeSel0 + 2][tid];

		float SimgaX0 = Ininf[FreeParaPerFluo * MoleculeSel0 + 3][tid];
		float SimgaY0 = Ininf[FreeParaPerFluo * MoleculeSel0 + 4][tid];


		float TotalPhoton0 = PeakPhoton0 * 2 * 3.141593f * SimgaX0 * SimgaY0;

		float CurSNR0 = sqrtf(QE) * PeakPhoton0 / sqrtf(PeakPhoton0 + Background);


		float SigmaDif0 = SimgaX0 * SimgaX0 - SimgaY0 * SimgaY0; //

		float CurZPos0 = 0;

		if (SigmaDif0 >= 0)
		{
			CurZPos0 = d_FitPara->p4_XGY*powf(SigmaDif0, 4) + d_FitPara->p3_XGY*powf(SigmaDif0, 3) + d_FitPara->p2_XGY*powf(SigmaDif0, 2) + d_FitPara->p1_XGY*SigmaDif0 + d_FitPara->p0_XGY; // nm
		}
		else
		{
			CurZPos0 = d_FitPara->p4_XLY*powf(SigmaDif0, 4) + d_FitPara->p3_XLY*powf(SigmaDif0, 3) + d_FitPara->p2_XLY*powf(SigmaDif0, 2) + d_FitPara->p1_XLY*SigmaDif0 + d_FitPara->p0_XLY; // nm
		}

		CurZPos0 = CurZPos0 * d_FitPara->ZDepthCorrFactor;


		float MeanSigma0 = (SimgaX0 + SimgaY0) / 2;
		
		// the second molecule
		float PeakPhoton1 = Ininf[FreeParaPerFluo * MoleculeSel1 + 0][tid];

		float XPos1_0 = Ininf[FreeParaPerFluo * MoleculeSel1 + 1][tid];
		float YPos1_0 = Ininf[FreeParaPerFluo * MoleculeSel1 + 2][tid];

		float SimgaX1 = Ininf[FreeParaPerFluo * MoleculeSel1 + 3][tid];
		float SimgaY1 = Ininf[FreeParaPerFluo * MoleculeSel1 + 4][tid];


		float TotalPhoton1 = PeakPhoton1 * 2 * 3.141593f * SimgaX1 * SimgaY1;

		float CurSNR1 = sqrtf(QE) * PeakPhoton1 / sqrtf(PeakPhoton1 + Background);


		float SigmaDif1 = SimgaX1 * SimgaX1 - SimgaY1 * SimgaY1; //

		float CurZPos1 = 0;

		if (SigmaDif1 >= 0)
		{
			CurZPos1 = d_FitPara->p4_XGY*powf(SigmaDif1, 4) + d_FitPara->p3_XGY*powf(SigmaDif1, 3) + d_FitPara->p2_XGY*powf(SigmaDif1, 2) + d_FitPara->p1_XGY*SigmaDif1 + d_FitPara->p0_XGY; // nm
		}
		else
		{
			CurZPos1 = d_FitPara->p4_XLY*powf(SigmaDif1, 4) + d_FitPara->p3_XLY*powf(SigmaDif1, 3) + d_FitPara->p2_XLY*powf(SigmaDif1, 2) + d_FitPara->p1_XLY*SigmaDif1 + d_FitPara->p0_XLY; // nm
		}


		CurZPos1 = CurZPos1 * d_FitPara->ZDepthCorrFactor;

		float MeanSigma1 = (SimgaX1 + SimgaY1) / 2;

		float Distance = sqrtf((XPos0_0 - XPos1_0)*(XPos0_0 - XPos1_0) + (YPos0_0 - YPos1_0)*(YPos0_0 - YPos1_0));


		// write the first molecule nearest from the center
		if (IsFittingValid && (CurSNR0 > 5.5f) && (MeanSigma0 > 0.6f) && (MeanSigma0 > Simga_Min_Th) && (MeanSigma0 <= Simga_Max_Th))
		{
			// for both online and offline localization

			pLocArry[CurFluoID][Pos_PPho] = PeakPhoton0; // peak photon
			pLocArry[CurFluoID][Pos_XPos] = XPos0_0 - ROISize_Half + XOffset; // the unit is pixel
			pLocArry[CurFluoID][Pos_YPos] = YPos0_0 - ROISize_Half + YOffset; // the unit is pixel
			pLocArry[CurFluoID][Pos_ZPos] = CurZPos0; // z depth, the unit is nm
			pLocArry[CurFluoID][Pos_SigX] = SimgaX0; // sigma x, the unit is pixel
			pLocArry[CurFluoID][Pos_SigY] = SimgaY0; // sigma y, the unit is pixel
			pLocArry[CurFluoID][Pos_TPho] = TotalPhoton0;   // total photon
			pLocArry[CurFluoID][Pos_Bakg] = Background; // background
			pLocArry[CurFluoID][Pos_PSNR] = CurSNR0;        // peak snr
			pLocArry[CurFluoID][Pos_Frme] = CurFrame; // frame


#if(!SPEED_TEST)

			// write the second
			if ((Distance > 1.5f) && (CurSNR1 > 5.2f) && (NearNeighborDistance > ROICenter) && (MeanSigma1 > 0.6f) && (MeanSigma1 > Simga_Min_Th) && (MeanSigma1 <= Simga_Max_Th) && (PeakPhoton1 > PeakPhoton0 / 3.0f) && (XPos1_0 > 1) && (YPos1_0 > 1) && (XPos1_0 < ROISize - 1) && (YPos1_0 < ROISize - 1))
			{
				int addpos = atomicAdd(d_MultiFit_AddedFluoNum, 1);

				if (addpos < MaxPointNum)
				{
					pLocArry[addpos][Pos_PPho] = PeakPhoton1; // peak photon
					pLocArry[addpos][Pos_XPos] = XPos1_0 - ROISize_Half + XOffset; // the unit is pixel
					pLocArry[addpos][Pos_YPos] = YPos1_0 - ROISize_Half + YOffset; // the unit is pixel
					pLocArry[addpos][Pos_ZPos] = CurZPos1; // z depth, the unit is nm
					pLocArry[addpos][Pos_SigX] = SimgaX1; // sigma x, the unit is pixel
					pLocArry[addpos][Pos_SigY] = SimgaY1; // sigma y, the unit is pixel
					pLocArry[addpos][Pos_TPho] = TotalPhoton1;   // total photon
					pLocArry[addpos][Pos_Bakg] = Background; // background
					pLocArry[addpos][Pos_PSNR] = CurSNR1;        // peak snr
					pLocArry[addpos][Pos_Frme] = CurFrame; // frame

				}

			}


#endif // SPEED_TEST
		}
		else
		{

			atomicAdd(d_RejectedFluoNum, 1);


			// invalid point
#pragma unroll
			for (rcnt = 0; rcnt < OutParaNumGS2D; rcnt++)
			{
				pLocArry[CurFluoID][rcnt] = 0;

			}


		}
	}


}

template <int ROISize, int ROIPixelNum, int FitParaNum, int IterateNum, int IterateNum_bs>
__device__ void MLELocalization_AS3D_2E(float ImageROI[][ROIPixelNum], float Ininf[][ThreadsPerBlock], float grad[][ThreadsPerBlock], float d0[][ThreadsPerBlock], float D0[], float* WLE_Weight, int tid)
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
	poissonfGradient_AS3D_2E<ROISize, ROIPixelNum, FitParaNum>(ImageROI, Ininf, grad, WLE_Weight, tid);

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

		ddat[0] = MLEFit_TargetF_AS3D_2E<ROISize, ROIPixelNum, FitParaNum>(ImageROI, &xd[0], WLE_Weight, tid);
		ddat[1] = MLEFit_TargetF_AS3D_2E<ROISize, ROIPixelNum, FitParaNum>(ImageROI, &xd[FitParaNum], WLE_Weight, tid);
	
		// bisection method to find optimal walk length
		for (int cnt = 0; cnt < IterateNum_bs; cnt++)
		{
			// which part shrink
			int xdsel = (ddat[0] < ddat[1]);


			dpos[xdsel] = (dpos[0] + dpos[1])*0.5f; //  /2.0f which one shift

			if (cnt < IterateNum_bs - 1)
			{
				VectorAddMul1<FitParaNum>(&xd[xdsel*FitParaNum], Ininf, d0, dpos[xdsel], tid);// xd=ininf+d0*scale
				ddat[xdsel] = MLEFit_TargetF_AS3D_2E<ROISize, ROIPixelNum, FitParaNum>(ImageROI, &xd[xdsel*FitParaNum], WLE_Weight, tid);
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

			poissonfGradient_AS3D_2E<ROISize, ROIPixelNum, FitParaNum>(ImageROI, Ininf, grad, WLE_Weight, tid);

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
__device__ void poissonfGradient_AS3D_2E(float ImageROI[][ROIPixelNum], float Ininf[][ThreadsPerBlock], float grad[][ThreadsPerBlock], float* WLE_Weight, int tid)
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
		tgradn = MLEFit_TargetF_AS3D_2E<ROISize, ROIPixelNum, FitParaNum>(ImageROI, tIninf, WLE_Weight, tid);

		tIninf[cnt] = tIninf[cnt] + 0.004f;
		tgradp = MLEFit_TargetF_AS3D_2E<ROISize, ROIPixelNum, FitParaNum>(ImageROI, tIninf, WLE_Weight, tid);

		grad[cnt][tid] = (tgradp - tgradn)*250.0f; // /0.004f;
		tIninf[cnt] = tIninf[cnt] - 0.002f;
	}
}



template <int ROISize, int ROIPixelNum, int FitParaNum>
__device__ float MLEFit_TargetF_AS3D_2E(float ImageROI[][ROIPixelNum], float Ininf[], float* WLE_Weight, int tid)
{
	int row, col;

	float rowpos, colpos;
	float TotalLoss = 0;

	float IVal0;
	float IVal1;
	float IVal_all;

	float pixval;
	float Loss_t;

	float(*pSubImg)[ROISize] = (float(*)[ROISize])&ImageROI[tid][0];

	float(*pWLE_Weight)[ROISize] = (float(*)[ROISize])WLE_Weight;



#pragma unroll
	for (row = 0; row < FitParaNum; row++)
	{
		if (Ininf[row]< 0.01f)Ininf[row] = 0.01f;
	}

#define Peak0_Scaled	Ininf[0] * AScalFactor
#define XPos0_Scaled	Ininf[1]
#define YPos0_Scaled	Ininf[2]
#define SigX0_Scaled	Ininf[3] * SScalFactor
#define SigY0_Scaled	Ininf[4] * SScalFactor

#define Peak1_Scaled	Ininf[5] * AScalFactor
#define XPos1_Scaled	Ininf[6]
#define YPos1_Scaled	Ininf[7]
#define SigX1_Scaled	Ininf[8] * SScalFactor
#define SigY1_Scaled	Ininf[9] * SScalFactor

#define Bakg_Scaled		Ininf[10] * BScalFactor


	for (row = 0; row<ROISize; row++)
	{
		for (col = 0; col<ROISize; col++)
		{
			// pf goal function
			pixval = pSubImg[row][col]; //ImageROI[tid][row*ROISize + col]; // coloffset+

			rowpos = row + 0.5f; // pixel center position
			colpos = col + 0.5f;

			// get model intensity

			IVal0 = Peak0_Scaled * __expf(-((colpos - XPos0_Scaled)*(colpos - XPos0_Scaled)*SigX0_Scaled + (rowpos - YPos0_Scaled)*(rowpos - YPos0_Scaled)*SigY0_Scaled));
			IVal1 = Peak1_Scaled * __expf(-((colpos - XPos1_Scaled)*(colpos - XPos1_Scaled)*SigX1_Scaled + (rowpos - YPos1_Scaled)*(rowpos - YPos1_Scaled)*SigY1_Scaled));

				
			IVal_all = IVal0 + IVal1 + Bakg_Scaled;

			Loss_t = IVal_all - pixval - pixval*(__logf(IVal_all) - __logf(pixval));


#if(WLE_ENABLE_MultiEmitterFit_3D == 1)
			// weighted MLE, improve localization precision
			TotalLoss = TotalLoss + pWLE_Weight[row][col] * Loss_t;
#else
			TotalLoss = TotalLoss + Loss_t;

#endif // WLE_ENABLE_MultiEmitterFit_3D

		}
	}

	return TotalLoss;
}


template <int ROISize, int ROIPixelNum, int FitParaNum>
__device__ void PreEstimation_AS3D_2E(float ImageROI[][ROIPixelNum], float Ininf[][ThreadsPerBlock], int tid)
{
	float(*pSubImg)[ROISize] = (float(*)[ROISize])&ImageROI[tid][0];

	int ROISize_Half = ROISize / 2;
	float ROICenter = ROISize / 2.0f;


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

	float Amp = pSubImg[ROISize_Half][ROISize_Half] - Bkg;

	// estimate sigma width
	float SigmaX = ROISize / 2.2f / 2.35f;
	float SigmaY = ROISize / 2.2f / 2.35f;

	float SigmaX1 = rSScalFactor*(0.5f / (SigmaX*SigmaX));
	float SigmaY1 = rSScalFactor*(0.5f / (SigmaY*SigmaY));

	//////////////// estimate two emitter fitting parameter

	// intensity of 9 region and find the ones with max intensity
	float SumData[9];

	int SumROILen = ROISize / 3 + 1;
	int SumROILen_Half = SumROILen / 2;

	int StartPos[3];
	StartPos[0] = 0;
	StartPos[1] = ROISize_Half - SumROILen_Half;
	StartPos[2] = ROISize - SumROILen;

	float RegionPos[3];
	RegionPos[0] = ROICenter - SumROILen / 2.0f;
	RegionPos[1] = ROICenter;
	RegionPos[2] = ROICenter + SumROILen / 2.0f;

	for (int cnt = 0; cnt < 9; cnt++)
	{
		float SumDat_t = 0;

		int XBias = StartPos[cnt % 3];
		int YBias = StartPos[cnt / 3];

		for (int r = 0; r < SumROILen; r++)
		{
			for (int c = 0; c < SumROILen; c++)
			{
				SumDat_t += pSubImg[r + YBias][c + XBias];
			}
		}
		SumData[cnt] = SumDat_t;
	}

	// estimate the first
	SumData[9 / 2] = 0;

	int MaxPos0 = 0;
	for (int cnt = 0; cnt < 9; cnt++)
	{
		if (SumData[cnt] > SumData[MaxPos0])
		{
			MaxPos0 = cnt;
		}
	}
	float XPos1_Ini = RegionPos[MaxPos0 % 3];
	float YPos1_Ini = RegionPos[MaxPos0 / 3];


	// emitter 1 2, share the same  background
	Ininf[FitAS3D_2E_Peak0][tid] = 0.7f *Amp *rAScalFactor; // /256.0f	// Amplitude
	Ininf[FitAS3D_2E_XPos0][tid] = ROICenter; // x
	Ininf[FitAS3D_2E_YPos0][tid] = ROICenter; // y
	Ininf[FitAS3D_2E_SigX0][tid] = SigmaX1;	// sigma
	Ininf[FitAS3D_2E_SigY0][tid] = SigmaY1;	// sigma

	Ininf[FitAS3D_2E_Peak1][tid] = 0.7f *Amp *rAScalFactor; // /256.0f	// Amplitude
	Ininf[FitAS3D_2E_XPos1][tid] = XPos1_Ini; // x
	Ininf[FitAS3D_2E_YPos1][tid] = YPos1_Ini; // y
	Ininf[FitAS3D_2E_SigX1][tid] = SigmaX1;	// sigma
	Ininf[FitAS3D_2E_SigY1][tid] = SigmaY1;	// sigma

	Ininf[FitAS3D_2E_Bakg][tid] = 0.8f * Bkg *rBScalFactor; // /32.0f	// bg

}


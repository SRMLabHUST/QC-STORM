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


// fitting parameters for astigmatism 3d  single emitter
#define FitParaNum_AS3D			6

#define IterateNum_AS3D			11  // total iteration number, fixed iteration number
#define IterateNum_AS3D_bs		11 // bisection iteration to find best walk length

// fitting parameters
#define FitAS3D_Peak				0
#define FitAS3D_XPos				1
#define FitAS3D_YPos				2
#define FitAS3D_SigX				3
#define FitAS3D_SigY				4
#define FitAS3D_Bakg				5


template <int ROISize, int SharedImageWidth1, int FitParaNum, int IterateNum, int IterateNum_bs>
__device__ void MLELocalization_AS3D(float subimg[][SharedImageWidth1], float Ininf[][ThreadsPerBlock], float grad[][ThreadsPerBlock], float d0[][ThreadsPerBlock], float D0[], float* WLE_Weight, int tid);

template <int ROISize, int SharedImageWidth1, int FitParaNum>
__device__ void poissonfGradient_AS3D(float subimg[][SharedImageWidth1], float Ininf[][ThreadsPerBlock], float grad[][ThreadsPerBlock], float* WLE_Weight, int tid);


template <int ROISize, int SharedImageWidth1, int FitParaNum>
__device__ float MLEFit_TargetF_AS3D(float subimg[][SharedImageWidth1], float Ininf[], float* WLE_Weight, int tid);



template <int ROISize, int SharedImageWidth1, int FitParaNum>
__device__ void PreEstimation_AS3D(float subimg[][SharedImageWidth1], float Ininf[][ThreadsPerBlock], float WLE_SigmaX, float WLE_SigmaY, int tid);


// algorithm core codes
template <int ROISize, int FitParaNum>
__global__ void bfgsMLELoc_AS3D(float *d_LocArry, unsigned short *d_SubRegion, float *d_WLEPara, float Offset, float kadc, float QE, float ZDepthCorrFactor, float p4, float p3, float p2, float p1, float p0, int FluoNum)
{
	enum {
		ROIPixelNum = (ROISize*ROISize),
		ROIWholeSize = (ROISize*(ROISize + 1)),

		// avoid bank conflict, if no conflict, +1 can be avoid
		SharedImageWidth1 = (ROIWholeSize + 1),

		HalfRegionLen = (int)(ROISize / 2),

		LoadLoopNum = ((ROIWholeSize + ThreadsPerBlock - 1) / ThreadsPerBlock)

	};


#define XYPosVaryBound			1.5f

	float XYLBound = ROISize / 2.0f - XYPosVaryBound; // 1.5f
	float XYUBound = ROISize / 2.0f + XYPosVaryBound; // 1.5f

	float SigmaLBound = 0.8f;
	float SigmaUBound = ROISize / 2.0f; // 2.35f


	//
	__shared__ float ImageRegion[ThreadsPerBlock][SharedImageWidth1];
	__shared__ float D0[ThreadsPerBlock][FitParaNum*FitParaNum];	// inv of matrix hessian 

	// avoid bank conflict
	__shared__ float Ininf[FitParaNum][ThreadsPerBlock]; // cur position
	__shared__ float grad[FitParaNum][ThreadsPerBlock];  // gradient
	__shared__ float d0[FitParaNum][ThreadsPerBlock];	// direction


	int gid = threadIdx.x + blockDim.x*blockIdx.x;
	int gid_0 = blockDim.x*blockIdx.x;

	int tid = threadIdx.x;


	//	if (gid >= FluoNum)return; // conflict with the __syncthreads below
	// it's not necessary since d_LocArry size is larger than total fluo number



	float(*pLocArry)[OutParaNumGS2D]; // for parameter array
	pLocArry = (float(*)[OutParaNumGS2D])d_LocArry;


	float(*pD0)[FitParaNum] = (float(*)[FitParaNum])&D0[tid][0];

	// WLE parameter array
	float(*pWLEPara)[WLE_ParaNumber] = (float(*)[WLE_ParaNumber])d_WLEPara;

	float WLE_Weight[ROIPixelNum];



	// load image region from global memory
	for (int fcnt = 0; fcnt < ThreadsPerBlock; fcnt++)
	{
		int CurFluoId = gid_0 + fcnt;
		int AddrOffset = CurFluoId*ROIWholeSize;

#pragma unroll
		for (int cnt = 0; cnt < LoadLoopNum; cnt++)
		{
			int CurLoadId = cnt*ThreadsPerBlock + tid;

			if (CurLoadId < ROIPixelNum)
			{
				ImageRegion[fcnt][CurLoadId] = fmaxf((d_SubRegion[AddrOffset + CurLoadId] - Offset)*kadc, 1.0f);
			}
			else
			{
				ImageRegion[fcnt][CurLoadId] = d_SubRegion[AddrOffset + CurLoadId];
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


	// why this function impact the time so much?
	__syncthreads();


	if (gid < FluoNum)
	{

		// estimated PSF width for mask of WLE, weighted MLE
		int MoleculeType = 0;

		float WLE_SigmaX = 0;
		float WLE_SigmaY = 0;


#if(WLE_ENABLE == 1)
		float SigmaL = pWLEPara[gid][WLE_Para_SigmaL];
		float SigmaR = pWLEPara[gid][WLE_Para_SigmaR];
		float SigmaU = pWLEPara[gid][WLE_Para_SigmaU];
		float SigmaD = pWLEPara[gid][WLE_Para_SigmaD];

		MoleculeType = pWLEPara[gid][WLE_Para_FluoType];


		SigmaL = min(SigmaL, SigmaR);
		SigmaU = min(SigmaU, SigmaD);

		WLE_SigmaX = SigmaL;
		WLE_SigmaY = SigmaU;

#endif // WLE_ENABLE


		// pre-estimation
		PreEstimation_AS3D<ROISize, SharedImageWidth1, FitParaNum>(ImageRegion, Ininf, WLE_SigmaX, WLE_SigmaY, tid);


		if (MoleculeType == MoleculeType_MLEFit)
		{
			// larger Width for uncontaminated molecules
			WLE_SigmaX = WLE_SigmaX*1.2f;
			WLE_SigmaY = WLE_SigmaY*1.2f;

			WLE_SigmaX = max(WLE_SigmaX, ROISize / 1.5f / 2.35f);
			WLE_SigmaY = max(WLE_SigmaY, ROISize / 1.5f / 2.35f);
		}
		else
		{
			// smaller Width for contaminated molecules
			WLE_SigmaX = WLE_SigmaX / 1.2f;
			WLE_SigmaY = WLE_SigmaY / 1.2f;

			WLE_SigmaX = min(WLE_SigmaX, ROISize / 2.0f / 2.35f / 1.2f);
			WLE_SigmaY = min(WLE_SigmaY, ROISize / 2.0f / 2.35f / 1.2f);

		}


		WLEWeightsCalc<ROISize>(WLE_Weight, WLE_SigmaX, WLE_SigmaY);


		poissonfGradient_AS3D<ROISize, SharedImageWidth1, FitParaNum>(ImageRegion, Ininf, grad, WLE_Weight, tid);

#pragma unroll
		for (rcnt = 0; rcnt < FitParaNum; rcnt++)
		{
			d0[rcnt][tid] = -grad[rcnt][tid];

		}


		MLELocalization_AS3D<ROISize, SharedImageWidth1, FitParaNum, IterateNum_AS3D, IterateNum_AS3D_bs>(ImageRegion, Ininf, grad, d0, &D0[tid][0], WLE_Weight, tid);

		// remove scaling factor
		float PeakPhoton = Ininf[FitAS3D_Peak][tid] * AScalFactor;
		float Background = Ininf[FitAS3D_Bakg][tid] * BScalFactor;

		float FittedXPos = Ininf[FitAS3D_XPos][tid];
		float FittedYPos = Ininf[FitAS3D_YPos][tid];

		float SimgaX = Ininf[FitAS3D_SigX][tid] * SScalFactor;
		float SimgaY = Ininf[FitAS3D_SigY][tid] * SScalFactor;

		SimgaX = sqrtf(0.5f / SimgaX); // pixel
		SimgaY = sqrtf(0.5f / SimgaY); // pixel

									   // remove nan
		FittedXPos = isnan(FittedXPos) ? -1 : FittedXPos;
		FittedYPos = isnan(FittedYPos) ? -1 : FittedYPos;
		SimgaX = isnan(SimgaX) ? -1 : SimgaX;
		SimgaY = isnan(SimgaY) ? -1 : SimgaY;

		//
		float CurSNR = PeakPhoton / sqrtf(PeakPhoton + Background); // e-

																	// convert e- to photon
		PeakPhoton = PeakPhoton / QE;
		Background = Background / QE;

		float TotalPhoton = PeakPhoton * 2 * 3.141593f*SimgaX*SimgaY;


		float SigmaDif = SimgaX*SimgaX - SimgaY*SimgaY; //



		float CurZPos = p4*powf(SigmaDif, 4) + p3*powf(SigmaDif, 3) + p2*powf(SigmaDif, 2) + p1*SigmaDif + p0; // nm


		float XPos = FittedXPos - HalfRegionLen + ImageRegion[tid][ROIPixelNum + 0]; // X
		float YPos = FittedYPos - HalfRegionLen + ImageRegion[tid][ROIPixelNum + 1]; // Y

		int CurFrame = ImageRegion[tid][ROIPixelNum + 3] * 65536 + ImageRegion[tid][ROIPixelNum + 2];

#if(WLE_TEST)
		if (1) //
#else
		if ((MoleculeType <= 1) && (FittedXPos > XYLBound) && (FittedXPos < XYUBound) && (FittedYPos > XYLBound) && (FittedYPos < XYUBound) && (SimgaX > SigmaLBound) && (SimgaX < SigmaUBound) && (SimgaY > SigmaLBound) && (SimgaY < SigmaUBound) && (abs(CurZPos) < 1000.0f))
#endif //WLE_TEST
		{
			// for both online and offline localization

			pLocArry[gid][Pos_PPho] = PeakPhoton; // peak photon
			pLocArry[gid][Pos_XPos] = XPos; // the unit is pixel
			pLocArry[gid][Pos_YPos] = YPos; // the unit is pixel
			pLocArry[gid][Pos_ZPos] = CurZPos*ZDepthCorrFactor; // z depth, the unit is nm
			pLocArry[gid][Pos_SigX] = SimgaX; // sigma x, the unit is pixel
			pLocArry[gid][Pos_SigY] = SimgaY; // sigma y, the unit is pixel
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
__device__ void MLELocalization_AS3D(float subimg[][SharedImageWidth1], float Ininf[][ThreadsPerBlock], float grad[][ThreadsPerBlock], float d0[][ThreadsPerBlock], float D0[], float* WLE_Weight, int tid)
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

		ddat[0] = MLEFit_TargetF_AS3D<ROISize, SharedImageWidth1, FitParaNum>(subimg, &xd[0], WLE_Weight, tid);
		ddat[1] = MLEFit_TargetF_AS3D<ROISize, SharedImageWidth1, FitParaNum>(subimg, &xd[FitParaNum], WLE_Weight, tid);
		// 
		for (cnt = 0; cnt < IterateNum_bs; cnt++)
		{
			// which part shrink
			xdsel = (ddat[0] < ddat[1]);


			dpos[xdsel] = (dpos[0] + dpos[1])*0.5f; //  /2.0f which one shift

			if (cnt < IterateNum_bs - 1)
			{
				VectorAddMul1<FitParaNum>(&xd[xdsel*FitParaNum], Ininf, d0, dpos[xdsel], tid);// xd=ininf+d0*scale
				ddat[xdsel] = MLEFit_TargetF_AS3D<ROISize, SharedImageWidth1, FitParaNum>(subimg, &xd[xdsel*FitParaNum], WLE_Weight, tid);
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

			poissonfGradient_AS3D<ROISize, SharedImageWidth1, FitParaNum>(subimg, Ininf, grad, WLE_Weight, tid);

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
__device__ void poissonfGradient_AS3D(float subimg[][SharedImageWidth1], float Ininf[][ThreadsPerBlock], float grad[][ThreadsPerBlock], float* WLE_Weight, int tid)
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
		tgradn = MLEFit_TargetF_AS3D<ROISize, SharedImageWidth1, FitParaNum>(subimg, tIninf, WLE_Weight, tid);

		tIninf[cnt] = tIninf[cnt] + 0.004f;
		tgradp = MLEFit_TargetF_AS3D<ROISize, SharedImageWidth1, FitParaNum>(subimg, tIninf, WLE_Weight, tid);

		grad[cnt][tid] = (tgradp - tgradn)*250.0f; // /0.004f;
		tIninf[cnt] = tIninf[cnt] - 0.002f;
	}
}



template <int ROISize, int SharedImageWidth1, int FitParaNum>
__device__ float MLEFit_TargetF_AS3D(float subimg[][SharedImageWidth1], float Ininf[], float* WLE_Weight, int tid)
{
	int row, col;

	float rowpos, colpos;
	float TotalLoss = 0;
	float IVal;
	float pixval;
	float tdat;

	float(*pSubImg)[ROISize] = (float(*)[ROISize])&subimg[tid][0];

	float(*pWLE_Weight)[ROISize] = (float(*)[ROISize])WLE_Weight;



#pragma unroll
	for (row = 0; row < FitParaNum; row++)
	{
		if (Ininf[row]< 0.01f)Ininf[row] = 0.01f;
	}

#define Peak_Scaled	Ininf[0] * AScalFactor
#define XPos_Scaled	Ininf[1]
#define YPos_Scaled	Ininf[2]
#define SigX_Scaled	Ininf[3] * SScalFactor
#define SigY_Scaled	Ininf[4] * SScalFactor
#define Bakg_Scaled	Ininf[5] * BScalFactor
	
		
	for (row = 0; row<ROISize; row++)
	{
		for (col = 0; col<ROISize; col++)
		{
			// pf goal function
			pixval = pSubImg[row][col]; //subimg[tid][row*ROISize + col]; // coloffset+

			rowpos = row + 0.5f; // pixel center position
			colpos = col + 0.5f;

			// get model intensity
			tdat = -((colpos - XPos_Scaled)*(colpos - XPos_Scaled)*SigX_Scaled + (rowpos - YPos_Scaled)*(rowpos - YPos_Scaled)*SigY_Scaled);

			IVal = Peak_Scaled * __expf(tdat) + Bakg_Scaled;

			tdat = IVal - pixval - pixval*(__logf(IVal) - __logf(pixval)); //to avoid negative value in logf function, that will be a disaster

#if(WLE_ENABLE == 1)
			// weighted MLE, improve localization precision
			TotalLoss = TotalLoss + pWLE_Weight[row][col] * tdat;
#else
			TotalLoss = TotalLoss + tdat;

#endif // WLE_ENABLE

		}
	}

	return TotalLoss;
}


template <int ROISize, int SharedImageWidth1, int FitParaNum>
__device__ void PreEstimation_AS3D(float subimg[][SharedImageWidth1], float Ininf[][ThreadsPerBlock], float WLE_SigmaX, float WLE_SigmaY, int tid)
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
	float SigmaY = WLE_SigmaY;

#else
	float SigmaX = ROISize / 2.0f / 2.35f;
	float SigmaY = ROISize / 2.0f / 2.35f;

#endif //WLE_ENABLE


	//
	float SigmaX1 = rSScalFactor*(0.5f / (SigmaX*SigmaX));
	float SigmaY1 = rSScalFactor*(0.5f / (SigmaY*SigmaY));


	// scaling
	Ininf[FitAS3D_Peak][tid] = Amp *rAScalFactor; // /256.0f	// Amplitude

	Ininf[FitAS3D_XPos][tid] = (float)CenterPos + 0.5f; // x
	Ininf[FitAS3D_YPos][tid] = (float)CenterPos + 0.5f; // y
	Ininf[FitAS3D_SigX][tid] = SigmaX1;	// sigma
	Ininf[FitAS3D_SigY][tid] = SigmaY1;	// sigma

	Ininf[FitAS3D_Bakg][tid] = Bkg *rBScalFactor; // /32.0f	// bg

}


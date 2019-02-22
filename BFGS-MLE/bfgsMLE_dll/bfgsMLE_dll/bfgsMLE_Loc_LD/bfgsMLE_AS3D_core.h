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


template <int RegionLen, int SharedImageWidth1, int FitParaNum, int IterateNum, int IterateNum_bs>
__device__ void MLELocalization_AS3D(float subimg[][SharedImageWidth1], float Ininf[][ThreadsPerBlock], float grad[][ThreadsPerBlock], float d0[][ThreadsPerBlock], float D0[], float SigmaX_Est, float SigmaY_Est, int tid);

template <int RegionLen, int SharedImageWidth1, int FitParaNum>
__device__ void poissonfGradient_AS3D(float subimg[][SharedImageWidth1], float Ininf[][ThreadsPerBlock], float grad[][ThreadsPerBlock], float SigmaX_Est, float SigmaY_Est, int tid);


template <int RegionLen, int SharedImageWidth1, int FitParaNum>
__device__ float MLEFit_TargetF_AS3D(float subimg[][SharedImageWidth1], float Ininf[], float SigmaX_Est, float SigmaY_Est, int tid);



template <int RegionLen, int SharedImageWidth1, int FitParaNum>
__device__ void PreEstimation_AS3D(float subimg[][SharedImageWidth1], float Ininf[][ThreadsPerBlock], float *SigmaX_Est, float *SigmaY_Est, int tid);


// algorithm core codes
template <int RegionLen, int FitParaNum>
__global__ void bfgsMLELoc_AS3D(unsigned short *d_SubRegion, float *d_LocArry, float Offset, float kadc, float QE, float ZDepthCorrFactor, float p4, float p3, float p2, float p1, float p0, int FluoNum)
{
	enum {
		ROIDataLen = (RegionLen*RegionLen),
		ROIWholeLen = (RegionLen*(RegionLen + 1)),

		LoadNum = (ROIDataLen / ThreadsPerBlock),
		LoadResidual = (ROIDataLen%ThreadsPerBlock),

		// avoid bank conflict, if no conflict, +1 can be avoid
		SharedImageWidth1 = (ROIWholeLen + 1),

		HalfRegionLen = (int)(RegionLen / 2),

		D0_Size = FitParaNum*FitParaNum,
		
	};


	// better
	// #define XYPosVaryBound		1.5f
// test WLE
#define XYPosVaryBound			3.0f


	float XYLBound = RegionLen / 2.0f - XYPosVaryBound; // 1.5f
	float XYUBound = RegionLen / 2.0f + XYPosVaryBound; // 1.5f

	float SigmaLBound = 0.8f;
	float SigmaUBound = RegionLen / 2.0f; // 2.35f



	__shared__ float ImageRegion[ThreadsPerBlock][SharedImageWidth1];
	__shared__ float D0[ThreadsPerBlock][D0_Size];	// inv of matrix hessian 

	// avoid bank conflict
	__shared__ float Ininf[FitParaNum][ThreadsPerBlock]; // cur position
	__shared__ float grad[FitParaNum][ThreadsPerBlock];  // gradient
	__shared__ float d0[FitParaNum][ThreadsPerBlock];	// direction


	int gid = threadIdx.x + blockDim.x*blockIdx.x;
	int tid = threadIdx.x;

	//	if (gid >= FluoNum)return; // conflict with the __syncthreads below


	int BlockOffset = blockDim.x*blockIdx.x;
	int gMemPos;
	int gMemPosOffset = BlockOffset*ROIWholeLen;

	float(*pLocArry)[OutParaNumGS2D]; // for parameter array
	pLocArry = (float(*)[OutParaNumGS2D])d_LocArry;

	// read sub image into shared memory
	int cnt = 0;

	float(*pD0)[FitParaNum] = (float(*)[FitParaNum])&D0[tid][0];



	// load image region from global memory, there two method to load
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
	if (tid < RegionLen)
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

	// estimated PSF width for mask of WLE, weighted MLE
	float SigmaX_Est;
	float SigmaY_Est;

	PreEstimation_AS3D<RegionLen, SharedImageWidth1, FitParaNum>(ImageRegion, Ininf, &SigmaX_Est, &SigmaY_Est, tid);

	poissonfGradient_AS3D<RegionLen, SharedImageWidth1, FitParaNum>(ImageRegion, Ininf, grad, SigmaX_Est, SigmaY_Est, tid);

#pragma unroll
	for (rcnt = 0; rcnt < FitParaNum; rcnt++)
	{
		d0[rcnt][tid] = -grad[rcnt][tid];

	}

	
	MLELocalization_AS3D<RegionLen, SharedImageWidth1, FitParaNum, IterateNum_AS3D, IterateNum_bs_AS3D>(ImageRegion, Ininf, grad, d0, &D0[tid][0], SigmaX_Est, SigmaY_Est, tid);


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
	

	float XPos = FittedXPos - HalfRegionLen + ImageRegion[tid][ROIDataLen + 0]; // X
	float YPos = FittedYPos - HalfRegionLen + ImageRegion[tid][ROIDataLen + 1]; // Y
	
	int CurFrame = ImageRegion[tid][ROIDataLen + 3] * 65536 + ImageRegion[tid][ROIDataLen + 2];

	if ((FittedXPos > XYLBound) && (FittedXPos < XYUBound) && (FittedYPos > XYLBound) && (FittedYPos < XYUBound) && (SimgaX > SigmaLBound) && (SimgaX < SigmaUBound) && (SimgaY > SigmaLBound) && (SimgaY < SigmaUBound) && (abs(CurZPos) < 1000.0f))
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

template <int RegionLen, int SharedImageWidth1, int FitParaNum, int IterateNum, int IterateNum_bs>
__device__ void MLELocalization_AS3D(float subimg[][SharedImageWidth1], float Ininf[][ThreadsPerBlock], float grad[][ThreadsPerBlock], float d0[][ThreadsPerBlock], float D0[], float SigmaX_Est, float SigmaY_Est, int tid)
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

		ddat[0] = MLEFit_TargetF_AS3D<RegionLen, SharedImageWidth1, FitParaNum>(subimg, &xd[0], SigmaX_Est, SigmaY_Est, tid);
		ddat[1] = MLEFit_TargetF_AS3D<RegionLen, SharedImageWidth1, FitParaNum>(subimg, &xd[FitParaNum], SigmaX_Est, SigmaY_Est, tid);
		// 
		for (cnt = 0; cnt < IterateNum_bs; cnt++)
		{
			// which part shrink
			xdsel = (ddat[0] < ddat[1]);


			dpos[xdsel] = (dpos[0] + dpos[1])*0.5f; //  /2.0f which one shift

			if (cnt < IterateNum_bs - 1)
			{
				VectorAddMul1<FitParaNum>(&xd[xdsel*FitParaNum], Ininf, d0, dpos[xdsel], tid);// xd=ininf+d0*scale
				ddat[xdsel] = MLEFit_TargetF_AS3D<RegionLen, SharedImageWidth1, FitParaNum>(subimg, &xd[xdsel*FitParaNum], SigmaX_Est, SigmaY_Est, tid);
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

			poissonfGradient_AS3D<RegionLen, SharedImageWidth1, FitParaNum>(subimg, Ininf, grad, SigmaX_Est, SigmaY_Est, tid);

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


template <int RegionLen, int SharedImageWidth1, int FitParaNum>
__device__ void poissonfGradient_AS3D(float subimg[][SharedImageWidth1], float Ininf[][ThreadsPerBlock], float grad[][ThreadsPerBlock], float SigmaX_Est, float SigmaY_Est, int tid)
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
		tgradn = MLEFit_TargetF_AS3D<RegionLen, SharedImageWidth1, FitParaNum>(subimg, tIninf, SigmaX_Est, SigmaY_Est, tid);

		tIninf[cnt] = tIninf[cnt] + 0.004f;
		tgradp = MLEFit_TargetF_AS3D<RegionLen, SharedImageWidth1, FitParaNum>(subimg, tIninf, SigmaX_Est, SigmaY_Est, tid);

		grad[cnt][tid] = (tgradp - tgradn)*250.0f; // /0.004f;
		tIninf[cnt] = tIninf[cnt] - 0.002f;
	}
}



template <int RegionLen, int SharedImageWidth1, int FitParaNum>
__device__ void PreEstimation_AS3D(float subimg[][SharedImageWidth1], float Ininf[][ThreadsPerBlock], float *SigmaX_Est, float *SigmaY_Est, int tid)
{
	float(*pSubImg)[RegionLen] = (float(*)[RegionLen])&subimg[tid][0];

	int CenterPos = (int)(RegionLen / 2);


	// get bkg
	float Bkg = 0;
	float Mean1 = 0;
	float Mean2 = 0;
	float Mean3 = 0;
	float Mean4 = 0;

#pragma unroll
	for (int cnt = 0; cnt < RegionLen; cnt++)
	{
		Mean1 += pSubImg[cnt][0];
		Mean2 += pSubImg[cnt][RegionLen-1];
		Mean3 += pSubImg[0][cnt];
		Mean4 += pSubImg[RegionLen - 1][cnt];
	}
	Mean1 = min(Mean1, Mean2);
	Mean3 = min(Mean3, Mean4);

	Bkg = (Mean1 + Mean3) / 2 / RegionLen;

	float Amp = pSubImg[CenterPos][CenterPos] - Bkg;

	// get mean col and row
	float MeanX[RegionLen];
	float MeanY[RegionLen];

#pragma unroll
	for (int cnt = 0; cnt < RegionLen; cnt++)
	{
		MeanX[cnt] = 0;
		MeanY[cnt] = 0;
	}

	// average 
#pragma unroll
	for (int row = 1; row < RegionLen - 1; row++)
	{
#pragma unroll
		for (int col = 0; col < RegionLen; col++)
		{
			MeanX[col] = MeanX[col] + pSubImg[row][col];
		}
	}
#pragma unroll
	for (int row = 0; row < RegionLen; row++)
	{
#pragma unroll
		for (int col = 1; col < RegionLen - 1; col++)
		{
			MeanY[row] = MeanY[row] + pSubImg[row][col];
		}
	}


#pragma unroll
	for (int cnt = 0; cnt < RegionLen; cnt++)
	{
		MeanX[cnt] = MeanX[cnt] / (RegionLen - 2);
		MeanY[cnt] = MeanY[cnt] / (RegionLen - 2);
	}

	// remove background for precise estimation
	RemoveMeanDatBkg<RegionLen>(MeanX);
	RemoveMeanDatBkg<RegionLen>(MeanY);

	// estimate center posotion
	float XPos_rough = GetMeanDatCenterPos_Rough<RegionLen>(MeanX);
	float YPos_rough = GetMeanDatCenterPos_Rough<RegionLen>(MeanY);

	// estimate sigma width
	float SigmaX = GetPSFSigmaEstimation<RegionLen>(MeanX, XPos_rough);
	float SigmaY = GetPSFSigmaEstimation<RegionLen>(MeanY, YPos_rough);


	*SigmaX_Est = SigmaX;
	*SigmaY_Est = SigmaY;

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

template <int RegionLen, int SharedImageWidth1, int FitParaNum>
__device__ float MLEFit_TargetF_AS3D(float subimg[][SharedImageWidth1], float Ininf[], float SigmaX_Est, float SigmaY_Est, int tid)
{
	int row, col;

	float rowpos, colpos;
	float TotalLoss = 0;
	float IVal;
	float pixval;
	float tdat;

	float(*pSubImg)[RegionLen] = (float(*)[RegionLen])&subimg[tid][0];

	// weighted MLE
	float WLE_Weight = 0;
	float WLE_SigmaX1 = 1.0f / 2.0f / SigmaX_Est / SigmaX_Est;
	float WLE_SigmaY1 = 1.0f / 2.0f / SigmaY_Est / SigmaY_Est;
	float ROICenter = RegionLen / 2.0f;


#pragma unroll
	for (row = 0; row < FitParaNum; row++)
	{
		if (Ininf[row]< 0.01f)Ininf[row] = 0.01f;
	}

#define Peak_Scale	Ininf[0] * AScalFactor
#define XPos_Scale	Ininf[1]
#define YPos_Scale	Ininf[2]
#define SigX_Scale	Ininf[3] * SScalFactor
#define SigY_Scale	Ininf[4] * SScalFactor
#define Bakg_Scale	Ininf[5] * BScalFactor

	for (row = 0; row<RegionLen; row++)
	{
		for (col = 0; col<RegionLen; col++)
		{
			// pf goal function
			pixval = pSubImg[row][col]; //subimg[tid][row*RegionLen + col]; // coloffset+

			rowpos = row + 0.5f; // pixel center position
			colpos = col + 0.5f;

			// get model intensity
			tdat = -((colpos - XPos_Scale)*(colpos - XPos_Scale)*SigX_Scale + (rowpos - YPos_Scale)*(rowpos - YPos_Scale)*SigY_Scale);

			IVal = Peak_Scale * __expf(tdat) + Bakg_Scale;

			tdat = IVal - pixval - pixval*(__logf(IVal) - __logf(pixval)); //to avoid negative value in logf function, that will be a disaster

#if(WLE_ENABLE == 1)
			// weighted MLE, improve localization precision, may slightly sacrifice speed
			WLE_Weight = __expf(-((colpos - ROICenter)*(colpos - ROICenter)*WLE_SigmaX1 + (rowpos - ROICenter)*(rowpos - ROICenter)*WLE_SigmaY1));

			TotalLoss = TotalLoss + WLE_Weight*tdat;
#else
			TotalLoss = TotalLoss + tdat;

#endif // WLE_ENABLE

		}
	}

	return TotalLoss;
}


/*
template <int RegionLen>
__device__ void FindMaxPos(float *CurDat, int *MaxPos)
{
int StartPos = RegionLen / 2 - 1;
int cnt = 0;

*MaxPos = RegionLen / 2;

#pragma unroll
for (cnt = 0; cnt < 3; cnt++)
{
if ((CurDat[StartPos + cnt] >= CurDat[StartPos + cnt - 1]) && (CurDat[StartPos + cnt] >= CurDat[StartPos + cnt + 1]))
{
*MaxPos = StartPos + cnt;
}
}
}
*/
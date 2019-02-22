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
__device__ void MLELocalization_DH3D(float subimg[][SharedImageWidth1], float Ininf[][ThreadsPerBlock], float grad[][ThreadsPerBlock], float d0[][ThreadsPerBlock], float D0[], int tid);
template <int RegionLen, int SharedImageWidth1, int FitParaNum>
__device__ void poissonfGradient_DH3D(float subimg[][SharedImageWidth1], float Ininf[][ThreadsPerBlock], float grad[][ThreadsPerBlock], int tid);

template <int RegionLen, int SharedImageWidth1, int FitParaNum>
__device__ void PreEstimation_DH3D(float subimg[][SharedImageWidth1], float Ininf[][ThreadsPerBlock], float *TotalPhoton, int tid);
template <int RegionLen, int SharedImageWidth1, int FitParaNum>
__device__ float MLEFit_TargetF_DH3D(float subimg[][SharedImageWidth1], float Ininf[], int tid);

template <int RegionLen>
__device__ float RestrictAngle(float angle);

// algorithm core codes

template <int RegionLen, int FitParaNum>
__global__ void bfgsMLELoc_DH3D(unsigned short *d_SubRegion, float *d_LocArry, float Offset, float kadc, float QE, int FluoNum)
{
	enum {
		ROIDataLen = (RegionLen*RegionLen),
		ROIWholeLen = (RegionLen*(RegionLen + 1)),

		LoadNum = (ROIDataLen / ThreadsPerBlock),
		LoadResidual = (ROIDataLen%ThreadsPerBlock),

		// avoid bank conflict, if no conflict, +1 can be avoid
		SharedImageWidth1 = (ROIWholeLen + 1),

		HalfRegionLen = (int)(RegionLen / 2),
	};

	float XYLBound = (((float)RegionLen) / 2.0f - 2.0f);
	float XYUBound = (((float)RegionLen) / 2.0f + 2.0f);

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

	float TotalPhoton = 0.0f;
	float CurSNR = 0;

	float CurAngle = 0;
	float tempuse;

	// load image region from global memory, there two method to load
	int loadcnt = 0;
#pragma unroll
	for (cnt = 0; cnt < ThreadsPerBlock; cnt++)
	{
		gMemPos = (BlockOffset + cnt)*ROIWholeLen;
#pragma unroll
		for (loadcnt = 0; loadcnt < LoadNum; loadcnt++)
		{
			ImageRegion[cnt][loadcnt*ThreadsPerBlock + tid] = (d_SubRegion[gMemPos + loadcnt*ThreadsPerBlock + tid] - Offset)*kadc;
		}
	}
	if (tid < LoadResidual)
	{
#pragma unroll
		for (cnt = 0; cnt < ThreadsPerBlock; cnt++)
		{
			gMemPos = (BlockOffset + cnt)*ROIWholeLen;
			ImageRegion[cnt][LoadNum*ThreadsPerBlock + tid] = (d_SubRegion[gMemPos + LoadNum*ThreadsPerBlock + tid] - Offset)*kadc;
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


	PreEstimation_DH3D<RegionLen, SharedImageWidth1, FitParaNum>(ImageRegion, Ininf, &TotalPhoton, tid);

	poissonfGradient_DH3D<RegionLen, SharedImageWidth1, FitParaNum>(ImageRegion, Ininf, grad, tid);

#pragma unroll
	for (rcnt = 0; rcnt < FitParaNum; rcnt++)
	{
		d0[rcnt][tid] = -grad[rcnt][tid];

	}


	MLELocalization_DH3D<RegionLen, SharedImageWidth1, FitParaNum, IterateNum_DH3D, IterateNum_bs_DH3D>(ImageRegion, Ininf, grad, d0, &D0[tid][0], tid);

	Ininf[FitDH3D_Peak][tid] = Ininf[FitDH3D_Peak][tid] * AScalFactor; // PeakPhoton
	Ininf[FitDH3D_SigX][tid] = Ininf[FitDH3D_SigX][tid] * SScalFactor;
	Ininf[FitDH3D_SigY][tid] = Ininf[FitDH3D_SigY][tid] * SScalFactor;
	Ininf[FitDH3D_Bakg][tid] = Ininf[FitDH3D_Bakg][tid] * BScalFactor; // bg

	Ininf[FitDH3D_SigX][tid] = sqrtf(0.5f / Ininf[FitDH3D_SigX][tid]); // sigma
	Ininf[FitDH3D_SigY][tid] = sqrtf(0.5f / Ininf[FitDH3D_SigY][tid]); // sigma


	TotalPhoton = TotalPhoton - Ininf[FitDH3D_Bakg][tid] * RegionLen*RegionLen; // e-
	CurSNR = Ininf[FitDH3D_Peak][tid] / sqrtf(Ininf[FitDH3D_Peak][tid] + Ininf[FitDH3D_Bakg][tid]); // e-

	// convert e- to photon
	Ininf[FitDH3D_Peak][tid] = Ininf[FitDH3D_Peak][tid] / QE;
	Ininf[FitDH3D_Bakg][tid] = Ininf[FitDH3D_Bakg][tid] / QE;
	TotalPhoton = TotalPhoton / QE;

	int CurFrame = ImageRegion[tid][ROIDataLen + 3] * 65536 + ImageRegion[tid][ROIDataLen + 2];

	// note the PeakPhoton and bg are electron
	if ((Ininf[FitDH3D_XPos][tid]>XYLBound) && (Ininf[FitDH3D_XPos][tid]<XYUBound) && (Ininf[FitDH3D_YPos][tid]>XYLBound) && (Ininf[FitDH3D_YPos][tid]<XYUBound)) //
	{
		// for both online and offline localization

		Ininf[FitDH3D_XPos][tid] = Ininf[FitDH3D_XPos][tid] + ImageRegion[tid][ROIDataLen + 0] - HalfRegionLen; // X
		Ininf[FitDH3D_YPos][tid] = Ininf[FitDH3D_YPos][tid] + ImageRegion[tid][ROIDataLen + 1] - HalfRegionLen; // Y

		//because there are two possible angle and sigma combination

		CurAngle = Ininf[FitDH3D_Angl][tid];
		CurAngle = RestrictAngle<RegionLen>(CurAngle);

		if (Ininf[FitDH3D_SigX][tid] < Ininf[FitDH3D_SigY][tid])
		{
			// change elliptical gaussian sigma, the first is longer one, and the angle is its angle
			tempuse = Ininf[FitDH3D_SigX][tid];
			Ininf[FitDH3D_SigX][tid] = Ininf[FitDH3D_SigY][tid];
			Ininf[FitDH3D_SigY][tid] = tempuse;

			// adjust the angle
			if (CurAngle > 0)CurAngle -= Math_PI / 2;
			if (CurAngle < 0)CurAngle += Math_PI / 2;
		}

		pLocArry[gid][Pos_PPho] = Ininf[FitDH3D_Peak][tid];	// peak photon
		pLocArry[gid][Pos_XPos] = Ininf[FitDH3D_XPos][tid];	// may have 0.5 or 1 pixel offset compared with other software
		pLocArry[gid][Pos_YPos] = Ininf[FitDH3D_YPos][tid];	// may have 0.5 or 1 pixel offset compared with other software
		pLocArry[gid][Pos_ZPos] = CurAngle; 	       // may have 0.5 or 1 pixel offset compared with other software
		pLocArry[gid][Pos_SigX] = Ininf[FitDH3D_SigX][tid];	// sigma x
		pLocArry[gid][Pos_SigY] = Ininf[FitDH3D_SigY][tid];	// sigma y
		pLocArry[gid][Pos_TPho] = TotalPhoton;				// total photon
		pLocArry[gid][Pos_Bakg] = Ininf[FitDH3D_Bakg][tid];	// background
		pLocArry[gid][Pos_PSNR] = CurSNR;						// peal snr
		pLocArry[gid][Pos_Frme] = CurFrame; // frame

	}
	else
	{
		// invalid point
		pLocArry[gid][Pos_PPho] = 0; // peak photon
		pLocArry[gid][Pos_XPos] = 0; // may have 0.5 or 1 pixel offset compared with other software
		pLocArry[gid][Pos_YPos] = 0; // may have 0.5 or 1 pixel offset compared with other software
		pLocArry[gid][Pos_ZPos] = 0; // may have 0.5 or 1 pixel offset compared with other software
		pLocArry[gid][Pos_SigX] = 0; // sigma
		pLocArry[gid][Pos_SigY] = 0; // angle
		pLocArry[gid][Pos_TPho] = 0; // total photon
		pLocArry[gid][Pos_Bakg] = 0; // background
		pLocArry[gid][Pos_PSNR] = 0; // peal snr
		pLocArry[gid][Pos_Frme] = 0; // frame
	}
}


template <int RegionLen, int SharedImageWidth1, int FitParaNum, int IterateNum, int IterateNum_bs>
__device__ void MLELocalization_DH3D(float subimg[][SharedImageWidth1], float Ininf[][ThreadsPerBlock], float grad[][ThreadsPerBlock], float d0[][ThreadsPerBlock], float D0[], int tid)
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

		ddat[0] = MLEFit_TargetF_DH3D<RegionLen, SharedImageWidth1, FitParaNum>(subimg, &xd[0], tid);
		ddat[1] = MLEFit_TargetF_DH3D<RegionLen, SharedImageWidth1, FitParaNum>(subimg, &xd[FitParaNum], tid);
		// 
		for (cnt = 0; cnt < IterateNum_bs; cnt++)
		{
			// which part shrink
			xdsel = (ddat[0] < ddat[1]);


			dpos[xdsel] = (dpos[0] + dpos[1])*0.5f; //  /2.0f which one shift

			if (cnt < IterateNum_bs - 1)
			{
				VectorAddMul1<FitParaNum>(&xd[xdsel*FitParaNum], Ininf, d0, dpos[xdsel], tid);// xd=ininf+d0*scale
				ddat[xdsel] = MLEFit_TargetF_DH3D<RegionLen, SharedImageWidth1, FitParaNum>(subimg, &xd[xdsel*FitParaNum], tid);
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

			poissonfGradient_DH3D<RegionLen, SharedImageWidth1, FitParaNum>(subimg, Ininf, grad, tid);

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
__device__ void poissonfGradient_DH3D(float subimg[][SharedImageWidth1], float Ininf[][ThreadsPerBlock], float grad[][ThreadsPerBlock], int tid)
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
		tgradn = MLEFit_TargetF_DH3D<RegionLen, SharedImageWidth1, FitParaNum>(subimg, tIninf, tid);

		tIninf[cnt] = tIninf[cnt] + 0.004f;
		tgradp = MLEFit_TargetF_DH3D<RegionLen, SharedImageWidth1, FitParaNum>(subimg, tIninf, tid);

		grad[cnt][tid] = (tgradp - tgradn)*250.0f; // /0.004f;
		tIninf[cnt] = tIninf[cnt] - 0.002f;
	}
}


template <int RegionLen, int SharedImageWidth1, int FitParaNum>
__device__ void PreEstimation_DH3D(float subimg[][SharedImageWidth1], float Ininf[][ThreadsPerBlock], float *TotalPhoton, int tid)
{
	int rcnt, ccnt;

	float sum_bg = 0;
	//	float PeakDat = 0;

	float(*pSubImg)[RegionLen] = (float(*)[RegionLen])&subimg[tid][0];

	int tPhotonNum = 0;
	// in uint e-
#pragma unroll
	for (rcnt = 0; rcnt < RegionLen; rcnt++)
	{
#pragma unroll
		for (ccnt = 0; ccnt < RegionLen; ccnt++)
		{
			tPhotonNum += pSubImg[rcnt][ccnt];
		}
	}

	*TotalPhoton = tPhotonNum;

	sum_bg = 0;
#pragma unroll
	for (rcnt = 0; rcnt < RegionLen; rcnt++)
	{
#pragma unroll
		for (ccnt = 0; ccnt < RegionLen; ccnt++)
		{
			if ((rcnt == 0) || (ccnt == 0) || (rcnt == RegionLen - 1) || (ccnt == RegionLen - 1))
			{
				sum_bg += pSubImg[rcnt][ccnt];
			}
		}
	}

	sum_bg = sum_bg / (4 * RegionLen - 4); // background



	Ininf[FitDH3D_Peak][tid] = (pSubImg[3][3] - sum_bg)*rAScalFactor; // Amplitude
	Ininf[FitDH3D_XPos][tid] = (float)RegionLen / 2.0f; // sumDatx / sumdat - 0.5f;	// x
	Ininf[FitDH3D_YPos][tid] = (float)RegionLen / 2.0f; // sumDaty / sumdat - 0.5f;	// y
	Ininf[FitDH3D_SigX][tid] = DHParaSigma*rSScalFactor;	// sigma
	Ininf[FitDH3D_SigY][tid] = DHParaSigma*rSScalFactor;	// sigma
	Ininf[FitDH3D_Bakg][tid] = sum_bg*rBScalFactor; // /128.0f;	// bg
	Ininf[FitDH3D_Angl][tid] = 0;	// sigma

}


template <int RegionLen, int SharedImageWidth1, int FitParaNum>
__device__ float MLEFit_TargetF_DH3D(float subimg[][SharedImageWidth1], float Ininf[], int tid)
{
	int row, col;

	float rowpos, colpos;
	float VPoss = 0;
	float IVal;
	float pixval;
	float tdat;

	float(*pSubImg)[RegionLen] = (float(*)[RegionLen])&subimg[tid][0];

#pragma unroll
	for (row = 0; row < FitParaNum - 1; row++)
	{
		// no restriction on angle
		if (Ininf[row]< 0.01f)Ininf[row] = 0.01f;
	}

#define Peak_Scale		Ininf[0] * AScalFactor
#define XPos_Scale		Ininf[1]
#define YPos_Scale		Ininf[2]
#define SigX_Scale		Ininf[3] * SScalFactor
#define SigY_Scale		Ininf[4] * SScalFactor
#define Bakg_Scale		Ininf[5] * BScalFactor
#define Angl_Scale		Ininf[6]


	for (row = 0; row<RegionLen; row++)
	{
		for (col = 0; col<RegionLen; col++)
		{
			// pf goal function
			pixval = pSubImg[row][col]; //subimg[tid][row*RegionLen + col]; // coloffset+

			rowpos = row + 0.5f; // -0.5
			colpos = col + 0.5f; // -0.5

			// get model intensity
			tdat = -(powf((colpos - XPos_Scale)* __cosf(Angl_Scale) + (rowpos - YPos_Scale)* __sinf(Angl_Scale), 2) * SigX_Scale + powf(-(colpos - XPos_Scale)* __sinf(Angl_Scale) + (rowpos - YPos_Scale)* __cosf(Angl_Scale), 2) * SigY_Scale);

			IVal = Peak_Scale*__expf(tdat) + Bakg_Scale;

			tdat = IVal - pixval - pixval*(__logf(IVal) - __logf(fmaxf(pixval, 1.0f))); //to avoid negative value in logf function, that will be a disaster

			VPoss = VPoss + tdat;
		}
	}

	return VPoss;
}

template <int RegionLen>
__device__ float RestrictAngle(float angle)
{
	while (angle > Math_PI)
	{
		// positive large
		angle -= Math_PI;
	}
	while (angle < -Math_PI)
	{
		// negative large
		angle += Math_PI;
	}
	return angle;
}


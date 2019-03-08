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

#include "DH3D_MoleculePair.h"



__device__ int IsValueSimilar(float Val1, float Val2, float DiffRatio);

__device__ void GetPairAngle(float *oAngle, float *oX0, float *oY0, float px1, float py1, float px2, float py2);


#define MaxSimilarFluoNum	5

__global__ void gpuPairMolecules(float *d_LocArry, float *d_oLocArry, float *d_oPosArry, int *d_ValidFluoNum, float MeanDistance, float DistanceTh, int FluoNum, int PosArryEnable, float p4, float p3, float p2, float p1, float p0, int RotateType)
{

	//    int tid = threadIdx.x;
	int gid = threadIdx.x + blockDim.x*blockIdx.x;


	float(*pLocArry)[OutParaNumGS2D]; // for 9 parameter array
	float(*poLocArry10)[OutParaNumGS2D]; // for 9 parameter array
	float(*poPosArry10)[OutParaNumGS2D]; // for 9 parameter array

	pLocArry = (float(*)[OutParaNumGS2D])d_LocArry;
	poLocArry10 = (float(*)[OutParaNumGS2D])d_oLocArry;
	poPosArry10 = (float(*)[OutParaNumGS2D])d_oPosArry;

	float PPhotonList[MaxSimilarFluoNum];
	float SXList[MaxSimilarFluoNum];
	float SYList[MaxSimilarFluoNum];
	float SigmaList[MaxSimilarFluoNum];
	float TotalPnList[MaxSimilarFluoNum];
	float BkgList[MaxSimilarFluoNum];
	float SNRList[MaxSimilarFluoNum];


	float PeakPhoton1;
	float PeakPhoton2;

	float XPos1;
	float XPos2;

	float YPos1;
	float YPos2;	

	float Sigma1;
	float Sigma2;
	
	float Background1;
	float Background2;

	float TotalPhoton1;
	float TotalPhoton2;

	float SNR1;
	float SNR2;

	int Frame1;
	int Frame2;

	int cnt = 0;

	float CurDistance = 0;

	int IsValid0 = 0;
	int IsValid1 = 0;
	int IsValid2 = 0;
	int IsValid3 = 0;

	float CurAngle;
	float CenterX0, CenterY0;

	int CurPos = 0;
	float CurZPos = 0;

	int SimilarFluoNum = 0;


	if (gid < FluoNum)
	{
		PeakPhoton1 = pLocArry[gid][Pos_PPho]; // peak photon 
		XPos1 = pLocArry[gid][Pos_XPos];
		YPos1 = pLocArry[gid][Pos_YPos];
		Sigma1 = pLocArry[gid][Pos_SigX];
		TotalPhoton1 = pLocArry[gid][Pos_TPho];
		Background1 = pLocArry[gid][Pos_Bakg];
		SNR1 = pLocArry[gid][Pos_PSNR];
		Frame1 = pLocArry[gid][Pos_Frme];


		for (cnt = 0; cnt < FluoNum; cnt++)
		{
			PeakPhoton2 = pLocArry[cnt][Pos_PPho]; // peak photon 
			XPos2 = pLocArry[cnt][Pos_XPos];
			YPos2 = pLocArry[cnt][Pos_YPos];
			Sigma2 = pLocArry[cnt][Pos_SigX];
			TotalPhoton2 = pLocArry[cnt][Pos_TPho];
			Background2 = pLocArry[cnt][Pos_Bakg];
			SNR2 = pLocArry[cnt][Pos_PSNR];
			Frame2 = pLocArry[cnt][Pos_Frme];


			CurDistance = sqrtf((XPos1 - XPos2)*(XPos1 - XPos2) + (YPos1 - YPos2)*(YPos1 - YPos2));

			IsValid0 = fabsf(CurDistance - MeanDistance) < DistanceTh;

			IsValid1 = IsValueSimilar(PeakPhoton1, PeakPhoton2, 1.0f);
			IsValid2 = IsValueSimilar(TotalPhoton1, TotalPhoton2, 2.0f);
			IsValid3 = IsValueSimilar(Sigma1, Sigma2, 0.8f);



			// this may be removed to the end of the for, to maximize performance, just record paired point first
			// a point could be paired to several other molecules, how to solve?
			// may be store maximul 3 pair pos, then choose the best similar


			// 
			if ((IsValid0) & (IsValid1 & IsValid2) & IsValid3 && (Frame1 == Frame2) && (gid < cnt))
			{
				// a&b, vs b&a, so judge a,b order can solve this
				PPhotonList[SimilarFluoNum] = PeakPhoton1 + PeakPhoton2;
				SXList[SimilarFluoNum] = XPos2;
				SYList[SimilarFluoNum] = YPos2;
				SigmaList[SimilarFluoNum] = (Sigma1 + Sigma2)*0.5f;
				TotalPnList[SimilarFluoNum] = TotalPhoton1 + TotalPhoton2;
				BkgList[SimilarFluoNum] = Background1 + Background2;
				SNRList[SimilarFluoNum] = (SNR1 + SNR2)*0.5f;

//				printf("dat:%f %f %f %d\n", PeakPhoton2, XPos2, YPos2, gid);

				SimilarFluoNum++;
				if (SimilarFluoNum >= MaxSimilarFluoNum)break;

			}
			if (Frame2 > Frame1)
			{
				break;
			}

		}
//		if (SimilarFluoNum > 1)printf("multi:%d %d--%f %f %f--%f %f %f--%f %f %f\n", gid, SimilarFluoNum, PeakPhoton1, TotalPhoton1, Sigma1, PPhotonList[0], TotalPnList[0], SigmaList[0], PPhotonList[1], TotalPnList[1], SigmaList[1]);

		for (cnt = 0; cnt < SimilarFluoNum; cnt++)
		{
			// remove it from the for, to reduce warp divergence
			// can be compared to choose the best matched molecule
			CurPos = atomicAdd(d_ValidFluoNum, 1); // cur molecules pos

			GetPairAngle(&CurAngle, &CenterX0, &CenterY0, XPos1, YPos1, SXList[cnt], SYList[cnt]);

			CurZPos = p4*powf(CurAngle, 4) + p3*powf(CurAngle, 3) + p2*powf(CurAngle, 2) + p1*CurAngle + p0;


			poLocArry10[CurPos][Pos_PPho] = PPhotonList[cnt];   // peak photon
			poLocArry10[CurPos][Pos_XPos] = CenterX0; 			// may have 0.5 or 1 pixel offset compared with other software
			poLocArry10[CurPos][Pos_YPos] = CenterY0; 			// may have 0.5 or 1 pixel offset compared with other software
			poLocArry10[CurPos][Pos_ZPos] = CurZPos;  			// may have 0.5 or 1 pixel offset compared with other software
			poLocArry10[CurPos][Pos_SigX] = SigmaList[cnt];    			// sigma x
			poLocArry10[CurPos][Pos_SigY] = SigmaList[cnt];    			// sigma y
			poLocArry10[CurPos][Pos_TPho] = TotalPnList[cnt]; 		// total photon
			poLocArry10[CurPos][Pos_Bakg] = BkgList[cnt]; 		// background
			poLocArry10[CurPos][Pos_PSNR] = SNRList[cnt]; 			// peal snr
			poLocArry10[CurPos][Pos_Frme] = Frame1;   			// frame

			// get the pair molecules information
			if (PosArryEnable)
			{
				CurDistance = sqrtf((XPos1 - SXList[cnt])*(XPos1 - SXList[cnt]) + (YPos1 - SYList[cnt])*(YPos1 - SYList[cnt]));

				poPosArry10[CurPos][Pos_PPho0] = PeakPhoton1; // peak photon
				poPosArry10[CurPos][Pos_XPos0] = XPos1; // may have 0.5 or 1 pixel offset compared with other software
				poPosArry10[CurPos][Pos_YPos0] = YPos1; // may have 0.5 or 1 pixel offset compared with other software
				poPosArry10[CurPos][Pos_PPho1] = PPhotonList[cnt];  // may have 0.5 or 1 pixel offset compared with other software
				poPosArry10[CurPos][Pos_XPos1] = SXList[cnt];    // molecule 1 x pos
				poPosArry10[CurPos][Pos_YPos1] = SYList[cnt];    // molecule 1 y pos
				poPosArry10[CurPos][Pos_ZPos0] = CurZPos; // molecule 2 x pos
				poPosArry10[CurPos][Pos_Dista] = CurDistance; // molecule 2 y pos
				poPosArry10[CurPos][Pos_Angle] = CurAngle; // molecules pair angle
				poPosArry10[CurPos][Pos_Frame] = Frame1;   // frame

			}
		}
	}
}



__global__ void gpuCalcDistanceDistribution(float *d_PosArry, int *d_DistanceDistrib, int FluoNum)
{
	//    int tid = threadIdx.x;
	int gid = threadIdx.x + blockDim.x*blockIdx.x;

	float(*pPosArry10)[OutParaNumGS2D]; // for 9 parameter array
	pPosArry10 = (float(*)[OutParaNumGS2D])d_PosArry;

	float CurDistance = 0;
	int HistPos = 0;

	if (gid < FluoNum)
	{
		CurDistance = pPosArry10[gid][Pos_Dista];

		HistPos = CurDistance / DistanceGap + 0.5f;
		if (HistPos < 0)HistPos = 0;
		if (HistPos > DistanceHistLen - 1)HistPos = DistanceHistLen - 1;

		atomicAdd(&d_DistanceDistrib[HistPos], 1);

	}

}




__device__ int IsValueSimilar(float Val1, float Val2, float DiffRatio)
{
	float MinVal = fminf(Val1, Val2);
	
	float DiffVal = fabsf(Val1 - Val2) / MinVal;

	int IsValid = 0;

	if (DiffVal < DiffRatio)
	{
		IsValid = 1;
	}

	return IsValid;
}

__device__ void GetPairAngle(float *oAngle, float *oX0, float *oY0, float px1, float py1, float px2, float py2)
{
	float upx, upy;
	float dpx, dpy;

	if (py1 > py2)
	{
		upx = px1;
		upy = py1;
		dpx = px2;
		dpy = py2;
	}
	else
	{
		upx = px2;
		upy = py2;
		dpx = px1;
		dpy = py1;
	}

	float Difx = upx - dpx;
	float Dify = upy - dpy;

	// radian is better than degree, since radian is small, the degree of angle is much bigger
	*oAngle = acosf(Difx*rsqrtf(Difx*Difx + Dify*Dify)); // arc difx/sqrt(difx^2+dify^2);  180 / 3.14159265f*
	*oX0 = (upx + dpx) / 2;
	*oY0 = (upy + dpy) / 2;
}


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

#include "ConsecutiveFit.h"


// switch between two array
__device__ float * GetPartSelID(int *PartSel, int *PartID, int CurPos, int LastFluoNum, float *d_LastLocArry, float *d_CurLocArry)
{
	float* oLocArry;

	*PartSel = CurPos >= LastFluoNum;

	if (*PartSel == 0)
	{
		*PartID = CurPos;
		oLocArry = d_LastLocArry;
	}
	else
	{
		*PartID = CurPos - LastFluoNum;
		oLocArry = d_CurLocArry;

	}
	return oLocArry;
}

__global__ void gpuFindConsecutiveMolecules(float *d_LastLocArry, int LastFluoNum, float *d_CurLocArry, int CurFluoNum, int *d_ForwardLinkID, int *d_BackwardLinkID, float DistanceTh_nm, float PixelSize, float QE)
{
	//    int tid = threadIdx.x;
	int gid = threadIdx.x + blockDim.x*blockIdx.x;
	int gid_0 = blockDim.x*blockIdx.x;

	float(*pLocArry)[OutParaNumGS2D]; // for 9 parameter array

	int TotalFluoNum = LastFluoNum + CurFluoNum;

	int PartSel;
	int PartID;


	int PairID;
	float XPos1;
	float XPos2;

	float YPos1;
	float YPos2;


	int Frame1 = 0;
	int Frame2 = 0;

	int cnt = 0;

	float DistanceTh2 = 0; // square of distance

	float CurDistance2 = 0;

	int SimilarFluoNum = 0;


	if (gid < TotalFluoNum)
	{
		
		pLocArry = (float(*)[OutParaNumGS2D])GetPartSelID(&PartSel, &PartID, gid, LastFluoNum, d_LastLocArry, d_CurLocArry);


		XPos1 = pLocArry[PartID][Pos_XPos];
		YPos1 = pLocArry[PartID][Pos_YPos];
		Frame1 = pLocArry[PartID][Pos_Frme];

		// can be revised to be loc precision

//		DistanceTh2 = 0.5f; // 0.5 pixel

		DistanceTh2 = DistanceTh_nm / PixelSize;

		DistanceTh2 = DistanceTh2*DistanceTh2;

		if ((XPos1 > 5) && (YPos1 > 5))
		{
			// molecule's position in next frame is larger
			for (cnt = gid_0; (cnt < TotalFluoNum) && (Frame2 <= Frame1 + 1); cnt++)
			{
				//
				pLocArry = (float(*)[OutParaNumGS2D])GetPartSelID(&PartSel, &PartID, cnt, LastFluoNum, d_LastLocArry, d_CurLocArry);

				//
				XPos2 = pLocArry[PartID][Pos_XPos];
				YPos2 = pLocArry[PartID][Pos_YPos];
				Frame2 = pLocArry[PartID][Pos_Frme];


				CurDistance2 = (XPos1 - XPos2)*(XPos1 - XPos2) + (YPos1 - YPos2)*(YPos1 - YPos2);

				// find the closest in next frame
				if ((CurDistance2 < DistanceTh2) && (Frame2 - Frame1 == 1))
				{
					PairID = cnt;
					SimilarFluoNum++;
					DistanceTh2 = CurDistance2; // update distance th for find the closest
				}
			}

			if (SimilarFluoNum > 0)
			{
				// if no molecule is paired, it's 0
				d_BackwardLinkID[gid] = PairID;
				d_ForwardLinkID[PairID] = 1; // has foward
			}
			// else is zero
		}
		else
		{
			// invalid molecules
			d_ForwardLinkID[gid] = -2;
		}
	}
}


__global__ void gpuConsecutiveFit(float *d_LastLocArry, int LastFluoNum, float *d_CurLocArry, int CurFluoNum, int *d_ForwardLinkID, int *d_BackwardLinkID, int ConsecFitFluoNum, float QE)
{
	//    int tid = threadIdx.x;
	int gid = threadIdx.x + blockDim.x*blockIdx.x;

	float(*pLocArry)[OutParaNumGS2D]; // for 9 parameter array

	int TotalFluoNum = LastFluoNum + CurFluoNum;

	int PartSel;
	int PartID;


	int CurId = 0;
	int NextValid = 0;


	float wXPos = 0;
	float wYPos = 0;
	float wZPos = 0;

	float wSigmaX = 0;
	float wSigmaY = 0;

	float SumPeakPhoton = 0;
	float SumTotalPhoton = 0;
	float SumBackground = 0;
	float SumSNR = 0;

	float SumPrec = 0;

	int ConsecNum = 0;

	if (gid < ConsecFitFluoNum)
	{
		// the first fluo is valid fluo and don't have a forward
		if (d_ForwardLinkID[gid] == 0)
		{
			CurId = gid;
			NextValid = 1;

			while (NextValid)
			{

				//
				pLocArry = (float(*)[OutParaNumGS2D])GetPartSelID(&PartSel, &PartID, CurId, LastFluoNum, d_LastLocArry, d_CurLocArry);
				//

				float LocError = pLocArry[PartID][Pos_CrbX]; // nm
				float LocPrec = 10 / LocError; // 1 / LocPrec; convert from loc error to precision

				// localization precision (1/loc err) weighted position, compared to nomal average
				SumPeakPhoton += pLocArry[PartID][Pos_PPho];
				SumTotalPhoton += pLocArry[PartID][Pos_TPho];
				SumBackground += pLocArry[PartID][Pos_Bakg];

				// loc precision weighted average
				wXPos += pLocArry[PartID][Pos_XPos] * LocPrec;
				wYPos += pLocArry[PartID][Pos_YPos] * LocPrec;
				wZPos += pLocArry[PartID][Pos_ZPos] * LocPrec;

				wSigmaX += pLocArry[PartID][Pos_SigX] * LocPrec;
				wSigmaY += pLocArry[PartID][Pos_SigY] * LocPrec;

				SumPrec += LocPrec;


				CurId = d_BackwardLinkID[CurId];
				NextValid = CurId > 0;

				ConsecNum++;
			}
			// weighted average
			wXPos /= SumPrec;
			wYPos /= SumPrec;
			wZPos /= SumPrec;

			wSigmaX /= SumPrec;
			wSigmaY /= SumPrec;


			// update consecutive fit
			SumSNR = SumPeakPhoton*QE / sqrtf(SumPeakPhoton*QE + SumBackground*QE);

			//
			pLocArry = (float(*)[OutParaNumGS2D])GetPartSelID(&PartSel, &PartID, gid, LastFluoNum, d_LastLocArry, d_CurLocArry);


			pLocArry[PartID][Pos_PPho] = SumPeakPhoton;
			pLocArry[PartID][Pos_XPos] = wXPos;
			pLocArry[PartID][Pos_YPos] = wYPos;
			pLocArry[PartID][Pos_ZPos] = wZPos;
			pLocArry[PartID][Pos_SigX] = wSigmaX;
			pLocArry[PartID][Pos_SigY] = wSigmaY;
			pLocArry[PartID][Pos_TPho] = SumTotalPhoton;
			pLocArry[PartID][Pos_Bakg] = SumBackground;
			pLocArry[PartID][Pos_PSNR] = SumSNR;
		}
	}
}


__global__ void gpuRemoveConsecutiveFluo(float *d_LastLocArry, int LastFluoNum, float *d_CurLocArry, int CurFluoNum, int *d_ForwardLinkID, int *d_BackwardLinkID, int ConsecFitFluoNum)
{
	//    int tid = threadIdx.x;
	int gid = threadIdx.x + blockDim.x*blockIdx.x;
	float(*pLocArry)[OutParaNumGS2D]; // for 9 parameter array

	int TotalFluoNum = LastFluoNum + CurFluoNum;

	int PartSel;
	int PartID;

	int cnt;

	int NextID;

	if (gid < ConsecFitFluoNum)
	{
		// the first fluo is valid fluo and don't have a forward
		if (d_ForwardLinkID[gid] == 0)
		{

			NextID = d_BackwardLinkID[gid];

			while (NextID > 0)
			{
				//
				pLocArry = (float(*)[OutParaNumGS2D])GetPartSelID(&PartSel, &PartID, NextID, LastFluoNum, d_LastLocArry, d_CurLocArry);

				//
#pragma unroll
				for (cnt = 0; cnt < OutParaNumGS2D; cnt++)
				{
					pLocArry[PartID][cnt] = 0;
				}
				NextID = d_BackwardLinkID[NextID];
			}
		}
	}
}


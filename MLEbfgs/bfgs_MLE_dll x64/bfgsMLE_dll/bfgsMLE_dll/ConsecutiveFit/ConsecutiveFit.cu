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


__global__ void gpuFindConsecutiveMolecules(float * d_LocArry_ConsecFit, int TotalFluoNum, int *d_ForwardLinkID, int *d_BackwardLinkID, float DistanceTh_Pixel)
{
	//    int tid = threadIdx.x;
	int gid = threadIdx.x + blockDim.x*blockIdx.x;
	int gid_0 = blockDim.x*blockIdx.x;

	float(*pLocArry)[OutParaNumGS2D] = (float(*)[OutParaNumGS2D])d_LocArry_ConsecFit;


	if (gid < TotalFluoNum)
	{
		int XPos1 = pLocArry[gid][Pos_XPos];
		int YPos1 = pLocArry[gid][Pos_YPos];
		int Frame1 = pLocArry[gid][Pos_Frme];

		int XPos2 = 0;
		int YPos2 = 0;
		int Frame2 = 0;

		int PairedID;

		int SimilarFluoNum = 0;

		float DistanceTh2 = DistanceTh_Pixel;

		DistanceTh2 = DistanceTh2*DistanceTh2;

		float CurDistance2 = 0;

		// molecule's position in next frame is larger
		for (int cnt = gid_0; (cnt < TotalFluoNum) && (Frame2 <= Frame1 + 1); cnt++)
		{

			XPos2 = pLocArry[cnt][Pos_XPos];
			YPos2 = pLocArry[cnt][Pos_YPos];
			Frame2 = pLocArry[cnt][Pos_Frme];


			CurDistance2 = (XPos1 - XPos2)*(XPos1 - XPos2) + (YPos1 - YPos2)*(YPos1 - YPos2);

			// find the closest in next frame
			if ((CurDistance2 < DistanceTh2) && (Frame2 - Frame1 == 1))
			{
				PairedID = cnt;

				SimilarFluoNum++;

				DistanceTh2 = CurDistance2; // update distance th for find the closest
			}
		}

		if (SimilarFluoNum > 0)
		{
			// if no molecule is paired, it's 0
			d_BackwardLinkID[gid] = PairedID;
			d_ForwardLinkID[PairedID] = 1; // has a foward
		}
		else
		{
			d_BackwardLinkID[gid] = -2;
		}
	}
}



__global__ void gpuConsecutiveFit(float * d_LocArry_ConsecFit, int TotalFluoNum, int FluoNum_LastGroup, int *d_ForwardLinkID, int *d_BackwardLinkID, int *d_OntimeDistrib, int *d_ValidFluoNum, float QE)
{
	//    int tid = threadIdx.x;
	int gid = threadIdx.x + blockDim.x*blockIdx.x;


	float(*pLocArry)[OutParaNumGS2D] = (float(*)[OutParaNumGS2D])d_LocArry_ConsecFit;


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


	if (gid < TotalFluoNum)
	{
		// the first fluo is valid fluo and don't have a forward
		if (d_ForwardLinkID[gid] == 0)
		{
			int CurId = gid;
			int NextValid = 1;

			int ConsecNum = 0;

			while (NextValid)
			{

				float LocError = pLocArry[CurId][Pos_CrbX]; // nm
				float LocPrec = 10 / LocError; // 1 / LocPrec; convert from loc error to precision


				// localization precision (1/loc err) weighted position, compared to nomal average
				SumPeakPhoton += pLocArry[CurId][Pos_PPho];
				SumTotalPhoton += pLocArry[CurId][Pos_TPho];
				SumBackground += pLocArry[CurId][Pos_Bakg];

				// loc precision weighted average
				wXPos += pLocArry[CurId][Pos_XPos] * LocPrec;
				wYPos += pLocArry[CurId][Pos_YPos] * LocPrec;
				wZPos += pLocArry[CurId][Pos_ZPos] * LocPrec;

				wSigmaX += pLocArry[CurId][Pos_SigX] * LocPrec;
				wSigmaY += pLocArry[CurId][Pos_SigY] * LocPrec;

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
			SumSNR = sqrtf(QE) * SumPeakPhoton / sqrtf(SumPeakPhoton + SumBackground);


			pLocArry[gid][Pos_PPho] = SumPeakPhoton;
			pLocArry[gid][Pos_XPos] = wXPos;
			pLocArry[gid][Pos_YPos] = wYPos;
			pLocArry[gid][Pos_ZPos] = wZPos;
			pLocArry[gid][Pos_SigX] = wSigmaX;
			pLocArry[gid][Pos_SigY] = wSigmaY;
			pLocArry[gid][Pos_TPho] = SumTotalPhoton;
			pLocArry[gid][Pos_Bakg] = SumBackground;
			pLocArry[gid][Pos_PSNR] = SumSNR;


			if (gid >= FluoNum_LastGroup)
			{
				int CurOntime = ConsecNum - 1;

				CurOntime = max(CurOntime, 0);
				CurOntime = min(CurOntime, MaxOnTimeConsecutiveNum - 1);
				

				atomicAdd(d_ValidFluoNum, 1);
				atomicAdd(&d_OntimeDistrib[CurOntime], 1);

			}
		}
	}
}


__global__ void gpuRemoveConsecutiveFluo(float * d_LocArry_ConsecFit, int TotalFluoNum, int *d_ForwardLinkID, int *d_BackwardLinkID)
{
	//    int tid = threadIdx.x;
	int gid = threadIdx.x + blockDim.x*blockIdx.x;


	float(*pLocArry)[OutParaNumGS2D] = (float(*)[OutParaNumGS2D])d_LocArry_ConsecFit;


	if (gid < TotalFluoNum)
	{
		// the first fluo is valid fluo and don't have a forward
		if (d_ForwardLinkID[gid] == 0)
		{

			int NextID = d_BackwardLinkID[gid];

			while (NextID > 0)
			{
				// remove consecutive molecules
#pragma unroll
				for (int cnt = 0; cnt < OutParaNumGS2D; cnt++)
				{
					pLocArry[NextID][cnt] = 0;
				}

				NextID = d_BackwardLinkID[NextID];
			}
		}
	}
}


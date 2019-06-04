#include "ConsecutiveFluoMerge.h"



__global__ void gpuFindConsecutiveFilterPair(float *d_LocArry, int *d_ForwardLinkID, int *d_BackwardLinkID, float Distance_th_pixel, int FluoNum)
{
	//    int tid = threadIdx.x;

	int gid = blockDim.x*blockIdx.x + threadIdx.x;

	// all thread in a warp read from the same memory to improve efficiency
	int gid_0 = gid / 32 * 32; 


	float(*pLocArry)[OutParaNumGS2D]; //
	pLocArry = (float(*)[OutParaNumGS2D])d_LocArry;


	float XPos1;
	float XPos2;

	float YPos1;
	float YPos2;

	int Frame1 = 0;
	int Frame2 = 0;

	int cnt = 0;

	float DistanceTh2; // square of distance

	float CurDistance2;

	int SimilarFluoNum = 0;

	int PairID = 0;

	if (gid < FluoNum)
	{

		XPos1 = pLocArry[gid][Pos_XPos];
		YPos1 = pLocArry[gid][Pos_YPos];
		Frame1 = pLocArry[gid][Pos_Frme];


		// square of distance, to avoid sqrt
		DistanceTh2 = Distance_th_pixel*Distance_th_pixel;

		if ((XPos1 > 5) && (YPos1 > 5))
		{
			// molecule's position in next frame is larger
			for (cnt = gid_0; (Frame2 <= Frame1 + 1) && (cnt < FluoNum); cnt++)
			{
				XPos2 = pLocArry[cnt][Pos_XPos];
				YPos2 = pLocArry[cnt][Pos_YPos];
				Frame2 = pLocArry[cnt][Pos_Frme];


				CurDistance2 = (XPos1 - XPos2)*(XPos1 - XPos2) + (YPos1 - YPos2)*(YPos1 - YPos2);

				// find the closest in next frame
				if ((CurDistance2 < DistanceTh2) && (Frame2 - Frame1 == 1))
				{
					PairID = cnt;
					DistanceTh2 = CurDistance2; // update distance th for find the closest

					SimilarFluoNum++;
				}
			}

			if (SimilarFluoNum > 0)
			{
				d_BackwardLinkID[gid] = PairID;
				d_ForwardLinkID[PairID] = 1; // has foward
			}
			else
			{
				// if don't have a backword consecutive fluo
				d_BackwardLinkID[gid] = -1;
			}
		}
		else
		{
			// invalid molecules
			d_ForwardLinkID[gid] = -2;
			d_BackwardLinkID[gid] = -2;

		}
	}
}



__global__ void gpuConsecutiveFit(float *d_LocArry, int *d_ForwardLinkID, int *d_BackwardLinkID, float QE, int FluoNum)
{
	//    int tid = threadIdx.x;
	int gid = threadIdx.x + blockDim.x*blockIdx.x;


	float(*pLocArry)[OutParaNumGS2D] = (float(*)[OutParaNumGS2D])d_LocArry;


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


	if (gid < FluoNum)
	{
		// the first fluo is valid fluo and don't have a forward
		if (d_ForwardLinkID[gid] == 0)
		{
			CurId = gid;
			NextValid = 1;

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

			}
			// weighted average
			wXPos /= SumPrec;
			wYPos /= SumPrec;
			wZPos /= SumPrec;

			wSigmaX /= SumPrec;
			wSigmaY /= SumPrec;


			// update consecutive fit
			SumSNR = SumPeakPhoton*QE / sqrtf(SumPeakPhoton*QE + SumBackground*QE);


			pLocArry[gid][Pos_PPho] = SumPeakPhoton;
			pLocArry[gid][Pos_XPos] = wXPos;
			pLocArry[gid][Pos_YPos] = wYPos;
			pLocArry[gid][Pos_ZPos] = wZPos;
			pLocArry[gid][Pos_SigX] = wSigmaX;
			pLocArry[gid][Pos_SigY] = wSigmaY;
			pLocArry[gid][Pos_TPho] = SumTotalPhoton;
			pLocArry[gid][Pos_Bakg] = SumBackground;
			pLocArry[gid][Pos_PSNR] = SumSNR;
		}
	}
}



__global__ void gpuRemoveConsecutiveFluo_KeepFirst(float *d_LocArry, int *d_ForwardLinkID, int *d_BackwardLinkID, int FluoNum)
{
	//    int tid = threadIdx.x;
	int gid = blockDim.x*blockIdx.x + threadIdx.x;


	float(*pLocArry)[OutParaNumGS2D]; // for 9 parameter array
	pLocArry = (float(*)[OutParaNumGS2D])d_LocArry;


	int cnt;

	int NextID;

	if (gid < FluoNum)
	{
		// the first fluo is valid fluo and don't have a forward, may have a backword
		if (d_ForwardLinkID[gid] == 0)
		{
			NextID = d_BackwardLinkID[gid];

			// remove fluo in a sequency if it not the first id
			while (NextID > 0)
			{
#pragma unroll
				for (cnt = 0; cnt < OutParaNumGS2D; cnt++)
				{
					pLocArry[NextID][cnt] = 0;
				}

				NextID = d_BackwardLinkID[NextID];
			}
		}
	}
}

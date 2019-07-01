#include "DH3D_MoleculePair.h"



__global__ void gpu_MoleculePair(float *d_LocArry, int FluoNum, int *d_PairID, int * d_ValidoFluoNum, float MeanDistance, float DistanceVaryTh)
{
	float(*pLocArry)[OutParaNumGS2D] = (float(*)[OutParaNumGS2D])d_LocArry;

	int(*pPairID)[PAIR_ID_LEN] = (int(*)[PAIR_ID_LEN])d_PairID;


	int gid = threadIdx.x + blockDim.x*blockIdx.x;
	int gid_0 = blockDim.x*blockIdx.x;

	if (gid < FluoNum)
	{
		float XPos1 = pLocArry[gid][Pos_XPos];
		float YPos1 = pLocArry[gid][Pos_YPos];

		float Photon1 = pLocArry[gid][Pos_TPho];

		int Frame1 = pLocArry[gid][Pos_Frme];

		if (Photon1 > 1.0f)
		{

			float XPos2 = 0;
			float YPos2 = 0;

			float Photon2 = 0;

			int Frame2 = 0;


			float LastDistance = DistanceVaryTh;
			int PairedId = -1;


			for (int cnt = gid_0; (cnt < FluoNum) && (Frame2 <= Frame1); cnt++)
			{

				XPos2 = pLocArry[cnt][Pos_XPos];
				YPos2 = pLocArry[cnt][Pos_YPos];

				Photon2 = pLocArry[cnt][Pos_TPho];
				Photon2 = max(Photon2, 1.0f);

				Frame2 = pLocArry[cnt][Pos_Frme];


				float MinPhoton = min(Photon1, Photon2);
				float MaxPhoton = max(Photon1, Photon2);

				float PhotonDiff = MaxPhoton / MinPhoton;


				if ((PhotonDiff <= 2.5f) && (cnt > gid) && (Frame2 == Frame1))
				{

					float CurDistance = __fsqrt_rz((XPos1 - XPos2)*(XPos1 - XPos2) + (YPos1 - YPos2)*(YPos1 - YPos2));
					float CurDistance1 = fabsf(CurDistance - MeanDistance);

					if (CurDistance1 < LastDistance)
					{
						LastDistance = CurDistance1;
						PairedId = cnt;
					}
				}
			}

			if (PairedId > 0)
			{
				int CurStorPos = atomicAdd(d_ValidoFluoNum, 1);
				pPairID[CurStorPos][PAIR_ID_1st] = gid;
				pPairID[CurStorPos][PAIR_ID_2nd] = PairedId;

//				printf("find pair: %d %d\n", gid, PairedId);

			}
		}
	}
}

/*
RotateMode : 0: angle from 0 to pi, pi/2 is the focal plane
RotateMode : 1: angle from -pi/2 to pi/2, 0 is the focal plane
p4:p0 calibration curve
*/

__global__ void gpu_MoleculeMerge(float *d_LocArry, float *d_oLocArry, int *d_PairID, int ValidFluoNum, int RotateMode, float p4, float p3, float p2, float p1, float p0)
{
	float(*pLocArry)[OutParaNumGS2D] = (float(*)[OutParaNumGS2D])d_LocArry;

	float(*poLocArry)[OutParaNumGS2D] = (float(*)[OutParaNumGS2D])d_oLocArry;

	int(*pPairID)[PAIR_ID_LEN] = (int(*)[PAIR_ID_LEN])d_PairID;


	int gid = threadIdx.x + blockDim.x*blockIdx.x;
	int gid_0 = blockDim.x*blockIdx.x;

	if (gid < ValidFluoNum)
	{
		int ID1 = pPairID[gid][PAIR_ID_1st];
		int ID2 = pPairID[gid][PAIR_ID_2nd];


		float x1 = pLocArry[ID1][Pos_XPos];
		float y1 = pLocArry[ID1][Pos_YPos];

		float x2 = pLocArry[ID2][Pos_XPos];
		float y2 = pLocArry[ID2][Pos_YPos];


		float dx = x1 - x2;
		float dy = y1 - y2;

		float distance = sqrtf(dx*dx + dy*dy);
		float SinA = dy / distance;
		float CosA = dx / distance;

		// don't directly use the coordinate system, but keep the same with eye seen
		float Section = -dx*dy;

		// measured in radian
		float Angle = asinf(fabsf(SinA));


		if (Section < 0)
		{
			if (RotateMode == 0)
			{
				Angle = 3.141592654f - Angle;
			}
			else
			{
				Angle = 0 - Angle;
			}
		}

		float ZPos = powf(Angle, 4)*p4 + powf(Angle, 3)* p3 + powf(Angle, 2)* p2 + Angle*p1 + p0;

//		printf("find pair: %d %d\n", ID1, ID2);

		for (int cnt = 0; cnt < OutParaNumGS2D; cnt++)
		{
			poLocArry[gid][cnt] = (pLocArry[ID1][cnt] + pLocArry[ID2][cnt]) / 2;
		}

		poLocArry[gid][Pos_ZPos] = ZPos;

		// get z from angle

//		printf("oxy: %f %f  %f %f %f %f\n", poLocArry[gid][Pos_XPos], poLocArry[gid][Pos_YPos], Angle, distance, Section, poLocArry[gid][Pos_Frme]);

	}

}



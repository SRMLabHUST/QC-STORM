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

#include "statisticsInfo.h"


__global__ void CalcLocResultStatInf(float *d_LocArry, int *d_Hist_TotalPhoton, int *d_Hist_LocPrecisionXY, int *d_Hist_PeakSNR, int *d_Hist_PSFWidth, float *d_ValidNum, float *d_LocPrecSum, float *d_SNRSum, float *d_PSFWSum, float *d_BgSum, int FluoNum)
{

	//    int tid = threadIdx.x;
	int gid = threadIdx.x + blockDim.x*blockIdx.x;


	float(*pLocArry)[OutParaNumGS2D] = (float(*)[OutParaNumGS2D])d_LocArry; // for parameter array

	int HistPos = 0;


	if (gid < FluoNum)
	{
		float TotalPhoton = pLocArry[gid][Pos_TPho]; // 
		float Background = pLocArry[gid][Pos_Bakg];  // 

		float PSFSigma = (pLocArry[gid][Pos_SigX] + pLocArry[gid][Pos_SigY]) / 2; // pixel
		float CurSNR = pLocArry[gid][Pos_PSNR];  //

		float LocPrecX = pLocArry[gid][Pos_CrbX]; // nm
		float LocPrecY = pLocArry[gid][Pos_CrbY]; // nm

		float LocPrecXY = sqrtf(LocPrecX*LocPrecX + LocPrecY*LocPrecY);
			

		// get statistic information distribution and  avoid failure point
		if ((TotalPhoton > 30.0f) && (Background > 2.0f) && (CurSNR >= 3.0f) && (PSFSigma >= 0.01f) && (LocPrecXY > 0.01f) && (LocPrecXY < 55.0))
		{

			// total photon
			HistPos = TotalPhoton / Hist_PhotonGap; 
			if (HistPos < 0)HistPos = 0;
			if (HistPos >= StatInf_Hist_DatLenMax)HistPos = StatInf_Hist_DatLenMax - 1;

			atomicAdd(&d_Hist_TotalPhoton[HistPos], 1);

			// localization precision
			HistPos = LocPrecXY / Hist_LocPrecisionGap;
			if (HistPos < 0)HistPos = 0;
			if (HistPos >= StatInf_Hist_DatLenMax)HistPos = StatInf_Hist_DatLenMax - 1;

			atomicAdd(&d_Hist_LocPrecisionXY[HistPos], 1);


			//  peak SNR
			HistPos = CurSNR / Hist_PeakSNRGap;
			if (HistPos < 0)HistPos = 0;
			if (HistPos >= StatInf_Hist_DatLenMax)HistPos = StatInf_Hist_DatLenMax - 1;

			atomicAdd(&d_Hist_PeakSNR[HistPos], 1);


			// PSF width
			HistPos = PSFSigma / Hist_PSFWidthGap;
			if (HistPos < 0)HistPos = 0;
			if (HistPos >= StatInf_Hist_DatLenMax)HistPos = StatInf_Hist_DatLenMax - 1;

			atomicAdd(&d_Hist_PSFWidth[HistPos], 1);

			
			// get average SNR and PSF width
			atomicAdd(d_ValidNum, 1.0f);

			atomicAdd(d_LocPrecSum, LocPrecXY);
			atomicAdd(d_SNRSum, CurSNR);
			atomicAdd(d_PSFWSum, PSFSigma);
			atomicAdd(d_BgSum, Background);

		}
	}
}


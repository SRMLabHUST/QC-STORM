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

#include "bfgs_top.h"
#include "bfgs_LocalizationFilter.h"



// filter localizations with low quality or invalid information (content is 0)
__global__ void BFGSLocalizationFilter(float *d_LocArry, float MinSNR, float PixelSize, float QE, int ImageWidth, int ImageHigh, int FluoNum)
{
	int gid = threadIdx.x + blockDim.x*blockIdx.x;


	float(*pLocArry)[OutParaNumGS2D]; // for parameter array

	pLocArry = (float(*)[OutParaNumGS2D])d_LocArry;

	float PeakPhoton, TotalPhoton, Background;
	float XPos, YPos;

	float CurSNR, LocPrec;

	int CurPos = 0;

	int cnt = 0;

	int Valid0 = 0;
	int Valid1 = 0;

	if (gid < FluoNum)
	{	
		PeakPhoton = pLocArry[gid][Pos_PPho] * QE;
		XPos = pLocArry[gid][Pos_XPos];
		YPos = pLocArry[gid][Pos_YPos];
		TotalPhoton = pLocArry[gid][Pos_TPho] * QE; // convert to e-
		Background = pLocArry[gid][Pos_Bakg] * QE; // convert to e-
		CurSNR = pLocArry[gid][Pos_PSNR];
		LocPrec = pLocArry[gid][Pos_CrbX];

		if ((PeakPhoton < 1.0f) || (TotalPhoton < 1.0f))
		{
			PeakPhoton = 1.0f;
			TotalPhoton = 1.0f;
		}

		if (LocPrec >= 0)
		{
			// not nan
		}
		else
		{
			// could be nan
			pLocArry[gid][Pos_CrbX] = 10.0f;
		}


		Valid0 = (XPos >= 5) && (YPos >= 5) && (XPos < ImageWidth - 5) && (YPos < ImageHigh - 5);
		Valid1 = (CurSNR > MinSNR)&& (TotalPhoton > LocFilter_TotalPhoton_MinTh);


		if (!(Valid0 & Valid1))
		{
#pragma unroll
			for (cnt = 0; cnt < OutParaNumGS2D; cnt++)
			{
				// copy valid array
				pLocArry[gid][cnt] = 0.0f;
			}
		}
	}
}


// filter localizations with low quality or invalid information (content is 0)
__global__ void BFGSLocalizationMeanSNRCalc(float *d_LocArry, float *d_SNRSumUp, int *d_ValidNum, float SNRth,  int ImageWidth, int ImageHigh, int FluoNum)
{
	int gid = threadIdx.x + blockDim.x*blockIdx.x;


	float(*pLocArry)[OutParaNumGS2D]; // for parameter array

	pLocArry = (float(*)[OutParaNumGS2D])d_LocArry;

	float XPos, YPos;

	float CurSNR;


	if (gid < FluoNum)
	{
		XPos = pLocArry[gid][Pos_XPos];
		YPos = pLocArry[gid][Pos_YPos];
		CurSNR = pLocArry[gid][Pos_PSNR];


		if ((CurSNR >= SNRth) && (XPos >= 5) && (YPos >= 5) && (XPos < ImageWidth - 5) && (YPos < ImageHigh - 5))
		{
			atomicAdd(d_SNRSumUp, CurSNR);
			atomicAdd(d_ValidNum, 1);

		}
	}
}

// cuda wrapper

void LDLocData_TypeDef::FilterBadFit(LocalizationPara & LocPara, int FluoNum, cudaStream_t cstream)
{

	int BlockDim = ThreadsPerBlock;
	int BlockNum = ((FluoNum + ThreadsPerBlock - 1) / ThreadsPerBlock);

	float MinSNR = LocFilter_SNR_MinTh;

	if (LocPara.BadFitFilterWithAutoThEn)
	{
		// get mean SNR
		cudaMemsetAsync(d_SNRSumUp, 0, sizeof(float), cstream);
		cudaMemsetAsync(d_ValidNum, 0, sizeof(int), cstream);

		float SNRSelTh = 4.5f;

		BFGSLocalizationMeanSNRCalc << <BlockNum, BlockDim, 0, cstream >> >(d_LocArry, d_SNRSumUp, d_ValidNum, SNRSelTh, LocPara.ImageWidth, LocPara.ImageHigh, FluoNum);
		
		cudaMemcpyAsync(h_SNRSumUp, d_SNRSumUp, sizeof(float), cudaMemcpyDeviceToHost, cstream);
		cudaMemcpyAsync(h_ValidNum, d_ValidNum, sizeof(int), cudaMemcpyDeviceToHost, cstream);
		cudaStreamSynchronize(cstream);


		// automatically find SNR filtering threshold
		MinSNR = (*h_SNRSumUp) / (*h_ValidNum) / 2;
		MinSNR = max(MinSNR, LocFilter_SNR_MinTh);

	}

	MinSNR = min(MinSNR, LocFilter_SNR_MaxTh);

	//	printf("MinSNR:%f \n", MinSNR);

	BFGSLocalizationFilter << <BlockNum, BlockDim, 0, cstream >> >(d_LocArry, MinSNR, LocPara.PixelSize, LocPara.QE, LocPara.ImageWidth, LocPara.ImageHigh, FluoNum);


//	cudaMemcpyAsync(h_LocArry, d_LocArry, FluoNum*OutParaNumGS2D*sizeof(float), cudaMemcpyDeviceToHost, cstream);
//	cudaStreamSynchronize(cstream);

}

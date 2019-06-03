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

#include "DriftCorrection_SRImageXCorr.h"


void SRDriftCorrData_TypeDef::CorrectSampleShift(string FileName, string oFileName, LocalizationPara & LocPara, int CorrFrameNum, cudaStream_t cstream)
{
	// initial information array

	TotalFrame = GetTotalFrame(FileName);


	int cnt = 0;

	CalcGroupNum(TotalFrame, CorrFrameNum);


	GetCorrGroupFluoPos(FileName);

//	printf("group num:%d %d\n", TotalFrame, CorrGroupNum);

	/*
	for (cnt; cnt < CorrGroupNum; cnt++)
	{
	printf("frame pos: %d-%d, %d-%d\n", GroupStartFrame[cnt], GroupEndFrame[cnt], GroupFrameStartPos[cnt], GroupFrameEndPos[cnt]);
	}
	*/

	float ShiftX = 0;
	float ShiftY = 0;
	int CorrShiftBiasX = 0;
	int CorrShiftBiasY = 0;

	// rend  super-resolution image
	// 0th is the reference image and the others are drifted image


	for (cnt = 0; cnt < CorrGroupNum; cnt += 1)
	{
		RenderSlice(FileName, cnt, LocPara.PixelSize, LocPara.QE, LocPara.SNR_th, cstream);

		if (cnt > 0)
		{
			GetSliceShift( &ShiftX, &ShiftY, CorrShiftBiasX, CorrShiftBiasY, cstream);
			CorrShiftBiasX = (int)ShiftX;
			CorrShiftBiasY = (int)ShiftY;
		}

		XSliceShift[cnt] = ShiftX / SRShiftCorr_PixelZoom;
		YSliceShift[cnt] = ShiftY / SRShiftCorr_PixelZoom;

		printf("x y shift(pixel):%f %f\n", XSliceShift[cnt], YSliceShift[cnt]);

	}

	ShiftInterpolation();


	ApplyShiftTop(FileName, oFileName, BIT0 | BIT1, cstream);


	cudaStreamSynchronize(cstream);

}


void SRDriftCorrData_TypeDef::Init(int RawImgWidth, int RawImgHigh, int MaxFrameNum)
{

	cudaError_t err;

	// image information
	ImageWidth = RawImgWidth;
	ImageHigh = RawImgHigh;

	SRImageWidth = RawImgWidth*SRShiftCorr_PixelZoom;
	SRImageHigh = RawImgHigh*SRShiftCorr_PixelZoom;
		

	// create inf array for each group
	GroupStartFrame = new int[MaxCorrGroupNum]; // start frame of each group
	GroupEndFrame = new int[MaxCorrGroupNum]; // start frame of each group

	GroupFrameStartPos = new int[MaxCorrGroupNum]; // start fluo position of each group
	GroupFrameEndPos = new int[MaxCorrGroupNum]; // start fluo position of each group



	XSliceShift = new float[MaxCorrGroupNum];
	YSliceShift = new float[MaxCorrGroupNum];


	cudaMallocHost((void **)&XFrameShift, MaxFrameNum*sizeof(float));
	cudaMallocHost((void **)&YFrameShift, MaxFrameNum*sizeof(float));

	cudaMalloc((void **)&d_XFrameShift, MaxFrameNum*sizeof(float));
	cudaMalloc((void **)&d_YFrameShift, MaxFrameNum*sizeof(float));

	cudaMallocHost((void **)&h_LocArry, PointNumTh * 2 * OutParaNumGS2D*sizeof(float));
	cudaMalloc((void **)&d_LocArry, PointNumTh * 2 * OutParaNumGS2D*sizeof(float));

	err = cudaMallocHost((void **)&h_FillImg1, SRImageWidth*SRImageHigh*sizeof(float));
	HandleErr(err, "cudaMallocHost shift corr h_FillImg1");
	err = cudaMallocHost((void **)&h_FillImg2, SRImageWidth*SRImageHigh*sizeof(float));
	HandleErr(err, "cudaMallocHost shift corr h_FillImg2");
	err = cudaMallocHost((void **)&h_SumLine, SRImageWidth*sizeof(float));


	err = cudaMalloc((void **)&d_FillImg1, SRImageWidth*SRImageHigh*sizeof(float));
	HandleErr(err, "cudaMalloc shift corr d_FillImg1");
	err = cudaMalloc((void **)&d_FillImg2, SRImageWidth*SRImageHigh*sizeof(float));
	HandleErr(err, "cudaMalloc shift corr d_FillImg2");
	err = cudaMalloc((void **)&d_MulImg, SRImageWidth*SRImageHigh*sizeof(float));
	HandleErr(err, "cudaMalloc shift corr d_MulImg");
	cudaMalloc((void **)&d_SumLine, SRImageWidth*sizeof(float));


}


void SRDriftCorrData_TypeDef::Deinit()
{
	cudaError_t err;

	delete[] GroupStartFrame;
	delete[] GroupEndFrame;

	delete[] GroupFrameStartPos;
	delete[] GroupFrameEndPos;


	delete[] XSliceShift;
	delete[] YSliceShift;

	cudaFreeHost(XFrameShift);
	cudaFreeHost(YFrameShift);

	cudaFree(d_XFrameShift);
	cudaFree(d_YFrameShift);


	cudaFreeHost(h_LocArry);
	cudaFree(d_LocArry);

	err = cudaFreeHost(h_FillImg1);
	HandleErr(err, "cudaFreeHost shift corr h_FillImg1");
	err = cudaFreeHost(h_FillImg2);
	HandleErr(err, "cudaFreeHost shift corr h_FillImg2");
	err = cudaFreeHost(h_SumLine);

	err = cudaFree(d_FillImg1);
	HandleErr(err, "cudaFree shift corr d_FillImg1");
	err = cudaFree(d_FillImg2);
	HandleErr(err, "cudaFree shift corr d_FillImg2");
	err = cudaFree(d_MulImg);
	err = cudaFree(d_SumLine);

}


void SRDriftCorrData_TypeDef::ResetFillImage(float *d_SRIntensityImg, int SRImageWidth, int SRImageHigh, cudaStream_t cstream)
{

	cudaMemsetAsync(d_SRIntensityImg, 0, SRImageWidth*SRImageHigh*sizeof(float), cstream);
	cudaStreamSynchronize(cstream);

}


void SRDriftCorrData_TypeDef::GetSliceShift(float *ShiftX, float *ShiftY, int CorrShiftBiasX, int CorrShiftBiasY, cudaStream_t cstream)
{
	// calculate cross correletion
	// don't calculate all, only calculate a useful region

	int CorrShiftX = 0;
	int CorrShiftY = 0;


	//	printf("cor size:%d %d\n", SRImageWidth, SRImageHigh);

	int cnt = 0;
	double CorrResult[CorrSize][CorrSize];

	int xcnt, ycnt;

	double MaxSumDat = 0;
	int MaxPosX = 0, MaxPosY = 0;


	float SumCorrX[2 * FittingRadius + 1]; // sum along x,y direction for fitting region
	float SumCorrY[2 * FittingRadius + 1]; // sum along x,y direction for fitting region

	int FitXS, FitXE, FitYS, FitYE;// X,Y start,end pos


	// cross-correlation calculation for a 51x51 region
	for (ycnt = 0; ycnt < CorrSize; ycnt++)
	{
		for (xcnt = 0; xcnt < CorrSize; xcnt++)
		{
			CorrShiftX = xcnt - HalfCorrSize;
			CorrShiftY = ycnt - HalfCorrSize;

			CorrResult[ycnt][xcnt] = CrossCorrelation(CorrShiftX, CorrShiftY, CorrShiftBiasX, CorrShiftBiasY, cstream);
			if (MaxSumDat < CorrResult[ycnt][xcnt])
			{
				MaxSumDat = CorrResult[ycnt][xcnt];
				MaxPosX = xcnt;
				MaxPosY = ycnt;
			}
			//			printf("%.1f ", CorrResult[ycnt][xcnt]);
		}
	}
	// normalize the correlation result
	for (ycnt = 0; ycnt < CorrSize; ycnt++)
	{
		for (xcnt = 0; xcnt < CorrSize; xcnt++)
		{

			CorrResult[ycnt][xcnt] = CorrResult[ycnt][xcnt] / MaxSumDat * 1000;
			//			printf("%.1f ", CorrResult[ycnt][xcnt]);
		}
	}


	//	printf("\nMaxPosX:%d %d\n", MaxPosX, MaxPosY);

	// find the fitting region centered with the max pos
	FitXS = max(MaxPosX - FittingRadius, 0);
	FitYS = max(MaxPosY - FittingRadius, 0);
	FitXE = min(MaxPosX + FittingRadius, CorrSize - 1);
	FitYE = min(MaxPosY + FittingRadius, CorrSize - 1);

	int xhlen = (FitXE - FitXS + 1) / 2;
	int yhlen = (FitYE - FitYS + 1) / 2;

	for (cnt = 0; cnt < 2 * FittingRadius + 1; cnt++)
	{
		SumCorrX[cnt] = 0;
		SumCorrY[cnt] = 0;
	}
	// center of mass fitting of gaussian shape cross-correlation
	for (ycnt = FitYS; ycnt <= FitYE; ycnt++)
	{
		for (xcnt = FitXS; xcnt <= FitXE; xcnt++)
		{
			SumCorrX[FittingRadius - xhlen + (xcnt - FitXS)] += CorrResult[ycnt][xcnt];
			SumCorrY[FittingRadius - yhlen + (ycnt - FitYS)] += CorrResult[ycnt][xcnt];
		}
	}

	float wSum = 0, cSum = 0; // weighted sum for center of mass
	float CenterX, CenterY;

	for (cnt = 0; cnt < 2 * FittingRadius + 1; cnt++)
	{
		if (SumCorrX[cnt] == 0)continue;

		wSum += SumCorrX[cnt] * (cnt + 1);
		cSum += SumCorrX[cnt];

	}
	CenterX = wSum / cSum - (FittingRadius + 1);

	wSum = 0;
	cSum = 0;

	for (cnt = 0; cnt < 2 * FittingRadius + 1; cnt++)
	{
		if (SumCorrY[cnt] == 0)continue;

		wSum += SumCorrY[cnt] * (cnt + 1);
		cSum += SumCorrY[cnt];

	}
	CenterY = wSum / cSum - (FittingRadius + 1);


	*ShiftX = CorrShiftBiasX + MaxPosX - HalfCorrSize + CenterX;
	*ShiftY = CorrShiftBiasY + MaxPosY - HalfCorrSize + CenterY;

	//	printf("shift pos:%f %f\n", *ShiftX, *ShiftY);

}


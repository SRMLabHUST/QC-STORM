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




void SRDriftCorrData_TypeDef::ShiftInterpolation()
{
	// set shift to 0
	memset(XFrameShift, 0, TotalFrame*sizeof(float));
	memset(YFrameShift, 0, TotalFrame*sizeof(float));

	int TotalFrame = TotalFrame;


	int CurPos = 0;

	// from 0 frame to the first center
	float LastCenterFrame = 0;
	float LastShiftX = 0;
	float LastShiftY = 0;

	float CurCenterFrame = 0;
	float CurShiftX = 0;
	float CurShiftY = 0;


	int cnt = 0;
	int fcnt = 0;
	float KX = 0, KY = 0;

	for (cnt = 0; cnt <= CorrGroupNum; cnt++)
	{
		if (cnt == CorrGroupNum)
		{
			// end group
			CurCenterFrame = GroupEndFrame[CorrGroupNum - 1]; // the end frame
			// KX KY keep the same with the last group
		}
		else
		{
			// normal group include 0 group
			CurCenterFrame = (GroupStartFrame[cnt] + GroupEndFrame[cnt]) / 2;
			CurShiftX = XSliceShift[cnt];
			CurShiftY = YSliceShift[cnt];

			KX = (CurShiftX - LastShiftX) / (CurCenterFrame - LastCenterFrame);
			KY = (CurShiftY - LastShiftY) / (CurCenterFrame - LastCenterFrame);
		}
		
		for (fcnt = LastCenterFrame; fcnt < CurCenterFrame; fcnt++)
		{
			XFrameShift[fcnt] = LastShiftX + (fcnt - LastCenterFrame) *KX;
			YFrameShift[fcnt] = LastShiftY + (fcnt - LastCenterFrame) *KY;
		}

		LastCenterFrame = CurCenterFrame;
		LastShiftX = CurShiftX;
		LastShiftY = CurShiftY;

	}
}

void SRDriftCorrData_TypeDef::ApplyShiftTop(string iFileName, string oFileName, int ShiftCorrEnable, cudaStream_t cstream)
{
	CppBinaryReadFile iLocFile(iFileName);

	CppBinaryWriteFile wfile(oFileName);



	cudaMemcpyAsync(d_XFrameShift, XFrameShift, TotalFrame*sizeof(float), cudaMemcpyHostToDevice, cstream);
	cudaMemcpyAsync(d_YFrameShift, YFrameShift, TotalFrame*sizeof(float), cudaMemcpyHostToDevice, cstream);
	cudaStreamSynchronize(cstream);


	int TotalFluoNum = iLocFile.GetLength() / sizeof(float) / OutParaNumGS2D;

	int ProcNum = TotalFluoNum / PointNumTh + 1;
	int CurFluoNum = PointNumTh;

	int cnt;

	// find image size first
	for (cnt = 0; cnt < ProcNum; cnt++)
	{
		if (cnt == ProcNum - 1)CurFluoNum = TotalFluoNum - cnt*PointNumTh;
		else CurFluoNum = PointNumTh;

		iLocFile.Read(h_LocArry, CurFluoNum*OutParaNumGS2D*sizeof(float));


		gpuApplyShift(h_LocArry, ShiftCorrEnable, CurFluoNum, cstream);

		wfile.Write(h_LocArry, CurFluoNum*OutParaNumGS2D*sizeof(float));
	}

	wfile.Close();
	iLocFile.Close();
}


void SRDriftCorrData_TypeDef::RenderSlice(string FileName, int RendGroup, float PixelSize, float QE, float SNR_th, cudaStream_t cstream)
{
	int FirstFluoPos, EndFluoPos;
	float *d_FillImage = NULL;

	FirstFluoPos = GroupFrameStartPos[RendGroup];
	EndFluoPos = GroupFrameEndPos[RendGroup];

	// render gap
	int TotalFluoNum = EndFluoPos - FirstFluoPos + 1;

//	printf("cur rend:%d-%d\n", FirstFluoPos, EndFluoPos);

	if (RendGroup == 0)
	{
		d_FillImage = d_FillImg1;
	}
	else
	{
		d_FillImage = d_FillImg2;
	}

	// reset this image before rendering
	ResetFillImage(d_FillImage, SRImageWidth, SRImageHigh, cstream);


	CppBinaryReadFile iLocFile(FileName);



	iLocFile.Seek(FirstFluoPos*OutParaNumGS2D*sizeof(float), ios::beg);

	int ProcNum = (TotalFluoNum + PointNumTh - 1) / PointNumTh;
	int CurFluoNum = PointNumTh;

	int cnt;

	// find image size first
	for (cnt = 0; cnt < ProcNum; cnt++)
	{
		if (cnt == ProcNum - 1)CurFluoNum = TotalFluoNum - cnt*PointNumTh;
		else CurFluoNum = PointNumTh;

		iLocFile.Read(h_LocArry, CurFluoNum*OutParaNumGS2D*sizeof(float));
	
		ImageRender(h_LocArry, d_LocArry, d_FillImage, QE, SNR_th, PixelSize, SRShiftCorr_PixelZoom, SRImageWidth, SRImageHigh, CurFluoNum, cstream);

	}

	iLocFile.Close();
}


void SRDriftCorrData_TypeDef::CalcGroupNum( int TotalFrame, int CorrFrameNum)
{
	// note frame is 1 to total frame
	int HalfCorrFrameNum = CorrFrameNum / 2;

	int GroupNum = TotalFrame / CorrFrameNum;

	GroupNum = GroupNum * 2 - 1; // shift with 1/2 overlaping group frames

	int ResidualFrameNum = TotalFrame % CorrFrameNum;

	int cnt;


	for (cnt = 0; cnt < GroupNum; cnt++)
	{
		// frame assign for each group
		GroupStartFrame[cnt] = cnt*HalfCorrFrameNum + 1; // start frame of each group
		GroupEndFrame[cnt] = cnt*HalfCorrFrameNum + CorrFrameNum; // start frame of each group
	}

	// if there are some residual frames
	if (ResidualFrameNum < HalfCorrFrameNum)
	{
		// residual grouped into last group
		GroupEndFrame[GroupNum - 1] = TotalFrame; // start frame of each group

	}
	else
	{
		// residual is a new group
		GroupNum++;
		GroupStartFrame[GroupNum - 1] = TotalFrame - CorrFrameNum + 1; // start frame of each group
		GroupEndFrame[GroupNum - 1] = TotalFrame; // start frame of each group
	}

	CorrGroupNum = GroupNum;
}


void SRDriftCorrData_TypeDef::GetCorrGroupFluoPos(string FileName)
{
	int cnt = 0;
	
	int BeginPos = 0;
	int EndPos = 0;
	GroupFrameStartPos[0] = 0;

	for (cnt = 0; cnt < CorrGroupNum; cnt++)
	{

		// end fluo position of each group
		EndPos = GetAFrameEndFluoPos(FileName, EndPos, GroupEndFrame[cnt]);

		GroupFrameEndPos[cnt] = EndPos; 

		if (cnt >= 1)
		{
			// start fluo position of each group
			BeginPos = GetAFrameEndFluoPos(FileName, BeginPos, GroupStartFrame[cnt] - 1);

			GroupFrameStartPos[cnt] = BeginPos + 1;
		}
	}
}

int SRDriftCorrData_TypeDef::GetAFrameEndFluoPos(string FileName, int OffsetFluoNum, int FindFrame)
{
	// find end molecular position for a frame
	int EndPos = 0;


	CppBinaryReadFile iLocFile(FileName);


	int TotalFluoNum = iLocFile.GetLength() / sizeof(float) / OutParaNumGS2D;
	int ResidFluoNum = TotalFluoNum - OffsetFluoNum;


	int ProcNum = (ResidFluoNum + PointNumTh - 1) / PointNumTh;
	int pcnt;

	//	printf("fluonum:%d %d %d\n", TotalFluoNum, ResidFluoNum, ProcNum);

	float *PointData = new float[PointNumTh*OutParaNumGS2D];

	iLocFile.Seek(OffsetFluoNum*OutParaNumGS2D*sizeof(float), ios::beg);


	int IsBreak = 0;
	int curFluoNum = PointNumTh;
	for (pcnt = 0; pcnt < ProcNum; pcnt++)
	{
		if (pcnt == ProcNum - 1)curFluoNum = ResidFluoNum - pcnt*PointNumTh;
		else curFluoNum = PointNumTh;

		iLocFile.Read(PointData, curFluoNum*OutParaNumGS2D*sizeof(float));

		int cnt = 0;
		int CurFrame = 0;
		float(*pLocArry)[OutParaNumGS2D] = (float(*)[OutParaNumGS2D])PointData;

		for (cnt = 0; cnt < curFluoNum; cnt++)
		{
			CurFrame = pLocArry[cnt][OutParaNumGS2D - 1];

			if ((CurFrame > 0) && (CurFrame > FindFrame))
			{
				EndPos = OffsetFluoNum + pcnt*PointNumTh + cnt - 1;
				IsBreak = 1;
				break;
			}

		}
		if (IsBreak)break;
	}

	iLocFile.Close();
	delete[] PointData;

	if (EndPos == 0)EndPos = TotalFluoNum - 1; // can't find until end

	return EndPos;
}


// code wrapper

void SRDriftCorrData_TypeDef::ImageRender(float *h_LocArry, float *d_LocArry, float *d_SRIntensityImg, float QE, float SNR_th, float PixelSize, float PixelZoom, int SRImageWidth, int SRImageHigh, int FluoNum, cudaStream_t cstream)
{
	cudaMemcpyAsync(d_LocArry, h_LocArry, FluoNum*OutParaNumGS2D*sizeof(float), cudaMemcpyHostToDevice, cstream);

	int BlockDim = ThreadsPerBlock;
	int BlockNum = ((FluoNum + ThreadsPerBlock - 1) / ThreadsPerBlock);

	PointFillRegionWithPrec1 << <BlockNum, BlockDim, 0, cstream >> >(d_LocArry, d_SRIntensityImg, QE, SNR_th, PixelSize, PixelZoom, SRImageWidth, SRImageHigh, FluoNum);

	cudaStreamSynchronize(cstream);
}

void SRDriftCorrData_TypeDef::gpuApplyShift(float *h_LocArry, int ShiftCorrEnable, int FluoNum, cudaStream_t cstream)
{
	int BlockDim = ThreadsPerBlock;
	int BlockNum = ((FluoNum + ThreadsPerBlock - 1) / ThreadsPerBlock);

	cudaMemcpyAsync(d_LocArry, h_LocArry, FluoNum*OutParaNumGS2D*sizeof(float), cudaMemcpyHostToDevice, cstream);



	gpuApplyShiftCorr << <BlockNum, BlockDim, 0, cstream >> >(d_LocArry, d_XFrameShift, d_YFrameShift, NULL, ShiftCorrEnable, FluoNum);

	cudaMemcpyAsync(h_LocArry, d_LocArry, FluoNum*OutParaNumGS2D*sizeof(float), cudaMemcpyDeviceToHost, cstream);

	cudaStreamSynchronize(cstream);
}

double SRDriftCorrData_TypeDef::CrossCorrelation(int ShiftX, int ShiftY, int CorrShiftBiasX, int CorrShiftBiasY, cudaStream_t cstream)
{
	double TotalSum = 0;
	int TotalPixel = SRImageWidth*SRImageHigh;

	cudaMemsetAsync(d_MulImg, 0, SRImageWidth*SRImageHigh*sizeof(float), cstream);


	int BlockDim = ThreadsPerBlock;
	int BlockNum = ((TotalPixel + ThreadsPerBlock - 1) / ThreadsPerBlock);

	// calculate multiply of two image
	gpuCrossCorr_Mul << <BlockNum, BlockDim, 0, cstream >> >(d_FillImg1, d_FillImg2, d_MulImg, ShiftX, ShiftY, CorrShiftBiasX, CorrShiftBiasY, SRImageWidth, SRImageHigh);


	TotalPixel = SRImageWidth;
	BlockNum = ((TotalPixel + ThreadsPerBlock - 1) / ThreadsPerBlock);
	// calculate the sum of the multipl
	gpuCrossCorr_SumOfMul << <BlockNum, BlockDim, 0, cstream >> >(d_MulImg, d_SumLine, SRImageWidth, SRImageHigh);

	cudaMemcpyAsync(h_SumLine, d_SumLine, SRImageWidth*sizeof(float), cudaMemcpyDeviceToHost, cstream);

	cudaStreamSynchronize(cstream);

	int cnt = 0;
	TotalSum = 0;
	for (cnt = 0; cnt < SRImageWidth; cnt++)
	{
		TotalSum += h_SumLine[cnt];
	}

	return TotalSum;
}



int SRDriftCorrData_TypeDef::GetTotalFrame(string FileName)
{
	// find total frame number
	int TotalFrame = 0;


	CppBinaryReadFile iLocFile(FileName);


	int FluoNum = iLocFile.GetLength() / OutParaNumGS2D / sizeof(float);
	int SelDataNum = min(FluoNum, 400);


	float *PointData = new float[SelDataNum * OutParaNumGS2D];


	iLocFile.Seek(-SelDataNum * OutParaNumGS2D*sizeof(float), ios::end);

	iLocFile.Read(PointData, SelDataNum * OutParaNumGS2D*sizeof(float));

	float(*pLocArry)[OutParaNumGS2D] = (float(*)[OutParaNumGS2D])PointData;
	int CurFrame;

	int cnt = 0;
	for (cnt = SelDataNum-1; cnt >= 0; cnt--)
	{
		CurFrame = pLocArry[cnt][Pos_Frme];
		if (CurFrame > 0)
		{
			TotalFrame = CurFrame;
			break;
		}
	}

	delete[]PointData;
	iLocFile.Close();

	return TotalFrame;
}

void SRDriftCorrData_TypeDef::GetMaxImgSize(string FileName, int *ImageWidth, int *ImageHigh)
{

	float *h_LocArry;
	float *d_LocArry;
	cudaStream_t cstream;

	cudaMallocHost((void **)&h_LocArry, PointNumTh * OutParaNumGS2D*sizeof(float));
	cudaMalloc((void **)&d_LocArry, PointNumTh * OutParaNumGS2D*sizeof(float));
	cudaStreamCreate(&cstream);


	CppBinaryReadFile iLocFile(FileName);


	int TotalFluoNum = iLocFile.GetLength() / sizeof(float) / OutParaNumGS2D;


	int FindImgWidth = 0, FindImgHigh = 0;
	*ImageWidth = 0;
	*ImageHigh = 0;

	const int ProcNum = (TotalFluoNum + PointNumTh - 1) / PointNumTh;

	int CurFluoNum;
	int cnt;

	// find image size first
	for (cnt = 0; cnt < ProcNum; cnt++)
	{
		if (cnt == ProcNum - 1)
		{
			CurFluoNum = TotalFluoNum % PointNumTh;
		}
		else
		{
			CurFluoNum = PointNumTh;
		}

		iLocFile.Read(h_LocArry, CurFluoNum*OutParaNumGS2D*sizeof(float));

		ImageRenderData_TypeDef::GetMaxImgSizeFromLocArry(h_LocArry, d_LocArry, &FindImgWidth, &FindImgHigh, CurFluoNum, cstream);

		*ImageWidth = max(*ImageWidth, FindImgWidth);
		*ImageHigh = max(*ImageHigh, FindImgHigh);

		if (cnt > 10)break;
	}

	*ImageWidth = (*ImageWidth + 8 - 1) / 8 * 8;
	*ImageHigh = (*ImageHigh + 8 - 1) / 8 * 8;


	cudaFreeHost(h_LocArry);
	cudaFree(d_LocArry);
	cudaStreamDestroy(cstream);

	iLocFile.Close();
}


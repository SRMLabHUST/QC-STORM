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


#include "LDROIExtraction.h"

#define JudgeROISize		3
#define JudgeROISize_Half				((int)(JudgeROISize/2))




__global__ void gpuROIFindingLD(unsigned short *d_RawImg_Smoothed, unsigned short * d_MoleculePosImage, int *d_ValidROINum, int *d_ROIMarkInfArray, float ImageVariance, int LocType, int MultiEmitterFitEn, int ROISize, int ImageWidth, int BatchedImageHigh, int ImageHigh)
{

	int gid = threadIdx.x + blockDim.x*blockIdx.x;
	int tid = threadIdx.x;


	// standard deviation
	ImageVariance = sqrtf(ImageVariance);

	float ImageTh = 9.5f * ImageVariance;
	float ImageTh1 = 3.8f * ImageVariance;

	if ((LocType == LocType_AS3D) || ((MultiEmitterFitEn == 0) && (LocType == LocType_GS2D)))
	{
		ImageTh = 8.5f * ImageVariance;
		ImageTh1 = 3.5f * ImageVariance;
	}

	// position in batched position
	int CurXPos = gid % ImageWidth;
	int CurYPos = gid / ImageWidth;

	// position in each image
	int RealXPos = CurXPos;
	int RealYPos = CurYPos % ImageHigh;


	const int BatchedImageNum = BatchedImageHigh / ImageHigh;


	int ImageROI[JudgeROISize*JudgeROISize];

	int(*pImageRegion)[JudgeROISize] = (int(*)[JudgeROISize])ImageROI;

	//
	int AddrOffset = 0;
	int ROISize_Half = ((int)(ROISize / 2));

	if ((RealXPos >= ROISize_Half) && (RealYPos >= ROISize_Half) && (RealXPos < ImageWidth - ROISize_Half) && (RealYPos < ImageHigh - ROISize_Half))
	{

		int CenterPixelVal = d_RawImg_Smoothed[CurYPos * ImageWidth + CurXPos];


		if (CenterPixelVal >= ImageTh)
		{

#pragma unroll
			for (int r = 0; r < JudgeROISize; r++)
			{
				AddrOffset = (CurYPos - JudgeROISize_Half + r)*ImageWidth + (CurXPos - JudgeROISize_Half);

#pragma unroll
				for (int c = 0; c < JudgeROISize; c++)
				{
					pImageRegion[r][c] = d_RawImg_Smoothed[AddrOffset + c];

				}
			}

#define JudgeBitNum_n8				8
#define JudgeBitNum_n4				4

			int JudgeBit1[JudgeBitNum_n8]; // 8 neighbor

			int JudgeBit_i4[JudgeBitNum_n4]; // judgement inner 4 neighbor
			int JudgeBit_o4[JudgeBitNum_n4]; // judgement outer 4 neighbor

			// subregion judgement
			JudgeBit1[0] = pImageRegion[0][0] >  ImageTh1;
			JudgeBit1[1] = pImageRegion[0][1] >  ImageTh1;
			JudgeBit1[2] = pImageRegion[0][2] >  ImageTh1;

			JudgeBit1[3] = pImageRegion[1][0] >  ImageTh1;
			JudgeBit1[4] = pImageRegion[1][2] >  ImageTh1;

			JudgeBit1[5] = pImageRegion[2][0] >  ImageTh1;
			JudgeBit1[6] = pImageRegion[2][1] >  ImageTh1;
			JudgeBit1[7] = pImageRegion[2][2] >  ImageTh1;

			// subregion judgement
			// inner 4 neighbor
			JudgeBit_i4[0]= CenterPixelVal >  pImageRegion[0][1];
			JudgeBit_i4[1]= CenterPixelVal >= pImageRegion[1][2];
			JudgeBit_i4[2]= CenterPixelVal >= pImageRegion[2][1];
			JudgeBit_i4[3]= CenterPixelVal >  pImageRegion[1][0];

			// outer 4 neighbor
			JudgeBit_o4[0] = CenterPixelVal > pImageRegion[0][0];
			JudgeBit_o4[1] = CenterPixelVal > pImageRegion[0][2];
			JudgeBit_o4[2] = CenterPixelVal > pImageRegion[2][0];
			JudgeBit_o4[3] = CenterPixelVal > pImageRegion[2][2];


			int SumJudge1 = 0;
			int SumJudge_i = 0;
			int SumJudge_o = 0;


			for (int jcnt = 0; jcnt < JudgeBitNum_n8; jcnt++)
			{
				SumJudge1 += JudgeBit1[jcnt];
			}

			for (int jcnt = 0; jcnt < JudgeBitNum_n4; jcnt++)
			{
				SumJudge_i += JudgeBit_i4[jcnt];
				SumJudge_o += JudgeBit_o4[jcnt];
			}

			int ROIType = 0;

			if ((SumJudge1 >= JudgeBitNum_n8 - 1) && (SumJudge_i == 4) && (SumJudge_o >= 3))
			{

				// SumJudge_o == 3 : potential
				// SumJudge_o == 4 : valid

				// mark a molecule in the image
				d_MoleculePosImage[gid] = 1;


				if (SumJudge_o == 4)
				{

					int CurROIInfPos = atomicAdd(d_ValidROINum, 1);

					if (CurROIInfPos < MaxPointNum)
					{
						d_ROIMarkInfArray[CurROIInfPos * ROIMarkInfNum + ROIMarkInf_XPos] = CurXPos;
						d_ROIMarkInfArray[CurROIInfPos * ROIMarkInfNum + ROIMarkInf_YPos] = CurYPos;

						//				printf("%d %d %d \n", CurXPos, CurYPos, ROIType);
					}
				}
			}
		}
	}
}


void ROIFindingLD(unsigned short *d_RawImg_Smoothed, unsigned short * d_MoleculePosImage, int *d_ValidROINum, int *d_ROIMarkInfArray,  float ImageVariance, int LocType, int MultiEmitterFitEn, int ROISize, int ImageWidth, int BatchedImageHigh, int ImageHigh, cudaStream_t cstream)
{
	int TotalThreadNum = ImageWidth * BatchedImageHigh;

	int BlockDim = ThreadsPerBlock;
	int BlockNum = (TotalThreadNum + ThreadsPerBlock - 1) / ThreadsPerBlock;


	gpuROIFindingLD << < BlockNum, BlockDim, 0, cstream >> > (d_RawImg_Smoothed, d_MoleculePosImage, d_ValidROINum, d_ROIMarkInfArray, ImageVariance, LocType, MultiEmitterFitEn, ROISize, ImageWidth, BatchedImageHigh, ImageHigh);

}


// note the ROI clasification are processed in gpuDetectedROIClasify and WLEParaEstimation
__global__ void gpuDetectedROIClasify(unsigned short * d_MoleculePosImage, int *d_ROIMarkInfArray, int ROINumber, int MultiEmitterFitEn, int ROISize, int ImageWidth)
{
	int gid = threadIdx.x + blockDim.x*blockIdx.x;

	int(*pROIMarkInfArray)[ROIMarkInfNum] = (int(*)[ROIMarkInfNum])d_ROIMarkInfArray;

	int NeighborFind_ROISize = ROISize;
	int NeighborFind_ROISize_Half = NeighborFind_ROISize / 2;
	 
	if (gid < ROINumber)
	{
		int CurX = pROIMarkInfArray[gid][ROIMarkInf_XPos];
		int CurY = pROIMarkInfArray[gid][ROIMarkInf_YPos];


		int NeighborNum = 0;
		int NeighborNum_i3 = 0;

		if (MultiEmitterFitEn)
		{
			// neighbor number in NeighborFind_ROISize
			for (int r = 0; r < NeighborFind_ROISize; r++)
			{
				for (int c = 0; c < NeighborFind_ROISize; c++)
				{
					int ReadY = CurY + r - NeighborFind_ROISize_Half;
					int ReadX = CurX + c - NeighborFind_ROISize_Half;
					int ReadPos = ReadY*ImageWidth + ReadX;

					NeighborNum += d_MoleculePosImage[ReadPos];
				}
			}

			// neighbor number in 3x3 ROI
			for (int r = 0; r < 3; r++)
			{
				for (int c = 0; c < 3; c++)
				{
					int ReadY = CurY + r - 1;
					int ReadX = CurX + c - 1;
					int ReadPos = ReadY*ImageWidth + ReadX;

					NeighborNum_i3 += d_MoleculePosImage[ReadPos];
				}
			}

			NeighborNum -= NeighborNum_i3;

			if (NeighborNum > 0)
			{
				pROIMarkInfArray[gid][ROIMarkInf_TypeFit] = ROIType_Fit_Multi;
			}
			else
			{
				pROIMarkInfArray[gid][ROIMarkInf_TypeFit] = ROIType_Fit_Single;
			}
		}
		else
		{
			pROIMarkInfArray[gid][ROIMarkInf_TypeFit] = ROIType_Fit_Single;
		}
	}
}


void DetectedROIClasify(unsigned short * d_MoleculePosImage, int *d_ROIMarkInfArray, int ROINumber, int MultiEmitterFitEn, int ROISize, int ImageWidth, cudaStream_t cstream)
{
	int BlockDim = ThreadsPerBlock;
	int BlockNum = (ROINumber + ThreadsPerBlock - 1) / ThreadsPerBlock;

	gpuDetectedROIClasify << < BlockNum, BlockDim, 0, cstream >> > (d_MoleculePosImage, d_ROIMarkInfArray, ROINumber, MultiEmitterFitEn, ROISize, ImageWidth);

}



__global__ void gpuSubregionExtractionLD(unsigned short *d_RawImg, unsigned short *d_ImageROI, int *d_ROIMarkInfArray, int ROISize, int ROINumber, int ImageWidth, int BatchedImageHigh, int ImageHigh, int StartFrame)
{
	int gid = threadIdx.x + blockDim.x*blockIdx.x;

	int(*pROIMarkInfArray_V)[ROIMarkInfNum] = (int(*)[ROIMarkInfNum])d_ROIMarkInfArray;

	if (gid < ROINumber)
	{
		int ROIWholeSize = ROISize*(ROISize + 1);
		int ROISize_Half = ROISize / 2;


		int HaflRegionSize = ROISize / 2;

		int RegionAddrOffset = gid*ROIWholeSize;

		int XPos = pROIMarkInfArray_V[gid][ROIMarkInf_XPos];
		int YPos = pROIMarkInfArray_V[gid][ROIMarkInf_YPos];

		int FitType = pROIMarkInfArray_V[gid][ROIMarkInf_TypeFit];

		if ((XPos >= ROISize_Half) && (YPos >= ROISize_Half) && (XPos < ImageWidth - ROISize_Half) && (YPos < BatchedImageHigh - ROISize_Half))
		{
			for (int ycnt = 0; ycnt < ROISize; ycnt++)
			{
				for (int xcnt = 0; xcnt < ROISize; xcnt++)
				{
					int XPos_GMem = XPos - ROISize_Half + xcnt;
					int YPos_GMem = YPos - ROISize_Half + ycnt;

					d_ImageROI[RegionAddrOffset + ycnt*ROISize + xcnt] = d_RawImg[YPos_GMem*ImageWidth + XPos_GMem];
				}
			}

			int AddrOffset = ROISize*ROISize;
			int CurFrame = StartFrame + (YPos / ImageHigh);

			d_ImageROI[RegionAddrOffset + AddrOffset + 0] = XPos;
			d_ImageROI[RegionAddrOffset + AddrOffset + 1] = (YPos % ImageHigh);
			d_ImageROI[RegionAddrOffset + AddrOffset + 2] = CurFrame % 65536; //
			d_ImageROI[RegionAddrOffset + AddrOffset + 3] = CurFrame / 65536;
			d_ImageROI[RegionAddrOffset + AddrOffset + 4] = FitType;

		}
	}
}



void SubregionExtractionLD(unsigned short *d_RawImg, unsigned short *d_ImageROI, int *d_ROIMarkInfArray, int ROISize, int ROINumber, int ImageWidth, int BatchedImageHigh, int ImageHigh, int StartFrame, cudaStream_t cstream)
{
	int BlockDim = ThreadsPerBlock;
	int BlockNum = (ROINumber + ThreadsPerBlock - 1) / ThreadsPerBlock;


	// region extraction with point position
	gpuSubregionExtractionLD << < BlockNum, BlockDim, 0, cstream >> > (d_RawImg, d_ImageROI, d_ROIMarkInfArray, ROISize, ROINumber, ImageWidth, BatchedImageHigh, ImageHigh, StartFrame);


}


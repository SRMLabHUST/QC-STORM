#include "LDROIExtraction.h"

#define JudgeROISize		3
#define JudgeROISize_Half				((int)(JudgeROISize/2))



__global__ void gpuROIFindingLD(unsigned short *d_RawImg_Smoothed, int *d_ROIPosArray, int *d_ROINumPerImage, float ImageVariance, int ROISize, int ImageWidth, int BatchedImageHigh, int ImageHigh);

__global__ void gpuSubregionExtractionLD(unsigned short *d_RawImg, unsigned short *d_ROIMem, int *d_ROIPosArray, int ROISize, int RegionNum_CurImage, int ImageWidth, int BatchedImageHigh, int ImageHigh, int StartFrame);



void ROIFindingLD(unsigned short *d_RawImg_Smoothed, int *d_ROIPosArray, int *d_ROINumPerImage, float ImageVariance, int ROISize, int ImageWidth, int BatchedImageHigh, int ImageHigh, cudaStream_t cstream)
{
	int TotalThreadNum = ImageWidth * BatchedImageHigh;

	int BlockDim = ThreadsPerBlock;
	int BlockNum = (TotalThreadNum + ThreadsPerBlock - 1) / ThreadsPerBlock;

	gpuROIFindingLD << < BlockNum, BlockDim, 0, cstream >> > (d_RawImg_Smoothed, d_ROIPosArray, d_ROINumPerImage, ImageVariance, ROISize, ImageWidth, BatchedImageHigh, ImageHigh);

}


__global__ void gpuROIFindingLD(unsigned short *d_RawImg_Smoothed, int *d_ROIPosArray, int *d_ROINumPerImage, float ImageVariance, int ROISize, int ImageWidth, int BatchedImageHigh, int ImageHigh)
{

	int gid = threadIdx.x + blockDim.x*blockIdx.x;
	int tid = threadIdx.x;


	// standard deviation
	ImageVariance = sqrtf(ImageVariance);

	float ImageTh = 10.0f * ImageVariance;
	float ImageTh1 = 7.0f * ImageVariance;


	int CurXPos = gid % ImageWidth;
	int CurYPos = gid / ImageWidth;

	int RealXPos = CurXPos;
	int RealYPos = CurYPos%ImageHigh;


	const int BatchedImageNum = BatchedImageHigh / ImageHigh;
	const int MaxFluoNumPerImage = GetMaxFluoNumPerImage(BatchedImageNum);


	int ImageRegion[JudgeROISize*JudgeROISize];

	int(*pImageRegion)[JudgeROISize] = (int(*)[JudgeROISize])ImageRegion;

	//
	int AddrOffset = 0;
	int ROISize_Half = ((int)(ROISize / 2));


	if ((RealXPos > ROISize_Half) && (RealYPos > ROISize_Half) && (RealXPos < ImageWidth - ROISize_Half) && (RealYPos < ImageHigh - ROISize_Half))
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

#define JudgeBitNum				8

			int JudgeBit[JudgeBitNum]; // judgement memory
			int JudgeBit1[JudgeBitNum]; // judgement memory

			// subregion judgement
			JudgeBit[0] = CenterPixelVal >  pImageRegion[0][0];
			JudgeBit[1] = CenterPixelVal >  pImageRegion[0][1];
			JudgeBit[2] = CenterPixelVal >  pImageRegion[0][2];

			JudgeBit[3] = CenterPixelVal >  pImageRegion[1][0];
			JudgeBit[4] = CenterPixelVal >= pImageRegion[1][2];

			JudgeBit[5] = CenterPixelVal >  pImageRegion[2][0];
			JudgeBit[6] = CenterPixelVal >= pImageRegion[2][1];
			JudgeBit[7] = CenterPixelVal >  pImageRegion[2][2];

			// subregion judgement
			JudgeBit1[0] = pImageRegion[0][0] >  ImageTh1;
			JudgeBit1[1] = pImageRegion[0][1] >  ImageTh1;
			JudgeBit1[2] = pImageRegion[0][2] >  ImageTh1;

			JudgeBit1[3] = pImageRegion[1][0] >  ImageTh1;
			JudgeBit1[4] = pImageRegion[1][2] >  ImageTh1;

			JudgeBit1[5] = pImageRegion[2][0] >  ImageTh1;
			JudgeBit1[6] = pImageRegion[2][1] >  ImageTh1;
			JudgeBit1[7] = pImageRegion[2][2] >  ImageTh1;


			int SumJudge = 0;
			int SumJudge1 = 0;

#pragma unroll
			for (int jcnt = 0; jcnt < JudgeBitNum; jcnt++)
			{
				SumJudge += JudgeBit[jcnt];
				SumJudge1 += JudgeBit1[jcnt];
			}

			if ((SumJudge == JudgeBitNum) && (SumJudge1 > 4))
			{
				int CurImageId = CurYPos / ImageHigh;

				// seperate the whole memory for each frame
				int curRegionPos = atomicAdd(&d_ROINumPerImage[CurImageId], 1);

				d_ROIPosArray[(CurImageId*MaxFluoNumPerImage + curRegionPos) * RegionPosInfNum + 0] = CurXPos;
				d_ROIPosArray[(CurImageId*MaxFluoNumPerImage + curRegionPos) * RegionPosInfNum + 1] = CurYPos;

//				printf("roi:%d %d\n", CurXPos, CurYPos);
			}
		}
	}

}

void SubregionExtractionLD(unsigned short *d_RawImg, unsigned short *d_ROIMem, int *d_ROIPosArray, int ROISize, int RegionNum_CurImage, int ImageWidth, int BatchedImageHigh, int ImageHigh, int StartFrame, cudaStream_t cstream)
{
	int BlockDim = ThreadsPerBlock;
	int BlockNum = (RegionNum_CurImage + ThreadsPerBlock - 1) / ThreadsPerBlock;


	// region extraction with point position
	gpuSubregionExtractionLD << < BlockNum, BlockDim, 0, cstream >> > (d_RawImg, d_ROIMem, d_ROIPosArray, ROISize, RegionNum_CurImage, ImageWidth, BatchedImageHigh, ImageHigh, StartFrame);


}

__global__ void gpuSubregionExtractionLD(unsigned short *d_RawImg, unsigned short *d_ROIMem, int *d_ROIPosArray, int ROISize, int RegionNum_CurImage, int ImageWidth, int BatchedImageHigh, int ImageHigh, int StartFrame)
{
	int gid = threadIdx.x + blockDim.x*blockIdx.x;

	if (gid < RegionNum_CurImage)
	{
		int ROIDataLen = ROISize*(ROISize + 1);
		int ROISize_Half = ROISize / 2;


		int HaflRegionSize = ROISize / 2;

		int RegionAddrOffset = gid*ROIDataLen;

		int XPos = d_ROIPosArray[gid * 2 + 0];
		int YPos = d_ROIPosArray[gid * 2 + 1];

		if ((XPos > ROISize_Half) && (YPos > ROISize_Half) && (XPos < ImageWidth - ROISize_Half) && (YPos < BatchedImageHigh - ROISize_Half))
		{
			for (int ycnt = 0; ycnt < ROISize; ycnt++)
			{
				for (int xcnt = 0; xcnt < ROISize; xcnt++)
				{
					int XPos_GMem = XPos - ROISize_Half + xcnt;
					int YPos_GMem = YPos - ROISize_Half + ycnt;

					d_ROIMem[RegionAddrOffset + ycnt*ROISize + xcnt] = d_RawImg[YPos_GMem*ImageWidth + XPos_GMem];
				}
			}


//			printf("roi:%d %d\n", XPos, YPos);


			int AddrOffset = ROISize*ROISize;
			int CurFrame = StartFrame + (YPos / ImageHigh);

			d_ROIMem[RegionAddrOffset + AddrOffset + 0] = XPos;
			d_ROIMem[RegionAddrOffset + AddrOffset + 1] = (YPos % ImageHigh);
			d_ROIMem[RegionAddrOffset + AddrOffset + 2] = CurFrame % 65536; //
			d_ROIMem[RegionAddrOffset + AddrOffset + 3] = CurFrame / 65536;
			d_ROIMem[RegionAddrOffset + AddrOffset + 4] = 0;

		}
	}
}



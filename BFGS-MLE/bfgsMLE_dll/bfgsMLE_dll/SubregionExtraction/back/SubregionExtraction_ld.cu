#include "SubregionExtraction_ld.h"


// raw to Gaussian filtered img, and raw to std image
__global__ void gpuImgFilleringLD1(unsigned short *d_RawImg, unsigned short *d_GaussImg, unsigned short *d_StdImg, int ImageWidth, int ImageHigh)
{
	// side pixels, valid pixels, side pixels
	__shared__ unsigned short ImageRegion[BatchLineNum + 6][ThreadsPerBlock]; 

	int gid = threadIdx.x + blockDim.x*blockIdx.x;
	int tid = threadIdx.x;

	int curx, cury;
	int cury_t;

	int ImgWidthBlockNum = (ImageWidth + ValidProcCore_7x7 - 1) / ValidProcCore_7x7;
	int ImgWidthThreadNum = ImgWidthBlockNum * 32;

	int CurColThread = gid % ImgWidthThreadNum;
	int CurRowThread = gid / ImgWidthThreadNum;

	int CurWidthBlock = CurColThread / 32;

	curx = CurWidthBlock * ValidProcCore_7x7 + tid - MarginSize_7x7;
	cury = CurRowThread*BatchLineNum - MarginSize_7x7;

	int lcnt = 0;

	//load image from global mem to shared mem

#pragma unroll
	for (lcnt = 0; lcnt < BatchLineNum + MarginSize_7x7*2; lcnt++)
	{
		cury_t = cury + lcnt;

		if ((curx >= 0) && (cury_t >= 0) && (curx<ImageWidth) && (cury_t<ImageHigh))
		{
			ImageRegion[lcnt][tid] = d_RawImg[cury_t*ImageWidth + curx];
		}
		else
		{
			ImageRegion[lcnt][tid] = 0;
		}
	}
	__syncthreads();


	int ProcLine0;
	int ProcLine1;
	int ProcLine2;
	int ProcLine3;
	int ProcLine4;

	float tdat1[6];
	float tdat2[6];

	float tempdat;
	float AvgDat;

	float curStd;

	// image filtering

	curx = CurWidthBlock * ValidProcCore_7x7 + tid;
	cury = CurRowThread*BatchLineNum;
#pragma unroll
	for (lcnt = 0; lcnt < BatchLineNum; lcnt++)
	{
		cury_t = cury + lcnt;
		if ((tid < ValidProcCore_7x7) && (curx<ImageWidth) && (cury_t<ImageHigh))
		{
			// 5x5 filter
			ProcLine0 = ImageRegion[lcnt + 1][tid + 1] * 78 + ImageRegion[lcnt + 2][tid + 1]  * 215 + ImageRegion[lcnt + 3][tid + 1] * 293 + ImageRegion[lcnt + 4][tid + 1]  * 215 + ImageRegion[lcnt + 5][tid + 1] * 78;
			ProcLine1 = ImageRegion[lcnt + 1][tid + 2] * 215 + ImageRegion[lcnt + 2][tid + 2] * 586 + ImageRegion[lcnt + 3][tid + 2] * 820 + ImageRegion[lcnt + 4][tid + 2]  * 586 + ImageRegion[lcnt + 5][tid + 2] * 215;
			ProcLine2 = ImageRegion[lcnt + 1][tid + 3] * 293 + ImageRegion[lcnt + 2][tid + 3] * 820 + ImageRegion[lcnt + 3][tid + 3] * 1172 + ImageRegion[lcnt + 4][tid + 3] * 820 + ImageRegion[lcnt + 5][tid + 3] * 293;
			ProcLine3 = ImageRegion[lcnt + 1][tid + 4] * 215 + ImageRegion[lcnt + 2][tid + 4] * 586 + ImageRegion[lcnt + 3][tid + 4] * 820 + ImageRegion[lcnt + 4][tid + 4]  * 586 + ImageRegion[lcnt + 5][tid + 4] * 215;
			ProcLine4 = ImageRegion[lcnt + 1][tid + 5] * 78 + ImageRegion[lcnt + 2][tid + 5]  * 215 + ImageRegion[lcnt + 3][tid + 5] * 293 + ImageRegion[lcnt + 4][tid + 5]  * 215 + ImageRegion[lcnt + 5][tid + 5] * 78;

			tempdat = (ProcLine0 + ProcLine1 + ProcLine2 + ProcLine3 + ProcLine4)/10000;

			d_GaussImg[cury_t*ImageWidth + curx] = tempdat;


			// std filter

			tdat1[0] = (ImageRegion[lcnt + 0][tid + 1] < ImageRegion[lcnt + 1][tid + 6]) ? ImageRegion[lcnt + 0][tid + 1] : ImageRegion[lcnt + 1][tid + 6];
			tdat1[1] = (ImageRegion[lcnt + 0][tid + 2] < ImageRegion[lcnt + 2][tid + 6]) ? ImageRegion[lcnt + 0][tid + 2] : ImageRegion[lcnt + 2][tid + 6];
			tdat1[2] = (ImageRegion[lcnt + 0][tid + 3] < ImageRegion[lcnt + 3][tid + 6]) ? ImageRegion[lcnt + 0][tid + 3] : ImageRegion[lcnt + 3][tid + 6];
			tdat1[3] = (ImageRegion[lcnt + 0][tid + 4] < ImageRegion[lcnt + 4][tid + 6]) ? ImageRegion[lcnt + 0][tid + 4] : ImageRegion[lcnt + 4][tid + 6];
			tdat1[4] = (ImageRegion[lcnt + 0][tid + 5] < ImageRegion[lcnt + 5][tid + 6]) ? ImageRegion[lcnt + 0][tid + 5] : ImageRegion[lcnt + 5][tid + 6];
			tdat1[5] = (ImageRegion[lcnt + 0][tid + 6] < ImageRegion[lcnt + 6][tid + 6]) ? ImageRegion[lcnt + 0][tid + 6] : ImageRegion[lcnt + 6][tid + 6];

			tdat2[0] = (ImageRegion[lcnt + 0][tid + 0] < ImageRegion[lcnt + 6][tid + 0]) ? ImageRegion[lcnt + 0][tid + 0] : ImageRegion[lcnt + 6][tid + 0];
			tdat2[1] = (ImageRegion[lcnt + 1][tid + 0] < ImageRegion[lcnt + 6][tid + 1]) ? ImageRegion[lcnt + 1][tid + 0] : ImageRegion[lcnt + 6][tid + 1];
			tdat2[2] = (ImageRegion[lcnt + 2][tid + 0] < ImageRegion[lcnt + 6][tid + 2]) ? ImageRegion[lcnt + 2][tid + 0] : ImageRegion[lcnt + 6][tid + 2];
			tdat2[3] = (ImageRegion[lcnt + 3][tid + 0] < ImageRegion[lcnt + 6][tid + 3]) ? ImageRegion[lcnt + 3][tid + 0] : ImageRegion[lcnt + 6][tid + 3];
			tdat2[4] = (ImageRegion[lcnt + 4][tid + 0] < ImageRegion[lcnt + 6][tid + 4]) ? ImageRegion[lcnt + 4][tid + 0] : ImageRegion[lcnt + 6][tid + 4];
			tdat2[5] = (ImageRegion[lcnt + 5][tid + 0] < ImageRegion[lcnt + 6][tid + 5]) ? ImageRegion[lcnt + 5][tid + 0] : ImageRegion[lcnt + 6][tid + 5];

			AvgDat = (tdat1[0] + tdat1[1] + tdat1[2] + tdat1[3] + tdat1[4] + tdat1[5] + tdat2[0] + tdat2[1] + tdat2[2] + tdat2[3] + tdat2[4] + tdat2[5]) / 12.0f;



			ProcLine0 =
				(tdat1[0] - AvgDat)*(tdat1[0] - AvgDat) +
				(tdat1[1] - AvgDat)*(tdat1[1] - AvgDat) +
				(tdat1[2] - AvgDat)*(tdat1[2] - AvgDat) +
				(tdat1[3] - AvgDat)*(tdat1[3] - AvgDat) +
				(tdat1[4] - AvgDat)*(tdat1[4] - AvgDat) +
				(tdat1[5] - AvgDat)*(tdat1[5] - AvgDat);
			ProcLine1 =
				(tdat2[0] - AvgDat)*(tdat2[0] - AvgDat) +
				(tdat2[1] - AvgDat)*(tdat2[1] - AvgDat) +
				(tdat2[2] - AvgDat)*(tdat2[2] - AvgDat) +
				(tdat2[3] - AvgDat)*(tdat2[3] - AvgDat) +
				(tdat2[4] - AvgDat)*(tdat2[4] - AvgDat) +
				(tdat2[5] - AvgDat)*(tdat2[5] - AvgDat);

				
			curStd = __fsqrt_rn((ProcLine0 + ProcLine1) / 12.0f);



			if (curStd < 10.0f)curStd = 10.0f;

			d_StdImg[cury_t*ImageWidth + curx] = curStd;

		}
	}	
}


// d_GaussImg to d_AnnuImg
__global__ void gpuImgFilleringLD2(unsigned short *d_GaussImg, unsigned short *d_AnnuImg, int ImageWidth, int ImageHigh)
{
	__shared__ unsigned short ImageRegion[BatchLineNum + 6][ThreadsPerBlock];

	int gid = threadIdx.x + blockDim.x*blockIdx.x;
	int tid = threadIdx.x;

	int curx, cury;
	int cury_t;

	int ImgWidthBlockNum = (ImageWidth + ValidProcCore_7x7 - 1) / ValidProcCore_7x7;
	int ImgWidthThreadNum = ImgWidthBlockNum * 32;

	int CurColThread = gid % ImgWidthThreadNum;
	int CurRowThread = gid / ImgWidthThreadNum;

	int CurWidthBlock = CurColThread / 32;

	curx = CurWidthBlock * ValidProcCore_7x7 + tid - MarginSize_7x7;
	cury = CurRowThread*BatchLineNum - MarginSize_7x7;

	int lcnt = 0;

	//load image from global mem to shared mem

#pragma unroll
	for (lcnt = 0; lcnt < BatchLineNum + 6; lcnt++)
	{
		cury_t = cury + lcnt;
		if ((curx >= 0) && (cury_t >= 0) && (curx<ImageWidth) && (cury_t<ImageHigh))
		{
			ImageRegion[lcnt][tid] = d_GaussImg[cury_t*ImageWidth + curx];
		}
		else
		{
			ImageRegion[lcnt][tid] = 0;
		}		
	}

	__syncthreads();

	int ProcLine0;
	int ProcLine1;
	int ProcLine2;
	int ProcLine3;

	float tdat1[6];
	float tdat2[6];


	int tempdat;
	int centdat;

	curx = CurWidthBlock * ValidProcCore_7x7 + tid;
	cury = CurRowThread*BatchLineNum;
#pragma unroll
	for (lcnt = 0; lcnt < BatchLineNum; lcnt++)
	{
		cury_t = cury + lcnt;
		if ((tid < ValidProcCore_7x7) && (curx < ImageWidth) && (cury_t < ImageHigh))
		{
			// 7x7 filter

			ProcLine0 = ImageRegion[lcnt + 0][tid + 1] + ImageRegion[lcnt + 0][tid + 2] + ImageRegion[lcnt + 0][tid + 3] + ImageRegion[lcnt + 0][tid + 4] + ImageRegion[lcnt + 0][tid + 5] + ImageRegion[lcnt + 0][tid + 6];
			ProcLine1 = ImageRegion[lcnt + 1][tid + 6] + ImageRegion[lcnt + 2][tid + 6] + ImageRegion[lcnt + 3][tid + 6] + ImageRegion[lcnt + 4][tid + 6] + ImageRegion[lcnt + 5][tid + 6] + ImageRegion[lcnt + 6][tid + 6];
			ProcLine2 = ImageRegion[lcnt + 0][tid + 0] + ImageRegion[lcnt + 1][tid + 0] + ImageRegion[lcnt + 2][tid + 0] + ImageRegion[lcnt + 3][tid + 0] + ImageRegion[lcnt + 4][tid + 0] + ImageRegion[lcnt + 5][tid + 0];
			ProcLine3 = ImageRegion[lcnt + 6][tid + 0] + ImageRegion[lcnt + 6][tid + 1] + ImageRegion[lcnt + 6][tid + 2] + ImageRegion[lcnt + 6][tid + 3] + ImageRegion[lcnt + 6][tid + 4] + ImageRegion[lcnt + 6][tid + 5];
			
			/**/
			// select 6 max pixels
			tdat1[0] = (ImageRegion[lcnt + 0][tid + 1] > ImageRegion[lcnt + 1][tid + 6]) ? ImageRegion[lcnt + 0][tid + 1] : ImageRegion[lcnt + 1][tid + 6];
			tdat1[1] = (ImageRegion[lcnt + 0][tid + 2] > ImageRegion[lcnt + 2][tid + 6]) ? ImageRegion[lcnt + 0][tid + 2] : ImageRegion[lcnt + 2][tid + 6];
			tdat1[2] = (ImageRegion[lcnt + 0][tid + 3] > ImageRegion[lcnt + 3][tid + 6]) ? ImageRegion[lcnt + 0][tid + 3] : ImageRegion[lcnt + 3][tid + 6];
			tdat1[3] = (ImageRegion[lcnt + 0][tid + 4] > ImageRegion[lcnt + 4][tid + 6]) ? ImageRegion[lcnt + 0][tid + 4] : ImageRegion[lcnt + 4][tid + 6];
			tdat1[4] = (ImageRegion[lcnt + 0][tid + 5] > ImageRegion[lcnt + 5][tid + 6]) ? ImageRegion[lcnt + 0][tid + 5] : ImageRegion[lcnt + 5][tid + 6];
			tdat1[5] = (ImageRegion[lcnt + 0][tid + 6] > ImageRegion[lcnt + 6][tid + 6]) ? ImageRegion[lcnt + 0][tid + 6] : ImageRegion[lcnt + 6][tid + 6];

			tdat2[0] = (ImageRegion[lcnt + 0][tid + 0] > ImageRegion[lcnt + 6][tid + 0]) ? ImageRegion[lcnt + 0][tid + 0] : ImageRegion[lcnt + 6][tid + 0];
			tdat2[1] = (ImageRegion[lcnt + 1][tid + 0] > ImageRegion[lcnt + 6][tid + 1]) ? ImageRegion[lcnt + 1][tid + 0] : ImageRegion[lcnt + 6][tid + 1];
			tdat2[2] = (ImageRegion[lcnt + 2][tid + 0] > ImageRegion[lcnt + 6][tid + 2]) ? ImageRegion[lcnt + 2][tid + 0] : ImageRegion[lcnt + 6][tid + 2];
			tdat2[3] = (ImageRegion[lcnt + 3][tid + 0] > ImageRegion[lcnt + 6][tid + 3]) ? ImageRegion[lcnt + 3][tid + 0] : ImageRegion[lcnt + 6][tid + 3];
			tdat2[4] = (ImageRegion[lcnt + 4][tid + 0] > ImageRegion[lcnt + 6][tid + 4]) ? ImageRegion[lcnt + 4][tid + 0] : ImageRegion[lcnt + 6][tid + 4];
			tdat2[5] = (ImageRegion[lcnt + 5][tid + 0] > ImageRegion[lcnt + 6][tid + 5]) ? ImageRegion[lcnt + 5][tid + 0] : ImageRegion[lcnt + 6][tid + 5];

			tdat1[0] = (tdat1[0] > tdat2[0]) ? tdat1[0] : tdat2[0];
			tdat1[1] = (tdat1[1] > tdat2[1]) ? tdat1[1] : tdat2[1];
			tdat1[2] = (tdat1[2] > tdat2[2]) ? tdat1[2] : tdat2[2];
			tdat1[3] = (tdat1[3] > tdat2[3]) ? tdat1[3] : tdat2[3];
			tdat1[4] = (tdat1[4] > tdat2[4]) ? tdat1[4] : tdat2[4];
			tdat1[5] = (tdat1[5] > tdat2[5]) ? tdat1[5] : tdat2[5];

			tdat1[0] = (tdat1[0] > tdat1[3]) ? tdat1[0] : tdat1[3];
			tdat1[1] = (tdat1[1] > tdat1[4]) ? tdat1[1] : tdat1[4];
			tdat1[2] = (tdat1[2] > tdat1[5]) ? tdat1[2] : tdat1[5];

			tempdat = tdat1[0] + tdat1[1] + tdat1[2];
			tempdat = (ProcLine0 + ProcLine1 + ProcLine2 + ProcLine3 - tempdat) / 21;
			

//			tempdat = (ProcLine0 + ProcLine1 + ProcLine2 + ProcLine3) / 24;

			centdat = ImageRegion[lcnt + 3][tid + 3];

			if (centdat > tempdat)tempdat = centdat - tempdat;
			else tempdat = 0;

			d_AnnuImg[cury_t*ImageWidth + curx] = tempdat;
		}
	}
}

__global__ void gpuSubregionFindingLD(unsigned short *d_AnnuImg, unsigned short *d_StdImg, int *d_RegionPosMem, int *d_RegionNum, int RawImageWidth, int GroupImageHigh, int RawImgHigh)
{
	__shared__ unsigned short AnnuImgRegion[BatchLineNum + 2][ThreadsPerBlock];
	__shared__ unsigned short StdImgRegion[BatchLineNum + 2][ThreadsPerBlock];

	int gid = threadIdx.x + blockDim.x*blockIdx.x;
	int tid = threadIdx.x;
	//
	const int CurGroupNum = GroupImageHigh / RawImgHigh;
	const int EachGroupFluoNum = MaxPointNum / CurGroupNum;
	int CurGroupId = 0;
	//
	int curx, cury;
	int cury_t;

	int ImgWidthBlockNum = (RawImageWidth + ValidProcCore_3x3 - 1) / ValidProcCore_3x3;
	int ImgWidthThreadNum = ImgWidthBlockNum * 32;

	int CurColThread = gid % ImgWidthThreadNum;
	int CurRowThread = gid / ImgWidthThreadNum;

	int CurWidthBlock = CurColThread / 32;

	curx = CurWidthBlock * ValidProcCore_3x3 + tid - MarginSize_3x3;
	cury = CurRowThread*BatchLineNum - MarginSize_3x3;

	int lcnt = 0;

	//load image from global mem to shared mem

#pragma unroll
	for (lcnt = 0; lcnt < BatchLineNum + 2; lcnt++)
	{
		cury_t = cury + lcnt;
		if ((curx >= 0) && (cury_t >= 0) && (curx<RawImageWidth) && (cury_t<GroupImageHigh))
		{
			AnnuImgRegion[lcnt][tid] = d_AnnuImg[cury_t*RawImageWidth + curx];
			StdImgRegion[lcnt][tid] = d_StdImg[cury_t*RawImageWidth + curx];
		}
		else
		{
			AnnuImgRegion[lcnt][tid] = 0;
			StdImgRegion[lcnt][tid] = 0;
		}
	}

	__syncthreads();


	int CurStd;
	// sum of pixel value of 3x3 region
	int sum9;
	int sum5; 

	int CenterPixel;

#define JudgeBitNum				11

	int JudgeBit[JudgeBitNum]; // judgement memory
	int SumJudge;

	int curRegionPos;
	int jcnt = 0;

	curx = CurWidthBlock * ValidProcCore_3x3 + tid;
	cury = CurRowThread*BatchLineNum;
#pragma unroll
	for (lcnt = 0; lcnt < BatchLineNum; lcnt++)
	{
		cury_t = cury + lcnt;
		if ((tid < ValidProcCore_3x3) && (curx < RawImageWidth) && (cury_t < GroupImageHigh))
		{


			CurStd = (StdImgRegion[lcnt + 0][tid + 0] + StdImgRegion[lcnt + 0][tid + 1] + StdImgRegion[lcnt + 0][tid + 2] +
				StdImgRegion[lcnt + 1][tid + 0] + StdImgRegion[lcnt + 1][tid + 1] + StdImgRegion[lcnt + 1][tid + 2] +
				StdImgRegion[lcnt + 2][tid + 0] + StdImgRegion[lcnt + 2][tid + 1] + StdImgRegion[lcnt + 2][tid + 2]) / 9;

			// sub region judgement

			sum9 = AnnuImgRegion[lcnt + 0][tid + 0] + AnnuImgRegion[lcnt + 0][tid + 1] + AnnuImgRegion[lcnt + 0][tid + 2] +
				AnnuImgRegion[lcnt + 1][tid + 0] + AnnuImgRegion[lcnt + 1][tid + 1] + AnnuImgRegion[lcnt + 1][tid + 2] +
				AnnuImgRegion[lcnt + 2][tid + 0] + AnnuImgRegion[lcnt + 2][tid + 1] + AnnuImgRegion[lcnt + 2][tid + 2];

			sum5 = AnnuImgRegion[lcnt + 0][tid + 1] +
				AnnuImgRegion[lcnt + 1][tid + 0] + AnnuImgRegion[lcnt + 1][tid + 1] + AnnuImgRegion[lcnt + 1][tid + 2] +
				AnnuImgRegion[lcnt + 2][tid + 1];

			CenterPixel = AnnuImgRegion[lcnt + 1][tid + 1];

			// subregion judgement
			JudgeBit[0] = CenterPixel >  AnnuImgRegion[lcnt + 0][tid + 0];
			JudgeBit[1] = CenterPixel >= AnnuImgRegion[lcnt + 1][tid + 0];
			JudgeBit[2] = CenterPixel >  AnnuImgRegion[lcnt + 2][tid + 0];

			JudgeBit[3] = CenterPixel >  AnnuImgRegion[lcnt + 0][tid + 1];
			JudgeBit[4] = CenterPixel >= AnnuImgRegion[lcnt + 2][tid + 1];

			JudgeBit[5] = CenterPixel >  AnnuImgRegion[lcnt + 0][tid + 2];
			JudgeBit[6] = CenterPixel >  AnnuImgRegion[lcnt + 1][tid + 2];
			JudgeBit[7] = CenterPixel >  AnnuImgRegion[lcnt + 2][tid + 2];



			JudgeBit[8] = CenterPixel >= CurStd * 3.0f;
			JudgeBit[9] = sum9 >= CurStd * 13.0f;
			JudgeBit[10] = sum5 >= CurStd * 9.0f;

			SumJudge = 0;

#pragma unroll
			for (jcnt = 0; jcnt < JudgeBitNum; jcnt++)
			{
				SumJudge += JudgeBit[jcnt];
			}

			if (SumJudge == JudgeBitNum)
			{
				if ((curx >= 4) && (cury_t >= 4) && (curx < RawImageWidth - 4) && (cury_t < GroupImageHigh - 4))
				{
					CurGroupId = cury_t / RawImgHigh;

					// seperate the whole memory for each frame
					curRegionPos = atomicAdd(&d_RegionNum[CurGroupId], 1);

					d_RegionPosMem[(CurGroupId*EachGroupFluoNum + curRegionPos) * RegionPosInfNum + 0] = curx;
					d_RegionPosMem[(CurGroupId*EachGroupFluoNum + curRegionPos) * RegionPosInfNum + 1] = cury_t;

				}
			}
		}
	}
}


__global__ void gpuSubregionExtractionLD(unsigned short *d_RawImg, unsigned short *d_RegionMem, int *d_RegionPosMem, int FluoNum, int ROISize, int ImageWidth, int GroupImgHigh, int RawImgHigh, int StartFrame)
{
	int gid = threadIdx.x + blockDim.x*blockIdx.x;

	if (gid < FluoNum)
	{
		int SubRegionWidth = ROISize*(ROISize + 1);
		int HalfRegionSize = ROISize / 2;

		int curRegionPos = 0;
		int RegionAddrOffset;

		int HaflRegionSize = ROISize / 2;

		RegionAddrOffset = gid*SubRegionWidth;

		int XPos = d_RegionPosMem[gid * 2 + 0];
		int YPos = d_RegionPosMem[gid * 2 + 1];

		if ((XPos > HalfRegionSize) && (YPos > HalfRegionSize) && (XPos < ImageWidth - HalfRegionSize) && (YPos < GroupImgHigh - HalfRegionSize))
		{
			int curx, cury;

			int xcnt, ycnt;

			int AddrInc;

			for (ycnt = 0; ycnt < ROISize; ycnt++)
			{
				AddrInc = ycnt*ROISize;
				cury = YPos - HaflRegionSize + ycnt;

				for (xcnt = 0; xcnt < ROISize; xcnt++)
				{
					curx = XPos - HaflRegionSize + xcnt;
					d_RegionMem[RegionAddrOffset + AddrInc + xcnt] = d_RawImg[cury*ImageWidth + curx];
				}
			}

			AddrInc = ROISize*ROISize;
			d_RegionMem[RegionAddrOffset + AddrInc + 0] = XPos;
			d_RegionMem[RegionAddrOffset + AddrInc + 1] = (YPos % RawImgHigh);
			d_RegionMem[RegionAddrOffset + AddrInc + 2] = StartFrame + (YPos / RawImgHigh); //
			d_RegionMem[RegionAddrOffset + AddrInc + 3] = 0;
			d_RegionMem[RegionAddrOffset + AddrInc + 4] = 0;

		}
	}
}


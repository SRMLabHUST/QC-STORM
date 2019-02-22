#include "SubregionExtraction_hd.h"

#include <stdio.h>
#include <string.h>



#define LineBlockNum			(ThreadsPerBlock)
#define LineBlockNum_7x7		(LineBlockNum + 7) // LineBlockNum +6, here to avoid bank conflict
#define ExtPixelNum_7x7			6

#define LineBlockNum_3x3		(LineBlockNum + 3) // LineBlockNum +6, here to avoid bank conflict
#define ExtPixelNum_3x3			2


__global__ void gpuSubregionExtraction(unsigned short *d_RawImg, int * d_RegionNum, int * d_RegionMem, unsigned short * d_RegionMem, int RegionSize, int LargeImgWidth, int LargeImgHigh, int ImgHighAdj, int FrameOffset)
{
	int gid = threadIdx.x + blockDim.x*blockIdx.x;
//	int tid = threadIdx.x;
	int CurFrame;


	if (gid < (*d_RegionNum))
	{
		int WholeImageWidth = RegionSize*(RegionSize + 1);

		int RegionAddrOffset = gid * WholeImageWidth;

		int HaflRegionSize = RegionSize / 2;

		int XPos = d_RegionMem[gid * 3 + 0];
		int YPos = d_RegionMem[gid * 3 + 1];
		int CurStd = d_RegionMem[gid * 3 + 2];

		CurFrame = YPos / ImgHighAdj + FrameOffset;
		int RealYPos = YPos % ImgHighAdj;

		int curx, cury;

		int xcnt, ycnt;

		int AddrInc;

		for (ycnt = 0; ycnt < RegionSize; ycnt++)
		{
			AddrInc = ycnt*RegionSize;
			cury = YPos - HaflRegionSize + ycnt;

			for (xcnt = 0; xcnt < RegionSize; xcnt++)
			{
				curx = XPos - HaflRegionSize + xcnt;
				d_RegionMem[RegionAddrOffset + AddrInc + xcnt] = d_RawImg[cury*LargeImgWidth + curx];
			}
		}
		AddrInc = RegionSize*RegionSize;
		d_RegionMem[RegionAddrOffset + AddrInc + 0] = XPos;
		d_RegionMem[RegionAddrOffset + AddrInc + 1] = RealYPos;
		d_RegionMem[RegionAddrOffset + AddrInc + 2] = CurFrame + 1;
		d_RegionMem[RegionAddrOffset + AddrInc + 3] = 0;
		d_RegionMem[RegionAddrOffset + AddrInc + 4] = 0;
		d_RegionMem[RegionAddrOffset + AddrInc + 5] = 0;
		d_RegionMem[RegionAddrOffset + AddrInc + 6] = 0;
		d_RegionMem[RegionAddrOffset + AddrInc + 7] = 0;
		d_RegionMem[RegionAddrOffset + AddrInc + 8] = CurStd;

	}
}

// gpuImgFillering
// calculate the dog filtered image and std image simultaneously

__global__ void gpuImgFiller_gauss_std(unsigned short * d_RawImg, unsigned short * d_GaussImg, unsigned short * d_StdImg1, int ImageWidth, int ImageHigh)

{
	// first, gaussian filter and Std Filler
	// side pixels, valid pixels, side pixels
	__shared__ float ImageRegion[ThreadsPerBlock + 6][LineBlockNum_7x7]; // 128-6 lines

	int gid = threadIdx.x + blockDim.x*blockIdx.x;
	int tid = threadIdx.x;
	int yoffset = blockDim.x*blockIdx.x;

	int curx;
	int cury;
	int xoffset;

	int LineProcNum = (ImageWidth + LineBlockNum - 1) / LineBlockNum;
	int lcnt = 0;
	int xcnt = 0;

	float ProcLine0;
	float ProcLine1;
	float ProcLine2;
	float ProcLine3;
	float ProcLine4;

	unsigned short minDat1[28];
	unsigned short(*pMinDat1)[7] = (unsigned short(*)[7])minDat1;

	unsigned short minDat2[16];
	unsigned short(*pMinDat2)[4] = (unsigned short(*)[4])minDat2;

	int AvgDat;

	float curStd;
	float tempdat;
	int pcnt;

	// gid and tid is y pos, lcnt is xpos
	for (lcnt = 0; lcnt < LineProcNum; lcnt++)
	{
		xoffset = lcnt*LineBlockNum;

		// read global image to shared image
		curx = xoffset + tid - 3;
		for (xcnt = 0; xcnt < ThreadsPerBlock + 6; xcnt++)
		{
			// load 32 pixel
			cury = yoffset + xcnt - 3;

			if ((curx >= 0) && (curx < ImageWidth) && (cury >= 0) && (cury < ImageHigh))
			{
				ImageRegion[xcnt][tid] = d_RawImg[cury*ImageWidth + curx];
			}
			else
			{
				ImageRegion[xcnt][tid] = 0;
			}
		}
		curx = xoffset + ThreadsPerBlock + tid - 3;
		for (xcnt = 0; xcnt < LineBlockNum + 6; xcnt++)
		{
			// load extra 6 margin pixel
			cury = yoffset + xcnt - 3;
			if(tid < 6)
			{
				if ((curx < ImageWidth) && (cury >= 0) && (cury < ImageHigh)) // (curx >= 0) && 
				{
					ImageRegion[xcnt][tid+ThreadsPerBlock] = d_RawImg[cury*ImageWidth + curx];
				}
				else
				{
					ImageRegion[xcnt][tid+ThreadsPerBlock] = 0;
				}
			}
		}

		__syncthreads();
		// image filtering

	
		cury = gid;

		for (xcnt = 0; xcnt < LineBlockNum; xcnt++)// y
		{
			curx = xoffset + xcnt;
			if ((curx < ImageWidth) && (cury < ImageHigh))
			{
				// gaussian filter
				// 5x5 filter

				ProcLine0 = ImageRegion[tid + 1][xcnt + 1]*0.0030f + ImageRegion[tid + 2][xcnt + 1]*0.0133f + ImageRegion[tid + 3][xcnt + 1]*0.0219f + ImageRegion[tid + 4][xcnt + 1]*0.0133f + ImageRegion[tid + 5][xcnt + 1]*0.0030f;
				ProcLine1 = ImageRegion[tid + 1][xcnt + 2]*0.0133f + ImageRegion[tid + 2][xcnt + 2]*0.0596f + ImageRegion[tid + 3][xcnt + 2]*0.0983f + ImageRegion[tid + 4][xcnt + 2]*0.0596f + ImageRegion[tid + 5][xcnt + 2]*0.0133f;
				ProcLine2 = ImageRegion[tid + 1][xcnt + 3]*0.0219f + ImageRegion[tid + 2][xcnt + 3]*0.0983f + ImageRegion[tid + 3][xcnt + 3]*0.1621f + ImageRegion[tid + 4][xcnt + 3]*0.0983f + ImageRegion[tid + 5][xcnt + 3]*0.0219f;
				ProcLine3 = ImageRegion[tid + 1][xcnt + 4]*0.0133f + ImageRegion[tid + 2][xcnt + 4]*0.0596f + ImageRegion[tid + 3][xcnt + 4]*0.0983f + ImageRegion[tid + 4][xcnt + 4]*0.0596f + ImageRegion[tid + 5][xcnt + 4]*0.0133f;
				ProcLine4 = ImageRegion[tid + 1][xcnt + 5]*0.0030f + ImageRegion[tid + 2][xcnt + 5]*0.0133f + ImageRegion[tid + 3][xcnt + 5]*0.0219f + ImageRegion[tid + 4][xcnt + 5]*0.0133f + ImageRegion[tid + 5][xcnt + 5]*0.0030f;

				tempdat = ProcLine0 + ProcLine1 + ProcLine2 + ProcLine3 + ProcLine4;


				if(tempdat < 0.0f) tempdat = 0.0f;

				d_GaussImg[cury*ImageWidth + curx] = tempdat;


				// std1 filter
				for(pcnt=0; pcnt<7; pcnt++)
				{
					if( ImageRegion[tid + 0][xcnt + pcnt] < ImageRegion[tid + 3][xcnt + pcnt]) 	pMinDat1[0][pcnt] = ImageRegion[tid + 0][xcnt + pcnt];
					else 																		pMinDat1[0][pcnt] = ImageRegion[tid + 3][xcnt + pcnt];

					if( ImageRegion[tid + 1][xcnt + pcnt] < ImageRegion[tid + 4][xcnt + pcnt]) 	pMinDat1[1][pcnt] = ImageRegion[tid + 1][xcnt + pcnt];
					else 																		pMinDat1[1][pcnt] = ImageRegion[tid + 4][xcnt + pcnt];

					if( ImageRegion[tid + 2][xcnt + pcnt] < ImageRegion[tid + 5][xcnt + pcnt]) 	pMinDat1[2][pcnt] = ImageRegion[tid + 2][xcnt + pcnt];
					else 																		pMinDat1[2][pcnt] = ImageRegion[tid + 5][xcnt + pcnt];

					pMinDat1[3][pcnt] = ImageRegion[tid + 6][xcnt + pcnt];
				}

				for(pcnt=0; pcnt<4; pcnt++)
				{
					if(pMinDat1[pcnt][0] < pMinDat1[pcnt][3])	pMinDat2[pcnt][0] = pMinDat1[pcnt][0];
					else 										pMinDat2[pcnt][0] = pMinDat1[pcnt][3];

					if(pMinDat1[pcnt][1] < pMinDat1[pcnt][4])	pMinDat2[pcnt][1] = pMinDat1[pcnt][1];
					else 										pMinDat2[pcnt][1] = pMinDat1[pcnt][4];

					if(pMinDat1[pcnt][2] < pMinDat1[pcnt][5])	pMinDat2[pcnt][2] = pMinDat1[pcnt][2];
					else 										pMinDat2[pcnt][2] = pMinDat1[pcnt][5];

					pMinDat2[pcnt][3] = pMinDat1[pcnt][6];

				}

				AvgDat = (  minDat2[ 0] + minDat2[ 1] + minDat2[ 2] + minDat2[ 3] + 
							minDat2[ 4] + minDat2[ 5] + minDat2[ 6] + minDat2[ 7] + 
							minDat2[ 8] + minDat2[ 9] + minDat2[10] + minDat2[11] + 
							minDat2[12] + minDat2[13] + minDat2[14] + minDat2[15])/16;


				ProcLine0 = abs(minDat2[ 0] - AvgDat) +  abs(minDat2[ 1] - AvgDat) +  abs(minDat2[ 2] - AvgDat) +  abs(minDat2[ 3] - AvgDat);
				ProcLine1 = abs(minDat2[ 4] - AvgDat) +  abs(minDat2[ 5] - AvgDat) +  abs(minDat2[ 6] - AvgDat) +  abs(minDat2[ 7] - AvgDat);
				ProcLine2 = abs(minDat2[ 8] - AvgDat) +  abs(minDat2[ 9] - AvgDat) +  abs(minDat2[10] - AvgDat) +  abs(minDat2[11] - AvgDat);
				ProcLine3 = abs(minDat2[12] - AvgDat) +  abs(minDat2[13] - AvgDat) +  abs(minDat2[14] - AvgDat) +  abs(minDat2[15] - AvgDat);


				curStd = (ProcLine0 + ProcLine1 + ProcLine2 + ProcLine3)/16;// / 16.0f;

				if (curStd < 10.0f)curStd = 10.0f;
				if (curStd > 80.0f)curStd = 80.0f;

				d_StdImg1[cury*ImageWidth + curx] = curStd;

			}
		}
	}
}

// gpuStdImgSmooth
// filter the std image by annular filter to remove light spots in it
__global__ void gpuImgFiller_annu_stdannu(unsigned short * d_GaussImg, unsigned short * d_StdImg1, unsigned short * d_AnnuImg, unsigned short * d_StdImg2, int ImageWidth, int ImageHigh)
{
	// 7x7 filter
	__shared__ unsigned short GaussImgRegion[ThreadsPerBlock + 6][LineBlockNum_7x7]; // 
	__shared__ unsigned short Std1ImgRegion[ThreadsPerBlock + 6][LineBlockNum_7x7];  // 

	int gid = threadIdx.x + blockDim.x*blockIdx.x;
	int tid = threadIdx.x;
	int yoffset = blockDim.x*blockIdx.x;


	int curx;
	int cury;
	int xoffset;

	int LineProcNum = (ImageWidth + LineBlockNum - 1) / LineBlockNum;
	int lcnt = 0;
	int xcnt = 0;

	int ProcLine0;
	int ProcLine1;
	int ProcLine2;
	int ProcLine3;
	int ProcLine4;


	int AvgDat;
	int CurDat;

	for (lcnt = 0; lcnt < LineProcNum; lcnt++)
	{
		xoffset = lcnt*LineBlockNum;

		// read global image to shared image
		curx = xoffset + tid - 3;
		for (xcnt = 0; xcnt < ThreadsPerBlock + 6; xcnt++)
		{
			// load 32 pixel
			cury = yoffset + xcnt - 3;

			if ((curx >= 0) && (curx < ImageWidth) && (cury >= 0) && (cury < ImageHigh))
			{
				GaussImgRegion[xcnt][tid] = d_GaussImg[cury*ImageWidth + curx];
				Std1ImgRegion[xcnt][tid]  = d_StdImg1[cury*ImageWidth + curx];
			}
			else
			{
				GaussImgRegion[xcnt][tid] = 0;
				Std1ImgRegion[xcnt][tid]  = 0;
			}
		}
		curx = xoffset + ThreadsPerBlock + tid - 3;
		for (xcnt = 0; xcnt < LineBlockNum + 6; xcnt++)
		{
			// load extra 6 margin pixel
			cury = yoffset + xcnt - 3;
			if(tid < 6)
			{
				if ((curx < ImageWidth) && (cury >= 0) && (cury < ImageHigh)) // (curx >= 0) && 
				{
					GaussImgRegion[xcnt][tid+ThreadsPerBlock]  = d_GaussImg[cury*ImageWidth + curx];
					Std1ImgRegion[xcnt][tid+ThreadsPerBlock]   = d_StdImg1[cury*ImageWidth + curx];
				}
				else
				{
					GaussImgRegion[xcnt][tid+ThreadsPerBlock] = 0;
					Std1ImgRegion[xcnt][tid+ThreadsPerBlock]  = 0;
				}
			}
		}

		__syncthreads();

		// image filtering 7x7

		cury = gid;
		for (xcnt = 0; xcnt < LineBlockNum; xcnt++)// y
		{
			curx = xoffset + xcnt;
			if ((curx < ImageWidth) && (cury < ImageHigh))
			{
				// common annular filter for fluo image to remove background
				ProcLine0 = GaussImgRegion[tid + 0][xcnt + 0] + GaussImgRegion[tid + 1][xcnt + 0] + GaussImgRegion[tid + 2][xcnt + 0] + GaussImgRegion[tid + 3][xcnt + 0] + GaussImgRegion[tid + 4][xcnt + 0] + GaussImgRegion[tid + 5][xcnt + 0] + GaussImgRegion[tid + 6][xcnt + 0];
				ProcLine1 = GaussImgRegion[tid + 0][xcnt + 1] + GaussImgRegion[tid + 0][xcnt + 2] + GaussImgRegion[tid + 0][xcnt + 3] + GaussImgRegion[tid + 0][xcnt + 4] + GaussImgRegion[tid + 0][xcnt + 5];
				ProcLine2 = GaussImgRegion[tid + 6][xcnt + 1] + GaussImgRegion[tid + 6][xcnt + 2] + GaussImgRegion[tid + 6][xcnt + 3] + GaussImgRegion[tid + 6][xcnt + 4] + GaussImgRegion[tid + 6][xcnt + 5];
				ProcLine3 = GaussImgRegion[tid + 0][xcnt + 6] + GaussImgRegion[tid + 1][xcnt + 6] + GaussImgRegion[tid + 2][xcnt + 6] + GaussImgRegion[tid + 3][xcnt + 6] + GaussImgRegion[tid + 4][xcnt + 6] + GaussImgRegion[tid + 5][xcnt + 6] + GaussImgRegion[tid + 6][xcnt + 6];

				AvgDat = (ProcLine0 + ProcLine1 + ProcLine2 + ProcLine3)/ 24;
				CurDat = GaussImgRegion[tid + 3][xcnt + 3] - AvgDat;

				if(CurDat < 0) CurDat = 0;
				
				d_AnnuImg[cury*ImageWidth + curx] = CurDat;


				// annular filter to remove bring point in the std image, and keeps background
				
				ProcLine0 =(Std1ImgRegion[tid + 0][xcnt + 0] + Std1ImgRegion[tid + 0][xcnt + 1] + Std1ImgRegion[tid + 0][xcnt + 2] + Std1ImgRegion[tid + 0][xcnt + 3] + Std1ImgRegion[tid + 0][xcnt + 4] + Std1ImgRegion[tid + 0][xcnt + 5] + Std1ImgRegion[tid + 0][xcnt + 6] +
							Std1ImgRegion[tid + 1][xcnt + 0] + Std1ImgRegion[tid + 2][xcnt + 0] + Std1ImgRegion[tid + 3][xcnt + 0] + Std1ImgRegion[tid + 4][xcnt + 0] + Std1ImgRegion[tid + 5][xcnt + 0])/12;
			
				ProcLine1 =(Std1ImgRegion[tid + 0][xcnt + 0] + Std1ImgRegion[tid + 0][xcnt + 1] + Std1ImgRegion[tid + 0][xcnt + 2] + Std1ImgRegion[tid + 0][xcnt + 3] + Std1ImgRegion[tid + 0][xcnt + 4] + Std1ImgRegion[tid + 0][xcnt + 5] + Std1ImgRegion[tid + 0][xcnt + 6] +
							Std1ImgRegion[tid + 1][xcnt + 6] + Std1ImgRegion[tid + 2][xcnt + 6] + Std1ImgRegion[tid + 3][xcnt + 6] + Std1ImgRegion[tid + 4][xcnt + 6] + Std1ImgRegion[tid + 5][xcnt + 6])/12;

				ProcLine2 =(Std1ImgRegion[tid + 6][xcnt + 0] + Std1ImgRegion[tid + 6][xcnt + 1] + Std1ImgRegion[tid + 6][xcnt + 2] + Std1ImgRegion[tid + 6][xcnt + 3] + Std1ImgRegion[tid + 6][xcnt + 4] + Std1ImgRegion[tid + 6][xcnt + 5] + Std1ImgRegion[tid + 6][xcnt + 6] +
							Std1ImgRegion[tid + 1][xcnt + 0] + Std1ImgRegion[tid + 2][xcnt + 0] + Std1ImgRegion[tid + 3][xcnt + 0] + Std1ImgRegion[tid + 4][xcnt + 0] + Std1ImgRegion[tid + 5][xcnt + 0])/12;

				ProcLine3 =(Std1ImgRegion[tid + 6][xcnt + 0] + Std1ImgRegion[tid + 6][xcnt + 1] + Std1ImgRegion[tid + 6][xcnt + 2] + Std1ImgRegion[tid + 6][xcnt + 3] + Std1ImgRegion[tid + 6][xcnt + 4] + Std1ImgRegion[tid + 6][xcnt + 5] + Std1ImgRegion[tid + 6][xcnt + 6] +
							Std1ImgRegion[tid + 1][xcnt + 6] + Std1ImgRegion[tid + 2][xcnt + 6] + Std1ImgRegion[tid + 3][xcnt + 6] + Std1ImgRegion[tid + 4][xcnt + 6] + Std1ImgRegion[tid + 5][xcnt + 6])/12;
				
				ProcLine4 =(Std1ImgRegion[tid + 2][xcnt + 2] + Std1ImgRegion[tid + 2][xcnt + 3] + Std1ImgRegion[tid + 2][xcnt + 4] +
							Std1ImgRegion[tid + 3][xcnt + 2] + Std1ImgRegion[tid + 3][xcnt + 3] + Std1ImgRegion[tid + 3][xcnt + 4] +
							Std1ImgRegion[tid + 4][xcnt + 2] + Std1ImgRegion[tid + 4][xcnt + 3] + Std1ImgRegion[tid + 4][xcnt + 4])/9;

				if(ProcLine0 > ProcLine1) ProcLine0 = ProcLine1;
				if(ProcLine2 > ProcLine3) ProcLine2 = ProcLine3;
				if(ProcLine0 > ProcLine2) ProcLine0 = ProcLine2;
				if(ProcLine0 > ProcLine4) ProcLine0 = ProcLine4;

				d_StdImg2[cury*ImageWidth + curx] = ProcLine0*1.2f;
			}
		}
	}	
}

__global__ void gpuSubregionFinding(unsigned short * d_AnnuImg, unsigned short * d_StdImg2, unsigned short * d_PointImg, int * d_RegionNum, int * d_RegionMem, int ImageWidth, int ImageHigh)
{
	// 3x3 processing
	__shared__ unsigned short AnnuImgRegion[ThreadsPerBlock + 2][LineBlockNum_3x3]; //
	__shared__ unsigned short StdImgRegion[ThreadsPerBlock + 2][LineBlockNum_3x3]; //

	int gid = threadIdx.x + blockDim.x*blockIdx.x;
	int tid = threadIdx.x;
	int yoffset = blockDim.x*blockIdx.x;

	int CurStd;

	int curx;
	int cury;
	int xoffset;

	int LineProcNum = (ImageWidth + LineBlockNum - 1) / LineBlockNum;
	int lcnt = 0;
	int xcnt = 0;


	int sum8, sum4; // pixel value sum of center pixels
	int JudgeBit[11]; // judgement memory
	int SumJudge;

	int curRegionPos=0;

	for (lcnt = 0; lcnt < LineProcNum; lcnt++)
	{
		xoffset = lcnt*LineBlockNum;

		// read global image to shared image
		curx = xoffset + tid - 1;
		for (xcnt = 0; xcnt < ThreadsPerBlock + 2; xcnt++)
		{
			// load 32 pixel
			cury = yoffset + xcnt - 1;

			if ((curx >= 0) && (curx < ImageWidth) && (cury >= 0) && (cury < ImageHigh))
			{
				AnnuImgRegion[xcnt][tid] = d_AnnuImg[cury*ImageWidth + curx];
				StdImgRegion[xcnt][tid] = d_StdImg2[cury*ImageWidth + curx];
			}
			else
			{
				AnnuImgRegion[xcnt][tid] = 0;
				StdImgRegion[xcnt][tid]  = 0;
			}
		}
		curx = xoffset + ThreadsPerBlock + tid - 1;
		for (xcnt = 0; xcnt < LineBlockNum + 2; xcnt++)
		{
			// load extra 2 margin pixel
			cury = yoffset + xcnt - 1;
			if(tid < 2)
			{
				if ((curx < ImageWidth) && (cury >= 0) && (cury < ImageHigh)) // (curx >= 0) && 
				{
					AnnuImgRegion[xcnt][tid+ThreadsPerBlock] = d_AnnuImg[cury*ImageWidth + curx];
					StdImgRegion[xcnt][tid+ThreadsPerBlock]  = d_StdImg2[cury*ImageWidth + curx];
				}
				else
				{
					AnnuImgRegion[xcnt][tid+ThreadsPerBlock] = 0;
					StdImgRegion[xcnt][tid+ThreadsPerBlock]  = 0;
				}
			}
		}
		__syncthreads();


		// image filtering

		cury = gid;
		for (xcnt = 0; xcnt < LineBlockNum; xcnt++)
		{
			curx = xoffset + xcnt;

			if ((curx < ImageWidth) && (cury < ImageHigh))
			{

				// average std
				/*
				CurStd = StdImgRegion[tid + 1][xcnt + 1];
			*/				
				
				CurStd = (	StdImgRegion[tid+0][xcnt+0] + StdImgRegion[tid+0][xcnt+1]	+ StdImgRegion[tid+0][xcnt+2] +
							StdImgRegion[tid+1][xcnt+0]									+ StdImgRegion[tid+1][xcnt+2] +
							StdImgRegion[tid+2][xcnt+0] + StdImgRegion[tid+2][xcnt+1]	+ StdImgRegion[tid+2][xcnt+2])/8;
	
				// + StdImgRegion[tid+1][xcnt+1] 



				// sub region judgement

			
				sum8 =	AnnuImgRegion[tid+0][xcnt+0] + AnnuImgRegion[tid+0][xcnt+1] + AnnuImgRegion[tid+0][xcnt+2] +
						AnnuImgRegion[tid+1][xcnt+0] + AnnuImgRegion[tid+1][xcnt+1] + AnnuImgRegion[tid+1][xcnt+2] +
						AnnuImgRegion[tid+2][xcnt+0] + AnnuImgRegion[tid+2][xcnt+1] + AnnuImgRegion[tid+2][xcnt+2];

				sum4 =	AnnuImgRegion[tid+0][xcnt+1] +
						AnnuImgRegion[tid+1][xcnt+0] + AnnuImgRegion[tid+1][xcnt+1] + AnnuImgRegion[tid+1][xcnt+2] +
						AnnuImgRegion[tid+2][xcnt+1] ;

				// subregion judgement
				JudgeBit[0] = AnnuImgRegion[tid + 1][xcnt + 1] > AnnuImgRegion[tid + 0][xcnt + 0];
				JudgeBit[1] = AnnuImgRegion[tid + 1][xcnt + 1] >=AnnuImgRegion[tid + 1][xcnt + 0];
				JudgeBit[2] = AnnuImgRegion[tid + 1][xcnt + 1] > AnnuImgRegion[tid + 2][xcnt + 0];

				JudgeBit[3] = AnnuImgRegion[tid + 1][xcnt + 1] > AnnuImgRegion[tid + 0][xcnt + 1];
				JudgeBit[4] = AnnuImgRegion[tid + 1][xcnt + 1] >=AnnuImgRegion[tid + 2][xcnt + 1];

				JudgeBit[5] = AnnuImgRegion[tid + 1][xcnt + 1] > AnnuImgRegion[tid + 0][xcnt + 2];
				JudgeBit[6] = AnnuImgRegion[tid + 1][xcnt + 1] > AnnuImgRegion[tid + 1][xcnt + 2];
				JudgeBit[7] = AnnuImgRegion[tid + 1][xcnt + 1] > AnnuImgRegion[tid + 2][xcnt + 2];

				JudgeBit[8] = AnnuImgRegion[tid + 1][xcnt + 1] >= 2.2f*CurStd;
				JudgeBit[9]  = sum8 >= 12*CurStd;
				JudgeBit[10] = sum4 >= 9*CurStd;

				SumJudge = JudgeBit[0] + JudgeBit[1] + JudgeBit[2] + JudgeBit[3] + JudgeBit[4] + JudgeBit[5] + JudgeBit[6] + JudgeBit[7] + JudgeBit[8] + JudgeBit[9] + JudgeBit[10];
				
				if(SumJudge==11)
				{
					if ((curx >= 4) && (cury >= 4) && (curx < ImageWidth - 4) && (cury < ImageHigh - 4))
					{
//						printf("cur std %d %d :%d-%d-%d\n", curx, cury, CurStd, sum8 / CurStd, sum4 / CurStd);

						d_PointImg[cury*ImageWidth + curx] = 1; // there is a valid subregion here

						curRegionPos = atomicAdd(d_RegionNum, 1);
						d_RegionMem[curRegionPos * 3 + 0] = curx;
						d_RegionMem[curRegionPos * 3 + 1] = cury;
						d_RegionMem[curRegionPos * 3 + 2] = CurStd;

					}
				}
			}
		}
	}
}

char IsValidClusterPos(unsigned short PosArry[49]);
void GetClusterPos(unsigned short PosArry[49], int *xoffset, int *yoffset);


// CPU program
void cpuPointCluster(unsigned short * h_PointImg, int *h_iPosCount, int *h_iPosMem, unsigned short * h_cPointImg, int *h_oPosCount, int *h_oPosMem, int ImageWidth, int ImageHigh)
{
	// the h_PointImg and h_cPointImg are the same size with ImageWidth x ImageHigh
	// the h_iPosMem and h_iPosCount indicate original subregion pozition and number
	//  the h_oPosMem and h_oPosCount indicate clustered subregion pozition and number
	// for each image frame, consider maximum 2048x2048 image size and maximum point density, 

	memset(h_cPointImg, 0, ImageWidth*ImageHigh*sizeof(unsigned short));
	*h_oPosCount = 0;

	int xpos, ypos; // original x,y position
	int cxpos, cypos; // clustered x,y position
	int CurStd;

	int tempx, tempy;
	int AddrOffset;

	int pcnt = 0;
	int tcnt = 0;

	unsigned short PointArry[49];
	unsigned short ClusterPointArry[49];

	for (pcnt = 0; pcnt < (*h_iPosCount); pcnt++)
	{
		// for each point
		xpos = h_iPosMem[pcnt * 3 + 0];
		ypos = h_iPosMem[pcnt * 3 + 1];
		CurStd = h_iPosMem[pcnt * 3 + 2];

		tempx = xpos - 3;
		tempy = ypos - 3;

		// load image region around cur point
		for (tcnt = 0; tcnt<7; tcnt++)
		{
			AddrOffset = (tempy + tcnt)*ImageWidth;

			PointArry[tcnt * 7 + 0] = h_PointImg[AddrOffset + tempx + 0];
			PointArry[tcnt * 7 + 1] = h_PointImg[AddrOffset + tempx + 1];
			PointArry[tcnt * 7 + 2] = h_PointImg[AddrOffset + tempx + 2];
			PointArry[tcnt * 7 + 3] = h_PointImg[AddrOffset + tempx + 3];
			PointArry[tcnt * 7 + 4] = h_PointImg[AddrOffset + tempx + 4];
			PointArry[tcnt * 7 + 5] = h_PointImg[AddrOffset + tempx + 5];
			PointArry[tcnt * 7 + 6] = h_PointImg[AddrOffset + tempx + 6];

			ClusterPointArry[tcnt * 7 + 0] = h_cPointImg[AddrOffset + tempx + 0];
			ClusterPointArry[tcnt * 7 + 1] = h_cPointImg[AddrOffset + tempx + 1];
			ClusterPointArry[tcnt * 7 + 2] = h_cPointImg[AddrOffset + tempx + 2];
			ClusterPointArry[tcnt * 7 + 3] = h_cPointImg[AddrOffset + tempx + 3];
			ClusterPointArry[tcnt * 7 + 4] = h_cPointImg[AddrOffset + tempx + 4];
			ClusterPointArry[tcnt * 7 + 5] = h_cPointImg[AddrOffset + tempx + 5];
			ClusterPointArry[tcnt * 7 + 6] = h_cPointImg[AddrOffset + tempx + 6];
		}

		if (IsValidClusterPos(ClusterPointArry))
		{
			// there is no extraction point around raw pos, else cur point is ignored

			GetClusterPos(PointArry, &cxpos, &cypos);
			cxpos += xpos;
			cypos += ypos;

			tempx = cxpos - 3;
			tempy = cypos - 3;

			// load image region around cur point
			for (tcnt = 0; tcnt<7; tcnt++)
			{
				AddrOffset = (tempy + tcnt)*ImageWidth;

				ClusterPointArry[tcnt * 7 + 0] = h_cPointImg[AddrOffset + tempx + 0];
				ClusterPointArry[tcnt * 7 + 1] = h_cPointImg[AddrOffset + tempx + 1];
				ClusterPointArry[tcnt * 7 + 2] = h_cPointImg[AddrOffset + tempx + 2];
				ClusterPointArry[tcnt * 7 + 3] = h_cPointImg[AddrOffset + tempx + 3];
				ClusterPointArry[tcnt * 7 + 4] = h_cPointImg[AddrOffset + tempx + 4];
				ClusterPointArry[tcnt * 7 + 5] = h_cPointImg[AddrOffset + tempx + 5];
				ClusterPointArry[tcnt * 7 + 6] = h_cPointImg[AddrOffset + tempx + 6];
			}

			if (IsValidClusterPos(ClusterPointArry))
			{
				// accept clustered point pos
				h_cPointImg[cypos*ImageWidth + cxpos] = 1;

				h_oPosMem[3 * (*h_oPosCount) + 0] = cxpos;
				h_oPosMem[3 * (*h_oPosCount) + 1] = cypos;
				h_oPosMem[3 * (*h_oPosCount) + 2] = CurStd;

				*h_oPosCount = *h_oPosCount + 1;
			}
			else
			{
				// accept original point pos
				h_cPointImg[ypos*ImageWidth + xpos] = 1;

				h_oPosMem[3 * (*h_oPosCount) + 0] = xpos;
				h_oPosMem[3 * (*h_oPosCount) + 1] = ypos;
				h_oPosMem[3 * (*h_oPosCount) + 2] = CurStd;
				*h_oPosCount = *h_oPosCount + 1;
			}
		}
	}
}

// CPU program
void GetClusterPos(unsigned short PosArry[49], int *xoffset, int *yoffset)
{
	// get cluster pos by center of mass
	unsigned short(*pArry)[7] = (unsigned short(*)[7])PosArry;

	int sumdat = 0;
	int cnt = 0;
	int xdat[7] = { 0, 0, 0, 0, 0, 0, 0 };
	int ydat[7] = { 0, 0, 0, 0, 0, 0, 0 };

	float xresult;
	float yresult;


	for (cnt = 0; cnt<7; cnt++)
	{
		xdat[0] += pArry[cnt][0];
		xdat[1] += pArry[cnt][1];
		xdat[2] += pArry[cnt][2];
		xdat[3] += pArry[cnt][3];
		xdat[4] += pArry[cnt][4];
		xdat[5] += pArry[cnt][5];
		xdat[6] += pArry[cnt][6];

		ydat[0] += pArry[0][cnt];
		ydat[1] += pArry[1][cnt];
		ydat[2] += pArry[2][cnt];
		ydat[3] += pArry[3][cnt];
		ydat[4] += pArry[4][cnt];
		ydat[5] += pArry[5][cnt];
		ydat[6] += pArry[6][cnt];

		sumdat += (pArry[cnt][0] + pArry[cnt][1] + pArry[cnt][2] + pArry[cnt][3] + pArry[cnt][4] + pArry[cnt][5] + pArry[cnt][6]);
	}
	// x,y cluster pos
	xresult = xdat[1] * 1 + xdat[2] * 2 + xdat[3] * 3 + xdat[4] * 4 + xdat[5] * 5 + xdat[6] * 6;
	yresult = ydat[1] * 1 + ydat[2] * 2 + ydat[3] * 3 + ydat[4] * 4 + ydat[5] * 5 + ydat[6] * 6;
	xresult = xresult / sumdat + 0.5f;
	yresult = yresult / sumdat + 0.5f;

	*xoffset = xresult - 3;
	*yoffset = yresult - 3;

}


// CPU program
char IsValidClusterPos(unsigned short PosArry[49])
{
	unsigned short(*pArry)[7] = (unsigned short(*)[7])PosArry;
	int cnt = 0;
	int sumdat = 0;
	char oval = 0;
	for (cnt = 0; cnt<7; cnt++)
	{
		sumdat += (pArry[cnt][0] + pArry[cnt][1] + pArry[cnt][2] + pArry[cnt][3] + pArry[cnt][4] + pArry[cnt][5] + pArry[cnt][6]);
	}
	if (sumdat)oval = 0;
	else oval = 1;

	return oval;
}


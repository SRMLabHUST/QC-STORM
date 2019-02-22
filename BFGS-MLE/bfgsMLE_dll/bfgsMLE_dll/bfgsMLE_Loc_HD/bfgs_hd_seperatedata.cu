#include "bfgs_hd_seperatedata.h"


#define RegionLen			9
#define SubRegionWidth		(RegionLen*RegionLen)	// 7*7

#define WholeImageWidth	(RegionLen*(RegionLen+1)) // 56 = 7*8, can gen bank conflict
// avoid bank conflict, if no conflict, +1 can be avoid
#define SharedImageWidth1	(WholeImageWidth+1)

__device__  char IsPoint(unsigned short* CurRegion, int stdth);


__global__ void GPUSubregionSeperate(unsigned short *d_SubRegion, int FluoNum, int *d_Point1Num, int *d_Point1PosArry, int *d_Point2Num, int *d_Point2PosArry, int *d_Point3Num, int *d_Point3PosArry)
{
	// raw image
	__shared__ float ImageRegion[ThreadsPerBlock][SharedImageWidth1];
	// filtered by a high pass filter to judge point number
	__shared__ float FillerdRegion[ThreadsPerBlock][SubRegionWidth];


	int gid =	threadIdx.x+blockDim.x*blockIdx.x;
	int tid =	threadIdx.x;

	if(gid>=FluoNum) return;
	
	int gMemPos = gid*WholeImageWidth;
	
	float (*pRawImg)[RegionLen]=(float(*)[RegionLen])&ImageRegion[tid][0];
	float (*pFillerdImg)[RegionLen]=(float(*)[RegionLen])&FillerdRegion[tid][0];

	unsigned short CurRegion[9];
	int CurStd;

	float tempdat;

	int cnt;
	int row,col;

	for(cnt=0; cnt<WholeImageWidth; cnt++)
	{
		// read sub image into shared memory
		ImageRegion[tid][cnt] = d_SubRegion[gMemPos + cnt];
	}

	CurStd = ImageRegion[tid][89];

	pFillerdImg[0][0]=0; pFillerdImg[0][1]=0; pFillerdImg[0][2]=0; pFillerdImg[0][3]=0; pFillerdImg[0][4]=0; pFillerdImg[0][5]=0; pFillerdImg[0][6]=0; pFillerdImg[0][7]=0; pFillerdImg[0][8]=0;
	pFillerdImg[8][0]=0; pFillerdImg[8][1]=0; pFillerdImg[8][2]=0; pFillerdImg[8][3]=0; pFillerdImg[8][4]=0; pFillerdImg[8][5]=0; pFillerdImg[8][6]=0; pFillerdImg[8][7]=0; pFillerdImg[8][8]=0;
	pFillerdImg[1][0]=0; pFillerdImg[2][0]=0; pFillerdImg[3][0]=0; pFillerdImg[4][0]=0; pFillerdImg[5][0]=0; pFillerdImg[6][0]=0; pFillerdImg[7][0]=0;
	pFillerdImg[1][8]=0; pFillerdImg[2][8]=0; pFillerdImg[3][8]=0; pFillerdImg[4][8]=0; pFillerdImg[5][8]=0; pFillerdImg[6][8]=0; pFillerdImg[7][8]=0;


	// image filter
/*
	// f2=[
    -0.1     -0.15   -0.1
    -0.15    1       -0.15
    -0.1     -0.15    -0.1];
*/
	for(row=1; row < RegionLen-1; row++)
	{
		for(col=1; col < RegionLen-1; col++)
		{
			tempdat={
				0 - pRawImg[row-1][col-1]*0.10f - pRawImg[row-1][col  ]*0.15f - pRawImg[row-1][col+1]*0.10f - 
					pRawImg[row  ][col-1]*0.15f + pRawImg[row  ][col  ]       - pRawImg[row  ][col+1]*0.15f - 
					pRawImg[row+1][col-1]*0.10f - pRawImg[row+1][col  ]*0.15f - pRawImg[row+1][col+1]*0.10f};

			if (tempdat < 0)tempdat = 0;
			pFillerdImg[row][col] = tempdat;
		}
	}

  	int PointNum=0;
  	int curpos;

	for(row=1; row<RegionLen-1; row++)
	{
		for(col=1; col<RegionLen-1; col++)
		{
			CurRegion[0] = pFillerdImg[row-1][col-1]; CurRegion[1] = pFillerdImg[row-1][col  ];  CurRegion[2] = pFillerdImg[row-1][col+1];
			CurRegion[3] = pFillerdImg[row  ][col-1]; CurRegion[4] = pFillerdImg[row  ][col  ];  CurRegion[5] = pFillerdImg[row  ][col+1];
			CurRegion[6] = pFillerdImg[row+1][col-1]; CurRegion[7] = pFillerdImg[row+1][col  ];  CurRegion[8] = pFillerdImg[row+1][col+1];

			if (IsPoint(CurRegion, CurStd))
			{
				PointNum++;
				if(PointNum<5)
				{
					pRawImg[9][2+PointNum] = col*256 + row; // high 8 bit and low 8 bit
				}
			}
		}
	}

	if(PointNum==1)
	{
		curpos = atomicAdd(d_Point1Num, 1);
		d_Point1PosArry[curpos] = gid;
	}
	else if(PointNum==2)
	{
		curpos = atomicAdd(d_Point2Num, 1);
		d_Point2PosArry[curpos] = gid;
	}
	else if(PointNum==3)
	{
		curpos = atomicAdd(d_Point3Num, 1);
		d_Point3PosArry[curpos] = gid;
	}



	for(cnt=SubRegionWidth; cnt<WholeImageWidth; cnt++)
	{
		// write useful information back to global memory
		d_SubRegion[gMemPos + cnt] = (unsigned short)ImageRegion[tid][cnt];
	}


}

__device__  char IsPoint(unsigned short* CurRegion, int stdth)
{
	if (stdth < 9)stdth = 9;

	char PointVal=0;
	unsigned short (*pRegion)[3]=(unsigned short(*)[3])&CurRegion[0];
	int Judge=0;

	int sum9;
	int sum5;
	int cdat;

	sum9 =	pRegion[0][0] + pRegion[0][1] + pRegion[0][2] +
			pRegion[1][0] + pRegion[1][1] + pRegion[1][2] +
			pRegion[2][0] + pRegion[2][1] + pRegion[2][2];

	sum5 =	pRegion[0][1] +
			pRegion[1][0] + pRegion[1][1] + pRegion[1][2] +
			pRegion[2][1];

	cdat = pRegion[1][1];

	Judge+= (pRegion[1][1] > pRegion[0][0]);
	Judge+= (pRegion[1][1] >=pRegion[0][1]);
	Judge+= (pRegion[1][1] > pRegion[0][2]);

	Judge+= (pRegion[1][1] > pRegion[1][0]);
	Judge+= (pRegion[1][1] >=pRegion[1][2]);

	Judge+= (pRegion[1][1] > pRegion[2][0]);
	Judge+= (pRegion[1][1] > pRegion[2][1]);
	Judge+= (pRegion[1][1] > pRegion[2][2]);

	Judge += (sum9 > stdth * 8);
	Judge += (sum5 > stdth * 6);
	Judge += (cdat > stdth * 2);



	if (Judge == 11)PointVal = 1;
	else PointVal = 0;

	return PointVal;
	
}


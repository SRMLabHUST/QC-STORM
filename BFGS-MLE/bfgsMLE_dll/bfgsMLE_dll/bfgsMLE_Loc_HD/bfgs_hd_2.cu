#include "bfgs_hd_2.h"



__device__ static void MatMulVector_2p(float D0[], float grad[][ThreadsPerBlock], float d0[][ThreadsPerBlock], int tid);
__device__ static void ConstructD0_2p(float D0[], float sk[], float yk[], int tid);
__device__ static void IninfConstrain_2p(float Ininf[][ThreadsPerBlock], int tid);
__device__ static void VectorAddMul1_2p(float oVector[], float Ininf[][ThreadsPerBlock], float d0[][ThreadsPerBlock], float coeff, int tid);

#define OutParaNum			8 // a,x,y,f,a,x,y,f

#define FitParaNum 			7
#define FittingNum			2

#define IterateNum			8  //8
#define IterateNum_bs		11 //11


// for region size			7

#define RegionLen			9
#define SubRegionWidth		(RegionLen*RegionLen)	// 7*7

#define WholeImageWidth	(RegionLen*(RegionLen+1)) // 56 = 7*8, can gen bank conflict
// avoid bank conflict, if no conflict, +1 can be avoid
#define SharedImageWidth1	(WholeImageWidth+1)


__device__ static void MLELocalization_2p(float subimg[][SharedImageWidth1], float Ininf[][ThreadsPerBlock], float grad[][ThreadsPerBlock], float d0[][ThreadsPerBlock], float D0[], float PsfWidth, int tid);

__device__ static void PreLocolization_2p(float subimg[][SharedImageWidth1], float Ininf[][ThreadsPerBlock], int tid);

__device__ static void poissonfGradient_2p(float subimg[][SharedImageWidth1], float Ininf[][ThreadsPerBlock], float grad[][ThreadsPerBlock], float PsfWidth, int tid);
__device__ static float poissonf_2p(float subimg[][SharedImageWidth1], float Ininf[], float PsfWidth, int tid);

// algorithm core codes

__global__ void MLEROILocTop_2p(unsigned short *d_SubRegion, float *d_LocArry,int *d_Point2PosArry, float Offset, float kadc, float PsfWidth, int FluoNum, int isFPGAProc)
{
	__shared__ float ImageRegion[ThreadsPerBlock][SharedImageWidth1];
	__shared__ float D0[ThreadsPerBlock][FitParaNum*FitParaNum];	// inv of matrix hessian 

	// avoid bank conflict
	__shared__ float Ininf[FitParaNum][ThreadsPerBlock]; // cur position
	__shared__ float grad[FitParaNum][ThreadsPerBlock];  // gradient
	__shared__ float d0[FitParaNum][ThreadsPerBlock];	// direction


	int gid =	threadIdx.x+blockDim.x*blockIdx.x;
	int tid =	threadIdx.x;
	
	if(gid>=FluoNum) return;

	// curid = d_Point2PosArry[gid]
	int BlockOffset = blockDim.x*blockIdx.x;
	int gMemPos;

	// read sub image into shared memory

	int cnt=0;
	int icnt = 0;

	int SensorId;
	float ypos;
	float (*pD0)[FitParaNum]=(float(*)[FitParaNum])&D0[tid][0];

	float MaxAmp=0;
	int PointTh1 = 0;
	int PointTh2 = 0;

	// load image region from global memory, there two method to load

	for (cnt = 0; cnt < ThreadsPerBlock; cnt++)
	{
		gMemPos = (BlockOffset + cnt)*WholeImageWidth;

		ImageRegion[cnt][tid] = (d_SubRegion[gMemPos + tid] - Offset)*kadc;
		ImageRegion[cnt][tid + 32] = (d_SubRegion[gMemPos + tid + 32] - Offset)*kadc;

		if (tid < (SubRegionWidth % 32))
		{
			ImageRegion[cnt][tid + 64] = (d_SubRegion[gMemPos + tid + 64] - Offset)*kadc;
		}
		if (tid < RegionLen)
		{
			ImageRegion[cnt][tid + SubRegionWidth] = d_SubRegion[gMemPos + tid + SubRegionWidth];
		}
	}

	//initial D0
#pragma unroll
	for (cnt = 0; cnt < FitParaNum; cnt++)
	{
#pragma unroll
		for (icnt = 0; icnt < FitParaNum; icnt++)
		{
			pD0[cnt][icnt] = 0.0f;
		}
	}
#pragma unroll
	for (icnt = 0; icnt < FitParaNum; icnt++)
	{
		pD0[icnt][icnt] = 1.0f; // set diagonal to 1
	}

	// why this function impact the time so much?
	__syncthreads();


	PreLocolization_2p(ImageRegion,Ininf,tid);

	poissonfGradient_2p(ImageRegion, Ininf, grad, PsfWidth, tid);

#pragma unroll
	for (cnt = 0; cnt < FitParaNum; cnt++)
	{
		d0[cnt][tid] = -grad[cnt][tid];
	}


	MLELocalization_2p(ImageRegion, Ininf, grad, d0, &D0[tid][0], PsfWidth, tid);
	
	// a, x, y, a, x, y, b
	// 0  1  2  3  4  5  6



	if(isFPGAProc==0)
	{
		// online localization

		// true	for camera link to usb
		Ininf[1][tid] = Ininf[1][tid] + ImageRegion[tid][SubRegionWidth + 1] - 3.0f; // X
		Ininf[2][tid] = Ininf[2][tid] + ImageRegion[tid][SubRegionWidth + 2] - 3.0f; // Y
		d_LocArry[gid*OutParaNum + 4] = ImageRegion[tid][SubRegionWidth + 3]; // frame inf


		d_LocArry[gid*OutParaNum + 0] = Ininf[0][tid] * 256.0f;	// A
		d_LocArry[gid*OutParaNum + 1] = Ininf[1][tid];			// X

		d_LocArry[gid*OutParaNum + 3] = Ininf[3][tid];			// sigma


		SensorId = ((int)ImageRegion[tid][SubRegionWidth + 5]) & 0x03;

		if (SensorId == 1) // up region sensor
		{
			// since only a image belt in the imge sensor center is valid, the region Offset is removed
			// 1030 = 1024 +6, 1018 = 1024 - 6,
			// remove the y pos Offset from image high

			// get true ypos
//			ypos = 1024.0f + 3.0f - Ininf[2][tid];
			ypos = 1024.0f + 6.0f - Ininf[2][tid];
			if (ypos > 1024.0f)
			{
				d_LocArry[gid*OutParaNum + 1] = 0.0f;			// X
				d_LocArry[gid*OutParaNum + 2] = 0.0f;
			}
			else
			{
				d_LocArry[gid*OutParaNum + 2] = ypos - (1024.0f - ImageRegion[tid][SubRegionWidth + 4]); // Y
			}
		}
		else  // bottom region sensor
		{
			// get true ypos
//			ypos = 1024.0f - 3.0f + Ininf[2][tid];
			ypos = 1024.0f - 6.0f + Ininf[2][tid];
			if (ypos < 1024.0f)
			{
				d_LocArry[gid*OutParaNum + 1] = 0.0f;			// X
				d_LocArry[gid*OutParaNum + 2] = 0.0f;
			}
			else
			{
				d_LocArry[gid*OutParaNum + 2] = ypos - (1024.0f - ImageRegion[tid][SubRegionWidth + 4]); // Y
			}
		}
	}else
	{
		// offline localization
		if (Ininf[0][tid]>Ininf[3][tid])MaxAmp = Ininf[0][tid];
		else MaxAmp = Ininf[3][tid];
		PointTh1 = (MaxAmp / Ininf[0][tid]) > 4;
		PointTh2 = (MaxAmp / Ininf[3][tid]) > 4;

		// discard some point
		if ((Ininf[1][tid]>1.5f) && (Ininf[1][tid]<7.5f) && (Ininf[2][tid]>1.5f) && (Ininf[2][tid] < 7.5f) && PointTh1)
		{
			Ininf[1][tid] = Ininf[1][tid] + ImageRegion[tid][SubRegionWidth + 0] - 4.0f; // X
			Ininf[2][tid] = Ininf[2][tid] + ImageRegion[tid][SubRegionWidth + 1] - 4.0f; // Y
		}
		else
		{
			Ininf[1][tid] = 0.0f; // X
			Ininf[2][tid] = 0.0f; // Y
		}
		if ((Ininf[4][tid]>1.5f) && (Ininf[4][tid]<7.5f) && (Ininf[5][tid]>1.5f) && (Ininf[5][tid] < 7.5f) && PointTh2)
		{
			Ininf[4][tid] = Ininf[4][tid] + ImageRegion[tid][SubRegionWidth + 0] - 4.0f; // X
			Ininf[5][tid] = Ininf[5][tid] + ImageRegion[tid][SubRegionWidth + 1] - 4.0f; // Y
		}
		else
		{
			Ininf[4][tid] = 0.0f; // X
			Ininf[5][tid] = 0.0f; // Y
		}

		d_LocArry[gid*OutParaNum + 0] = Ininf[0][tid] * 128.0f;	// A
		d_LocArry[gid*OutParaNum + 1] = Ininf[1][tid];	// X
		d_LocArry[gid*OutParaNum + 2] = Ininf[2][tid];	// Y
		d_LocArry[gid*OutParaNum + 3] = ImageRegion[tid][SubRegionWidth + 2];	// f

		d_LocArry[gid*OutParaNum + 4] = Ininf[3][tid] * 128.0f;	// A
		d_LocArry[gid*OutParaNum + 5] = Ininf[4][tid];	// X
		d_LocArry[gid*OutParaNum + 6] = Ininf[5][tid];	// Y
		d_LocArry[gid*OutParaNum + 7] = ImageRegion[tid][SubRegionWidth + 2];	// f

	}

}



__device__ static void MLELocalization_2p(float subimg[][SharedImageWidth1], float Ininf[][ThreadsPerBlock], float grad[][ThreadsPerBlock], float d0[][ThreadsPerBlock], float D0[], float PsfWidth, int tid)
{
	// adjust d0
	float td0[FitParaNum];
	float td0_total;

	float sk[FitParaNum]; // bfgs quasi-newton method
	float yk[FitParaNum];
	float tgrad[FitParaNum];

	int cnt = 0;

	int itcnt = 0; // iteration number

	// find work length, divide 2 method


	float scale;

	float xd[FitParaNum * 2];

	float ddat[2];
	float dpos[2];

	int xdsel = 0;

	for (itcnt = 0; itcnt< IterateNum; itcnt++)
	{
		// adjust d0
		td0[0] = abs(d0[0][tid]);
		td0[1] = abs(d0[1][tid]);
		td0[2] = abs(d0[2][tid]);
		td0[3] = abs(d0[3][tid]);
		td0[4] = abs(d0[4][tid]);
		td0[5] = abs(d0[5][tid]);
		td0[6] = abs(d0[6][tid]);

		td0_total = __fdividef( 8.4f, (td0[0] + td0[1] + td0[2] + td0[3] + td0[4] + td0[5] + td0[6])); // /7/1.2


		d0[0][tid] = d0[0][tid] * td0_total;
		d0[1][tid] = d0[1][tid] * td0_total;
		d0[2][tid] = d0[2][tid] * td0_total;
		d0[3][tid] = d0[3][tid] * td0_total;
		d0[4][tid] = d0[4][tid] * td0_total;
		d0[5][tid] = d0[5][tid] * td0_total;
		d0[6][tid] = d0[6][tid] * td0_total;


		dpos[0] = 0.0001f; // scale factor left limit, should not equal to 0 and smaller
		dpos[1] = 1.0f; // scale factor right limit, should not lager than 2

		VectorAddMul1_2p(&xd[0], Ininf, d0, 0.0001f, tid);
		VectorAddMul1_2p(&xd[FitParaNum], Ininf, d0, 1.0f, tid);

		ddat[0] = poissonf_2p(subimg, &xd[0], PsfWidth, tid);
		ddat[1] = poissonf_2p(subimg, &xd[FitParaNum], PsfWidth, tid);
		// 
		for (cnt = 0; cnt<IterateNum_bs; cnt++)
		{
			if (ddat[0]<ddat[1])
			{
				xdsel = 1; // right shift
			}
			else
			{
				xdsel = 0; // left shift
			}

			dpos[xdsel] = (dpos[0] + dpos[1])*0.5f; //  /2.0f which one shift

			if (cnt<IterateNum_bs - 1)
			{
				VectorAddMul1_2p(&xd[xdsel*FitParaNum], Ininf, d0, dpos[xdsel], tid);// xd=ininf+d0*scale
				ddat[xdsel] = poissonf_2p(subimg, &xd[xdsel*FitParaNum], PsfWidth, tid);


			}

		}
		scale = (dpos[0] + dpos[1])*0.5f; //  /2.0f calculated direction

		// sk 		= d0*base;
		// inInf 	= inInf+d0*base;
		sk[0] = d0[0][tid] * scale;
		sk[1] = d0[1][tid] * scale;
		sk[2] = d0[2][tid] * scale;
		sk[3] = d0[3][tid] * scale;
		sk[4] = d0[4][tid] * scale;
		sk[5] = d0[5][tid] * scale;
		sk[6] = d0[6][tid] * scale;

		Ininf[0][tid] = Ininf[0][tid] + sk[0];
		Ininf[1][tid] = Ininf[1][tid] + sk[1];
		Ininf[2][tid] = Ininf[2][tid] + sk[2];
		Ininf[3][tid] = Ininf[3][tid] + sk[3];
		Ininf[4][tid] = Ininf[4][tid] + sk[4];
		Ininf[5][tid] = Ininf[5][tid] + sk[5];
		Ininf[6][tid] = Ininf[6][tid] + sk[6];


		IninfConstrain_2p(Ininf, tid); // 

		if (itcnt<IterateNum - 1)
		{
			// tgrad = grad;
			tgrad[0] = grad[0][tid];
			tgrad[1] = grad[1][tid];
			tgrad[2] = grad[2][tid];
			tgrad[3] = grad[3][tid];
			tgrad[4] = grad[4][tid];
			tgrad[5] = grad[5][tid];
			tgrad[6] = grad[6][tid];

			poissonfGradient_2p(subimg, Ininf, grad, PsfWidth, tid);

			yk[0] = grad[0][tid] - tgrad[0];
			yk[1] = grad[1][tid] - tgrad[1];
			yk[2] = grad[2][tid] - tgrad[2];
			yk[3] = grad[3][tid] - tgrad[3];
			yk[4] = grad[4][tid] - tgrad[4];
			yk[5] = grad[5][tid] - tgrad[5];
			yk[6] = grad[6][tid] - tgrad[6];



			ConstructD0_2p(D0, sk, yk, tid);
			MatMulVector_2p(D0, grad, d0, tid);
		}
	}
}


__device__ static float poissonf_2p(float subimg[][SharedImageWidth1], float Ininf[], float PsfWidth, int tid)
{
	int row, col;

	float rowpos, colpos;
	float VPoss = 0;
	float IVal=0;
	float pixval=0;
	float tdat=0;

	float(*pSubImg)[RegionLen] = (float(*)[RegionLen])&subimg[tid][0];

	if (Ininf[0]<0.0f) Ininf[0] = 0.0f; // A
	if (Ininf[1]<0.0f) Ininf[1] = 0.0f; // X
	if (Ininf[2]<0.0f) Ininf[2] = 0.0f; // Y
	if (Ininf[3]<0.0f) Ininf[3] = 0.0f; // A
	if (Ininf[4]<0.0f) Ininf[4] = 0.0f; // X
	if (Ininf[5]<0.0f) Ininf[5] = 0.0f; // Y
	if (Ininf[6]<0.0f) Ininf[6] = 0.0f; // bg

	// a, x, y, a, x, y, b
	// 0  1  2  3  4  5  6

	for (row = 0; row<RegionLen; row++)
	{
		for (col = 0; col<RegionLen; col++)
		{
			rowpos = row + 0.5f; // -0.5
			colpos = col + 0.5f; // -0.5

			// get model intensity
			tdat = ((colpos - Ininf[1])*(colpos - Ininf[1]) + (rowpos - Ininf[2])*(rowpos - Ininf[2]))*PsfWidth;
			IVal = 128.0f*Ininf[0] * __expf(-tdat) + Ininf[6] * 32.0f; // point 1

			tdat = ((colpos - Ininf[4])*(colpos - Ininf[4]) + (rowpos - Ininf[5])*(rowpos - Ininf[5]))*PsfWidth;
			IVal = IVal+ 128.0f*Ininf[3] * __expf(-tdat); // point 2

			// pf goal function
			pixval = pSubImg[row][col]; //subimg[tid][row*RegionLen + col]; // coloffset+

			//			tdat = IVal - pixval - pixval*(__logf(IVal) - __logf(pixval + 80.12345f)); //to avoid negative value in logf function, that will be a disaster
			tdat = IVal - pixval - pixval*(__logf(IVal) - __logf(fmaxf(pixval, 1.0f))); //to avoid negative value in logf function, that will be a disaster

			VPoss = VPoss + IVal;
		}
	}
	return VPoss;
}

__device__ static void poissonfGradient_2p(float subimg[][SharedImageWidth1], float Ininf[][ThreadsPerBlock], float grad[][ThreadsPerBlock], float PsfWidth, int tid)
{
	float tIninf[FitParaNum];
	float tgradn;
	float tgradp;
	int cnt;
	//lv 0

	for (cnt = 0; cnt<FitParaNum; cnt++)
	{
		tIninf[cnt] = Ininf[cnt][tid];
	}

	for (cnt = 0; cnt<FitParaNum; cnt++)
	{
		tIninf[cnt] = tIninf[cnt] - 0.002f;
		tgradn = poissonf_2p(subimg, tIninf, PsfWidth, tid);

		tIninf[cnt] = tIninf[cnt] + 0.004f;
		tgradp = poissonf_2p(subimg, tIninf, PsfWidth, tid);

		tIninf[cnt] = tIninf[cnt] - 0.002f;

		grad[cnt][tid] = (tgradp - tgradn)*250.0f; // /0.004f;
	}
}


__device__ static void PreLocolization_2p(float subimg[][SharedImageWidth1], float Ininf[][ThreadsPerBlock], int tid)
{
	//	int rowoff=tid*8;

	float sum_bg = 0;

	unsigned int curx, cury;

	float(*pSubImg)[RegionLen] = (float(*)[RegionLen])&subimg[tid][0];

	sum_bg = 	pSubImg[0][0] + pSubImg[0][1] + pSubImg[0][2] + pSubImg[0][3] + pSubImg[0][4] + pSubImg[0][5] + pSubImg[0][6] + pSubImg[0][7] + pSubImg[0][8] +\
				pSubImg[8][0] + pSubImg[8][1] + pSubImg[8][2] + pSubImg[8][3] + pSubImg[8][4] + pSubImg[8][5] + pSubImg[8][6] + pSubImg[8][7] + pSubImg[8][8] +\
				pSubImg[1][0] + pSubImg[2][0] + pSubImg[3][0] + pSubImg[4][0] + pSubImg[5][0] + pSubImg[6][0] + pSubImg[7][0] +\
				pSubImg[1][8] + pSubImg[2][8] + pSubImg[3][8] + pSubImg[4][8] + pSubImg[5][8] + pSubImg[6][8] + pSubImg[7][8];

	sum_bg = sum_bg / 32.0f;// *0.041667f; // /24.0f; // background

	curx = (((int)pSubImg[9][3]) / 256);
	cury = (((int)pSubImg[9][3]) % 256);

	Ininf[0][tid] = (pSubImg[cury][curx] - sum_bg)/128.0f; // /256.0f;	// Amplitude
	Ininf[1][tid] = curx + 0.5f; // sumDatx / sumdat - 0.5f;	// x
	Ininf[2][tid] = cury + 0.5f; // sumDaty / sumdat - 0.5f;	// y

	curx = (((int)pSubImg[9][4]) / 256);
	cury = (((int)pSubImg[9][4]) % 256);

	Ininf[3][tid] = (pSubImg[cury][curx] - sum_bg)/128.0f; // /256.0f;	// Amplitude
	Ininf[4][tid] = curx + 0.5f; // sumDatx / sumdat - 0.5f;	// x
	Ininf[5][tid] = cury + 0.5f; // sumDaty / sumdat - 0.5f;	// y

	Ininf[6][tid] = sum_bg/32.0f; // /128.0f;	// bg

}


// bfgs core algorithm
// 7 parameters


__device__ static void VectorAddMul1_2p(float oVector[], float Ininf[][ThreadsPerBlock], float d0[][ThreadsPerBlock], float coeff, int tid)
{
	oVector[0] = d0[0][tid] * coeff + Ininf[0][tid];
	oVector[1] = d0[1][tid] * coeff + Ininf[1][tid];
	oVector[2] = d0[2][tid] * coeff + Ininf[2][tid];
	oVector[3] = d0[3][tid] * coeff + Ininf[3][tid];
	oVector[4] = d0[4][tid] * coeff + Ininf[4][tid];
	oVector[5] = d0[5][tid] * coeff + Ininf[5][tid];
	oVector[6] = d0[6][tid] * coeff + Ininf[6][tid];

}

__device__ static void IninfConstrain_2p(float Ininf[][ThreadsPerBlock], int tid)
{
	if (Ininf[0][tid]<0.001f) Ininf[0][tid] = 0.001f; // A
	if (Ininf[1][tid]<0.001f) Ininf[1][tid] = 0.001f; // X
	if (Ininf[2][tid]<0.001f) Ininf[2][tid] = 0.001f; // Y
	if (Ininf[3][tid]<0.001f) Ininf[3][tid] = 0.001f; // sigma
	if (Ininf[4][tid]<0.001f) Ininf[4][tid] = 0.001f; // bg
	if (Ininf[5][tid]<0.001f) Ininf[5][tid] = 0.001f; // sigma
	if (Ininf[6][tid]<0.001f) Ininf[6][tid] = 0.001f; // bg

}

__device__ static void ConstructD0_2p(float D0[], float sk[FitParaNum], float yk[FitParaNum], int tid)
{
	float divdat;
	float skyk[FitParaNum*FitParaNum]; //  I-(sk*yk')/(yk'*sk)
	//	float yksk[25]; // the same with skyk but transposition
	float sksk[FitParaNum*FitParaNum];
	float tD0[FitParaNum*FitParaNum];

	float(*pskyk)[FitParaNum] = (float(*)[FitParaNum])skyk;
	float(*psksk)[FitParaNum] = (float(*)[FitParaNum])sksk;
	float(*ptD0)[FitParaNum] = (float(*)[FitParaNum])tD0;
	float(*pD0)[FitParaNum] = (float(*)[FitParaNum])&D0[0];

	int row = 0;
	int col = 0;

	// = sk.*yk
	divdat = sk[0] * yk[0] + sk[1] * yk[1] + sk[2] * yk[2] + sk[3] * yk[3] + sk[4] * yk[4] + sk[5] * yk[5] + sk[6] * yk[6];
	divdat = 1.0f / divdat;
	// divdat = __fdividef(1.0f , divdat); // 

	float tdat10[FitParaNum];
	float tdat20[FitParaNum];

	tdat10[0] = yk[0] * divdat;
	tdat10[1] = yk[1] * divdat;
	tdat10[2] = yk[2] * divdat;
	tdat10[3] = yk[3] * divdat;
	tdat10[4] = yk[4] * divdat;
	tdat10[5] = yk[5] * divdat;
	tdat10[6] = yk[6] * divdat;

	tdat20[0] = sk[0] * divdat;
	tdat20[1] = sk[1] * divdat;
	tdat20[2] = sk[2] * divdat;
	tdat20[3] = sk[3] * divdat;
	tdat20[4] = sk[4] * divdat;
	tdat20[5] = sk[5] * divdat;
	tdat20[6] = sk[6] * divdat;

	for (row = 0; row<FitParaNum; row++)
	{
		for (col = 0; col<FitParaNum; col++)
		{
			// I-(sk*yk')/(yk'*sk)
			if (row == col) pskyk[row][col] = 1.0f - sk[row] * tdat10[col];
			else            pskyk[row][col] = 0.0f - sk[row] * tdat10[col];

			// (sk*sk')/(yk'*sk)
			psksk[row][col] = sk[row] * tdat20[col];
		}
	}
	// 
	for (row = 0; row<FitParaNum; row++)
	{
		for (col = 0; col<FitParaNum; col++)
		{
			// tD0 = skyk*D0
			ptD0[row][col] = pskyk[row][0] * pD0[0][col] + \
							 pskyk[row][1] * pD0[1][col] + \
							 pskyk[row][2] * pD0[2][col] + \
							 pskyk[row][3] * pD0[3][col] + \
							 pskyk[row][4] * pD0[4][col] + \
							 pskyk[row][5] * pD0[5][col] + \
							 pskyk[row][6] * pD0[6][col];		
		}
	}

	for (row = 0; row<FitParaNum; row++)
	{
		for (col = 0; col<FitParaNum; col++)
		{
			// D0 = D0*yksk
			pD0[row][col] = ptD0[row][0] * pskyk[col][0] + \
							ptD0[row][1] * pskyk[col][1] + \
							ptD0[row][2] * pskyk[col][2] + \
							ptD0[row][3] * pskyk[col][3] + \
							ptD0[row][4] * pskyk[col][4] + \
							ptD0[row][5] * pskyk[col][5] + \
							ptD0[row][6] * pskyk[col][6];	
		}
	}
	// D0=D0+sksk;
	for (row = 0; row<FitParaNum; row++)
	{
		pD0[row][0] = pD0[row][0] + psksk[row][0];
		pD0[row][1] = pD0[row][1] + psksk[row][1];
		pD0[row][2] = pD0[row][2] + psksk[row][2];
		pD0[row][3] = pD0[row][3] + psksk[row][3];
		pD0[row][4] = pD0[row][4] + psksk[row][4];
		pD0[row][5] = pD0[row][5] + psksk[row][5];
		pD0[row][6] = pD0[row][6] + psksk[row][6];
	}
}

__device__ static void MatMulVector_2p(float D0[], float grad[][ThreadsPerBlock], float d0[][ThreadsPerBlock], int tid)
{
	int row;
	float(*pD0)[FitParaNum] = (float(*)[FitParaNum])&D0[0];

	// d0=-D0*grad;    %search direction
	for (row = 0; row< FitParaNum; row++)
	{
		d0[row][tid] = -(pD0[row][0] * grad[0][tid] + 
						 pD0[row][1] * grad[1][tid] + 
						 pD0[row][2] * grad[2][tid] + 
						 pD0[row][3] * grad[3][tid] + 
						 pD0[row][4] * grad[4][tid] + 
						 pD0[row][5] * grad[5][tid] + 
						 pD0[row][6] * grad[6][tid] );
	}

}


#include "bfgs_hd_1.h"



#define OutParaNum			4 // a,x,y,f

#define FitParaNum 			5
#define FittingNum			1

#define IterateNum			8  //8
#define IterateNum_bs		11 //11

// for region size			7

#define RegionLen			9
#define SubRegionWidth		(RegionLen*RegionLen)	// 9x9

#define WholeImageWidth	(RegionLen*(RegionLen+1)) // 56 = 7*8, can gen bank conflict
// avoid bank conflict, if no conflict, +1 can be avoid
#define SharedImageWidth1	(WholeImageWidth+1)


__device__ static void MLELocalization_1p(float subimg[][SharedImageWidth1], float Ininf[][ThreadsPerBlock], float grad[][ThreadsPerBlock], float d0[][ThreadsPerBlock], float D0[], int tid);

__device__ static void PreLocolization_1p(float subimg[][SharedImageWidth1], float Ininf[][ThreadsPerBlock], int tid);
__device__ static float poissonf_1p(float subimg[][SharedImageWidth1], float Ininf[], int tid);
__device__ static void poissonfGradient_1p(float subimg[][SharedImageWidth1], float Ininf[][ThreadsPerBlock], float grad[][ThreadsPerBlock], int tid);

// BFGS core
__device__ static void MatMulVector_1p(float D0[], float grad[][ThreadsPerBlock], float d0[][ThreadsPerBlock], int tid);
__device__ static void ConstructD0_1p(float D0[], float sk[], float yk[], int tid);
__device__ static void IninfConstrain_1p(float Ininf[][ThreadsPerBlock], int tid);
__device__ static void VectorAddMul1_1p(float oVector[], float Ininf[][ThreadsPerBlock], float d0[][ThreadsPerBlock], float coeff, int tid);


// algorithm core codes

__global__ void MLEROILocTop_1p(unsigned short *d_SubRegion, float *d_LocArry,int *d_Point1PosArry, float Offset, float kadc, int FluoNum, int isFPGAProc)
{
	__shared__ float ImageRegion[ThreadsPerBlock][SharedImageWidth1];
	__shared__ float D0[ThreadsPerBlock][FitParaNum*FitParaNum];	// inv of matrix hessian 

	// avoid bank conflict
	__shared__ float Ininf[FitParaNum][ThreadsPerBlock]; // cur position
	__shared__ float grad[FitParaNum][ThreadsPerBlock];  // gradient
	__shared__ float d0[FitParaNum][ThreadsPerBlock];	// direction

	__shared__ float prec[ThreadsPerBlock]; // prec for discard
	__shared__ float prec1[ThreadsPerBlock]; // prec for discard

	int gid =	threadIdx.x+blockDim.x*blockIdx.x;
	int tid =	threadIdx.x;
	if(gid>=FluoNum) return;

	// curid = d_Point1PosArry[gid]
	int BlockOffset = blockDim.x*blockIdx.x;
	int gMemPos;

	// read sub image into shared memory

	int cnt = 0;
	int icnt = 0;

	int SensorId;
	float ypos;
	float prec_th;
	float (*pD0)[FitParaNum]=(float(*)[FitParaNum])&D0[tid][0];

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


	PreLocolization_1p(ImageRegion,Ininf,tid);

	poissonfGradient_1p(ImageRegion,Ininf,grad,tid);

#pragma unroll
	for (cnt = 0; cnt < FitParaNum; cnt++)
	{
		d0[cnt][tid] = -grad[cnt][tid];
	}


	MLELocalization_1p(ImageRegion,Ininf,grad,d0,&D0[tid][0],tid);
	

	Ininf[3][tid]= 0.5f/Ininf[3][tid];

	// used to reject bad points
	prec[tid]  = Ininf[3][tid];
	prec1[tid] = Ininf[3][tid];

	__syncthreads();
	prec[tid] = (prec[tid]> prec[(tid+16)%32]) ? prec[(tid+16)%32] : prec[tid];
	__syncthreads();
	prec[tid] = prec[tid] + prec[(tid+8)%32];
	__syncthreads();
	prec[tid] = prec[tid] + prec[(tid+4)%32];
	__syncthreads();
	prec[tid] = prec[tid] + prec[(tid+2)%32];
	__syncthreads();
	prec[tid] = prec[tid] + prec[(tid+1)%32];
	__syncthreads();
	// if the sigma is too large, then discard it
	prec_th=prec[0]*2.3f/16.0f; // /32*2.2

	if((prec_th>0.0f)&&(prec_th<4.0f))
	{
		prec_th=prec_th;
	}else
	{
		prec_th=4.0f;
	}
	// note the PeakPhoton and bg are electron
	if((Ininf[1][tid]>3.5f)&&(Ininf[1][tid]<5.5f)&&(Ininf[2][tid]>3.5f)&&(Ininf[2][tid]<5.5f)&&(prec1[tid]<prec_th)) //
	{
		if(isFPGAProc==0)
		{
			// online localization

			// true	for camera link to usb
			Ininf[1][tid] = Ininf[1][tid] + ImageRegion[tid][SubRegionWidth + 1] - 3.0f; // X
			Ininf[2][tid] = Ininf[2][tid] + ImageRegion[tid][SubRegionWidth + 2] - 3.0f; // Y
			d_LocArry[gid*OutParaNum + 4] = ImageRegion[tid][SubRegionWidth + 3]; // frame inf


			d_LocArry[gid*OutParaNum + 0] = Ininf[0][tid] * 128.0f;	// A
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

			// for usb to usb only
			Ininf[1][tid] = Ininf[1][tid] + ImageRegion[tid][SubRegionWidth + 0] - 4.0f; // X
			Ininf[2][tid] = Ininf[2][tid] + ImageRegion[tid][SubRegionWidth + 1] - 4.0f; // Y

			d_LocArry[gid*OutParaNum + 0] = Ininf[0][tid] * 128.0f;	// A
			d_LocArry[gid*OutParaNum + 1] = Ininf[1][tid];	// X
			d_LocArry[gid*OutParaNum + 2] = Ininf[2][tid];	// Y
			d_LocArry[gid*OutParaNum + 3] = ImageRegion[tid][SubRegionWidth + 2]; 	// frame inf

//			d_LocArry[gid*OutParaNum + 3] = Ininf[4][tid]; 	// frame inf
//			d_LocArry[gid*OutParaNum + 3] = Ininf[3][tid];	// sigma
		}
	}
	else
	{
		d_LocArry[gid*OutParaNum + 0] = 0.0f; // A
		d_LocArry[gid*OutParaNum + 1] = 0.0f; // X
		d_LocArry[gid*OutParaNum + 2] = 0.0f; // Y
		d_LocArry[gid*OutParaNum + 3] = 0.0f; // sigma


	}
}


__device__ static void MLELocalization_1p(float subimg[][SharedImageWidth1], float Ininf[][ThreadsPerBlock], float grad[][ThreadsPerBlock], float d0[][ThreadsPerBlock], float D0[], int tid)
{
	// adjust d0
	float td0[FitParaNum];
	float td0_total;
	
	float sk[FitParaNum]; // bfgs quasi-newton method
	float yk[FitParaNum];
	float tgrad[FitParaNum];

	int cnt=0;

	int itcnt=0; // iteration number
	
	// find work length, divide 2 method


	float scale;

	float xd[FitParaNum*2];

	float ddat[2];
	float dpos[2];

	int xdsel=0;

	for (itcnt=0;itcnt<IterateNum;itcnt++)
	{
		// adjust d0
		td0[0]=abs(d0[0][tid]);
		td0[1]=abs(d0[1][tid]);
		td0[2]=abs(d0[2][tid]);
		td0[3]=abs(d0[3][tid]);
		td0[4]=abs(d0[4][tid]);

		td0_total=__fdividef((td0[0]+td0[1]+td0[2]+td0[3]+td0[4]),5.0f);

		/*
		if(td0_total<0.1f)
		{
			td0_total=0.1f;
		}
		else if(td0_total<1.0f)
		{
			td0_total=1.0f;
		}
		else
		{
			td0_total=td0_total;
		}
		*/

		td0_total=__fdividef(1.0f,td0_total); // 1.0f/td0_total;

		d0[0][tid]=d0[0][tid]*td0_total;
		d0[1][tid]=d0[1][tid]*td0_total;
		d0[2][tid]=d0[2][tid]*td0_total;
		d0[3][tid]=d0[3][tid]*td0_total;
		d0[4][tid]=d0[4][tid]*td0_total;

		dpos[0]=0.00001f; // scale factor left limit, should not equal to 0 and smaller
		dpos[1]=1.0f; // scale factor right limit, should not lager than 2

		VectorAddMul1_1p(&xd[0],       Ininf, d0, 0.0001f, tid);
		VectorAddMul1_1p(&xd[FitParaNum], Ininf, d0, 1.0f,    tid);

		ddat[0]= poissonf_1p(subimg, &xd[0], tid);
		ddat[1]= poissonf_1p(subimg, &xd[FitParaNum], tid);
		// 
		for(cnt=0;cnt<IterateNum_bs;cnt++)
		{		
			if(ddat[0]<ddat[1])
			{
				xdsel=1; // right shift
			}
			else
			{
				xdsel=0; // left shift
			}

			dpos[xdsel]=(dpos[0]+dpos[1])*0.5f; //  /2.0f which one shift

			if(cnt<IterateNum_bs-1)
			{
				VectorAddMul1_1p(&xd[xdsel*FitParaNum], Ininf, d0, dpos[xdsel], tid);// xd=ininf+d0*scale
				ddat[xdsel]= poissonf_1p(subimg, &xd[xdsel*FitParaNum], tid);
			}

		}
		scale=(dpos[0]+dpos[1])*0.5f; //  /2.0f calculated direction

		// sk 		= d0*base;
		// inInf 	= inInf+d0*base;
		sk[0] = d0[0][tid]*scale;
		sk[1] = d0[1][tid]*scale;
		sk[2] = d0[2][tid]*scale;
		sk[3] = d0[3][tid]*scale;
		sk[4] = d0[4][tid]*scale;

		Ininf[0][tid] = Ininf[0][tid] + sk[0];
		Ininf[1][tid] = Ininf[1][tid] + sk[1];
		Ininf[2][tid] = Ininf[2][tid] + sk[2];
		Ininf[3][tid] = Ininf[3][tid] + sk[3];
		Ininf[4][tid] = Ininf[4][tid] + sk[4];

		IninfConstrain_1p(Ininf,tid); // 


		if(itcnt<IterateNum-1)
		{
			// tgrad = grad;
			tgrad[0] = grad[0][tid];
			tgrad[1] = grad[1][tid];
			tgrad[2] = grad[2][tid];
			tgrad[3] = grad[3][tid];
			tgrad[4] = grad[4][tid];

			poissonfGradient_1p(subimg,Ininf,grad,tid);

			yk[0] = grad[0][tid] - tgrad[0];
			yk[1] = grad[1][tid] - tgrad[1];
			yk[2] = grad[2][tid] - tgrad[2];
			yk[3] = grad[3][tid] - tgrad[3];
			yk[4] = grad[4][tid] - tgrad[4];

			ConstructD0_1p(D0, sk, yk, tid);
			MatMulVector_1p(D0,grad,d0,tid);
		}
	}
}


__device__ static float poissonf_1p(float subimg[][SharedImageWidth1],float Ininf[],int tid)
{
	int row,col;

	float rowpos,colpos;
	float IVal;
	float VPoss=0;
	float pixval;
	float tdat;

	float (*pSubImg)[RegionLen]=(float(*)[RegionLen])&subimg[tid][0];

	if(Ininf[0]< 0.0f)Ininf[0]=0.0f;
	if(Ininf[1]< 0.0f)Ininf[1]=0.0f;
	if(Ininf[2]< 0.0f)Ininf[2]=0.0f;
	if(Ininf[3]< 0.0f)Ininf[3]=0.0f;
	if(Ininf[4]< 0.0f)Ininf[4]=0.0f;

	for(row=1; row < 8; row++) // RegionLen-1
	{
		for(col=1; col < 8; col++)
		{
			rowpos = row + 0.5f; // -0.5
			colpos = col + 0.5f; // -0.5
			

			tdat = -((colpos-Ininf[1])*(colpos-Ininf[1]) + (rowpos-Ininf[2])*(rowpos-Ininf[2]))*Ininf[3];

			// get model intensity
			IVal = 128.0f*Ininf[0] * __expf(tdat) + Ininf[4] * 32.0f;

			// pf goal function
			pixval = pSubImg[row][col]; //subimg[tid][row*RegionLen + col]; // coloffset+

			//			tdat = IVal - pixval - pixval*(__logf(IVal) - __logf(pixval + 80.12345f)); //to avoid negative value in logf function, that will be a disaster
			tdat = IVal - pixval - pixval*(__logf(IVal) - __logf(fmaxf(pixval, 1.0f))); //to avoid negative value in logf function, that will be a disaster

			VPoss = VPoss + IVal;
		}
	}
	return VPoss;
}

__device__ static void PreLocolization_1p(float subimg[][SharedImageWidth1],float Ininf[][ThreadsPerBlock],int tid)
{
//	int rowoff=tid*8;

	float sum_bg=0;


	float (*pSubImg)[RegionLen]=(float(*)[RegionLen])&subimg[tid][0];

	sum_bg = 	pSubImg[1][1] + pSubImg[1][2] + pSubImg[1][3] + pSubImg[1][4] + pSubImg[1][5] + pSubImg[1][6] + pSubImg[1][7] +\
				pSubImg[7][1] + pSubImg[7][2] + pSubImg[7][3] + pSubImg[7][4] + pSubImg[7][5] + pSubImg[7][6] + pSubImg[7][7] +\
				pSubImg[2][1] + pSubImg[3][1] + pSubImg[4][1] + pSubImg[5][1] + pSubImg[6][1] +\
				pSubImg[2][7] + pSubImg[3][7] + pSubImg[4][7] + pSubImg[5][7] + pSubImg[6][7];

	sum_bg = sum_bg / 24.0f;// *0.041667f; // /24.0f; // background
	

	Ininf[0][tid] = (pSubImg[4][4]-sum_bg)/128.0f; // /256.0f;	// Amplitude

	Ininf[1][tid] = 4.5f; // sumDatx / sumdat - 0.5f;	// x
	Ininf[2][tid] = 4.5f; // sumDaty / sumdat - 0.5f;	// y
	Ininf[3][tid] = ParaSigma;	// sigma

	Ininf[4][tid] = sum_bg/32.0f; // /128.0f;	// bg
	
}


__device__ static void poissonfGradient_1p(float subimg[][SharedImageWidth1],float Ininf[][ThreadsPerBlock],float grad[][ThreadsPerBlock],int tid)
{
	float tIninf[FitParaNum];
	float tgradn;
	float tgradp;
	//lv 0
	tIninf[0]=Ininf[0][tid];
	tIninf[1]=Ininf[1][tid];
	tIninf[2]=Ininf[2][tid];
	tIninf[3]=Ininf[3][tid];
	tIninf[4]=Ininf[4][tid];

	int cnt;
	for(cnt=0;cnt<FitParaNum;cnt++)
	{
		tIninf[cnt]=tIninf[cnt]-0.002f;
		tgradn=poissonf_1p(subimg,tIninf,tid);

		tIninf[cnt]=tIninf[cnt]+0.004f;
		tgradp=poissonf_1p(subimg,tIninf,tid);

		grad[cnt][tid]=(tgradp-tgradn)*250.0f; // /0.004f;
		tIninf[cnt]=tIninf[cnt]-0.002f;
	}
}


__device__ static void VectorAddMul1_1p(float oVector[], float Ininf[][ThreadsPerBlock], float d0[][ThreadsPerBlock], float coeff, int tid)
{
	oVector[0] = d0[0][tid] * coeff + Ininf[0][tid];
	oVector[1] = d0[1][tid] * coeff + Ininf[1][tid];
	oVector[2] = d0[2][tid] * coeff + Ininf[2][tid];
	oVector[3] = d0[3][tid] * coeff + Ininf[3][tid];
	oVector[4] = d0[4][tid] * coeff + Ininf[4][tid];
}

__device__ static void IninfConstrain_1p(float Ininf[][ThreadsPerBlock], int tid)
{
	if (Ininf[0][tid]<0.01f) Ininf[0][tid] = 0.01f; // A
	if (Ininf[1][tid]<0.01f) Ininf[1][tid] = 0.01f; // X
	if (Ininf[2][tid]<0.01f) Ininf[2][tid] = 0.01f; // Y
	if (Ininf[3][tid]<0.01f) Ininf[3][tid] = 0.01f; // sigma
	if (Ininf[4][tid]<0.00f) Ininf[4][tid] = 0.00f; // bg

}

__device__ static void ConstructD0_1p(float D0[], float sk[FitParaNum], float yk[FitParaNum], int tid)
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
	divdat = sk[0] * yk[0] + sk[1] * yk[1] + sk[2] * yk[2] + sk[3] * yk[3] + sk[4] * yk[4];
	divdat = 1.0f / divdat;
	// divdat = __fdividef(1.0f , divdat); // 

	float tdat10[FitParaNum];
	float tdat20[FitParaNum];

	tdat10[0] = yk[0] * divdat;
	tdat10[1] = yk[1] * divdat;
	tdat10[2] = yk[2] * divdat;
	tdat10[3] = yk[3] * divdat;
	tdat10[4] = yk[4] * divdat;

	tdat20[0] = sk[0] * divdat;
	tdat20[1] = sk[1] * divdat;
	tdat20[2] = sk[2] * divdat;
	tdat20[3] = sk[3] * divdat;
	tdat20[4] = sk[4] * divdat;

	for (row = 0; row<FitParaNum; row++)
	{
		for (col = 0; col<FitParaNum; col++)
		{
			// I-(sk*yk')/(yk'*sk)
			if (row == col) pskyk[row][col] = 1.0f - sk[row] * tdat10[col];
			else         pskyk[row][col] = 0.0f - sk[row] * tdat10[col];

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
				pskyk[row][4] * pD0[4][col];
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
				ptD0[row][4] * pskyk[col][4];
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
	}
}

__device__ static void MatMulVector_1p(float D0[], float grad[][ThreadsPerBlock], float d0[][ThreadsPerBlock], int tid)
{
	int row;
	float(*pD0)[FitParaNum] = (float(*)[FitParaNum])&D0[0];

	// d0=-D0*grad;    %search direction
	for (row = 0; row<FitParaNum; row++)
	{
		d0[row][tid] = -(pD0[row][0] * grad[0][tid] + pD0[row][1] * grad[1][tid] + pD0[row][2] * grad[2][tid] + pD0[row][3] * grad[3][tid] + pD0[row][4] * grad[4][tid]);
	}

}


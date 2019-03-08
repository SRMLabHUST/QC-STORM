#include "LDROIExtraction.h"



// calculate background intensity
#define BkgBatchLineNum					8

#define BkgFilterSize					9			
#define BkgFilterSize_Half				((int)(BkgFilterSize/2))

// seperable convolution
#define LineFilterSize				5



// calculate background intensity
__global__ void gpuBackgroundIntensityCalc(unsigned short *d_RawImg, unsigned short *d_BackgroundImage, int ImageWidth, int ImageHigh);


// seperable convolution horizontal
__global__ void gpuImageConv_Seperable_H(unsigned short *d_iImage, unsigned short *d_oImage, float *d_LineFilter, int ImageWidth, int ImageHigh);

// seperable convolution vertical
__global__ void gpuImageConv_Seperable_V(unsigned short *d_iImage, unsigned short *d_oImage, float *d_LineFilter, int ImageWidth, int ImageHigh);

// ImageSubstract, d_ImageC = d_ImageA - d_ImageB
__global__ void gpuImageSubstract(unsigned short *d_ImageA, unsigned short *d_ImageB, unsigned short *d_ImageC, int ImageWidth, int ImageHigh);


__global__ void gpuImageVarianceCalcOnSelectRegion(unsigned short *d_RawImg, float *d_MeanDataX, float *d_MeanDataX2, int ImageCenterX, int ImageCenterY, int StdImageSize, int ImageWidth, int ImageHigh);



// calculate background intensity
void BackgroundIntensityCalc(unsigned short *d_RawImg, unsigned short *d_BackgroundImage, int ImageWidth, int ImageHigh, cudaStream_t cstream)
{
	int ImageWidth_Thread = GetImageWidth_Thread(ImageWidth);

	int TotalThreadNum = ((ImageHigh + BkgBatchLineNum - 1) / BkgBatchLineNum) * ImageWidth_Thread;

	int BlockDim = ThreadsPerBlock;
	int BlockNum = (TotalThreadNum + ThreadsPerBlock - 1) / ThreadsPerBlock;
	
	gpuBackgroundIntensityCalc << < BlockNum, BlockDim, 0, cstream >> > (d_RawImg, d_BackgroundImage, ImageWidth, ImageHigh);

}


// calculate background intensity
__global__ void gpuBackgroundIntensityCalc(unsigned short *d_RawImg, unsigned short *d_BackgroundImage, int ImageWidth, int ImageHigh)
{
	enum {
		SharedMemHigh = BkgBatchLineNum + 2 * BkgFilterSize_Half,
		SharedMemWidth = ThreadsPerBlock + 2 * BkgFilterSize_Half,

		BkgCmpArrayLen = BkgFilterSize - 1,
	};

	// side pixels, valid pixels, side pixels
	__shared__ unsigned short ImageRegion[SharedMemHigh][SharedMemWidth];


	int gid = threadIdx.x + blockDim.x*blockIdx.x;
	int tid = threadIdx.x;


	int ImageWidth_Thread = GetImageWidth_Thread(ImageWidth);


	// x,y pos of current thread
	int CurXPos = gid % ImageWidth_Thread;
	int CurYPos = (gid / ImageWidth_Thread)*BkgBatchLineNum;


	int XPos_GMem = 0;
	int YPos_GMem = 0;

	// load global memory to shared memory, first part
#pragma unroll
	for (int cnt = 0; cnt < SharedMemHigh; cnt++)
	{
		XPos_GMem = CurXPos - BkgFilterSize_Half;
		YPos_GMem = CurYPos - BkgFilterSize_Half + cnt;

		XPos_GMem = max(XPos_GMem, 0);
		YPos_GMem = max(YPos_GMem, 0);

		XPos_GMem = min(XPos_GMem, ImageWidth - 1);
		YPos_GMem = min(YPos_GMem, ImageHigh - 1);

		int Addr_GMem = YPos_GMem*ImageWidth + XPos_GMem;

		ImageRegion[cnt][tid] = d_RawImg[Addr_GMem];
	}

	// load global memory to shared memory, second part
#pragma unroll
	for (int cnt = 0; cnt < SharedMemHigh; cnt++)
	{
		if (tid < SharedMemWidth - ThreadsPerBlock)
		{
			XPos_GMem = CurXPos - BkgFilterSize_Half;
			YPos_GMem = CurYPos - BkgFilterSize_Half + cnt;

			// right part
			XPos_GMem = XPos_GMem + ThreadsPerBlock;

			XPos_GMem = max(XPos_GMem, 0);
			YPos_GMem = max(YPos_GMem, 0);

			XPos_GMem = min(XPos_GMem, ImageWidth - 1);
			YPos_GMem = min(YPos_GMem, ImageHigh - 1);

			int Addr_GMem = YPos_GMem*ImageWidth + XPos_GMem;

			ImageRegion[cnt][tid + ThreadsPerBlock] = d_RawImg[Addr_GMem];
		}
	}

	__syncthreads();

	// filter part

	// margin pixel values, up, right, bottom, left
	float PixelNumL1[BkgCmpArrayLen];
	float PixelNumL2[BkgCmpArrayLen];
	float PixelNumL3[BkgCmpArrayLen];
	float PixelNumL4[BkgCmpArrayLen];


	for (int cnt = 0; cnt < BkgBatchLineNum; cnt++)
	{
		int XCenterPos_SMem = tid + BkgFilterSize_Half;
		int YCenterPos_SMem = cnt + BkgFilterSize_Half;

#pragma unroll
		for (int i = 0; i < BkgCmpArrayLen; i++)
		{
			PixelNumL1[i] = ImageRegion[YCenterPos_SMem - BkgFilterSize_Half][XCenterPos_SMem - BkgFilterSize_Half + i];
			PixelNumL3[i] = ImageRegion[YCenterPos_SMem + BkgFilterSize_Half][XCenterPos_SMem - BkgFilterSize_Half + i + 1];
			PixelNumL2[i] = ImageRegion[YCenterPos_SMem - BkgFilterSize_Half + i][XCenterPos_SMem + BkgFilterSize_Half];
			PixelNumL4[i] = ImageRegion[YCenterPos_SMem - BkgFilterSize_Half + i + 1][XCenterPos_SMem - BkgFilterSize_Half];
		}

#pragma unroll
		for (int i = 0; i < BkgCmpArrayLen; i++)
		{
			PixelNumL1[i] = min(PixelNumL1[i], PixelNumL2[i]);
			PixelNumL3[i] = min(PixelNumL3[i], PixelNumL4[i]);
		}

		float MeanData = 0;

#pragma unroll
		for (int i = 0; i < BkgCmpArrayLen; i++)
		{
			MeanData += min(PixelNumL1[i], PixelNumL3[BkgCmpArrayLen - 1 - i]);
		}

		MeanData /= BkgCmpArrayLen;

		XPos_GMem = CurXPos;
		YPos_GMem = CurYPos + cnt;

		if ((XPos_GMem < ImageWidth) && (YPos_GMem < ImageHigh))
		{
			d_BackgroundImage[YPos_GMem*ImageWidth + XPos_GMem] = MeanData;
		}
	}
}



// seperable convolution horizontal
void ImageConv_Seperable_H(unsigned short *d_iImage, unsigned short *d_oImage, float *d_LineFilter, int ImageWidth, int ImageHigh, cudaStream_t cstream)
{

	int ImageWidth_Thread = GetImageWidth_Thread(ImageWidth);

	int TotalThreadNum = ImageHigh * ImageWidth_Thread;

	int BlockDim = ThreadsPerBlock;
	int BlockNum = (TotalThreadNum + ThreadsPerBlock - 1) / ThreadsPerBlock;

	gpuImageConv_Seperable_H << < BlockNum, BlockDim, 0, cstream >> > (d_iImage, d_oImage, d_LineFilter, ImageWidth, ImageHigh);


}


// seperable convolution horizontal
__global__ void gpuImageConv_Seperable_H(unsigned short *d_iImage, unsigned short *d_oImage, float *d_LineFilter, int ImageWidth, int ImageHigh)
{
	enum {
		FilterSize_Half = (int)(LineFilterSize / 2),

		SharedMemWidth = ThreadsPerBlock + 2 * FilterSize_Half,
	};

	// side pixels, valid pixels, side pixels
	__shared__ unsigned short ImageRegion[SharedMemWidth];


	int gid = threadIdx.x + blockDim.x*blockIdx.x;
	int tid = threadIdx.x;

	// gaurantee pixels in a thread block are in the same line
	int ImageWidth_Thread = GetImageWidth_Thread(ImageWidth);


	int CurXPos = gid % ImageWidth_Thread;
	int CurYPos = gid / ImageWidth_Thread;


	// read global memory to shared memory, first part

	int XPos_GMem = CurXPos - FilterSize_Half;
	int YPos_GMem = CurYPos;

	XPos_GMem = max(XPos_GMem, 0);
	YPos_GMem = max(YPos_GMem, 0);

	XPos_GMem = min(XPos_GMem, ImageWidth - 1);
	YPos_GMem = min(YPos_GMem, ImageHigh - 1);


	ImageRegion[tid] = d_iImage[YPos_GMem*ImageWidth + XPos_GMem];


	// read global memory to shared memory, second part
	if (tid < SharedMemWidth - ThreadsPerBlock)
	{
		XPos_GMem = CurXPos - FilterSize_Half;
		YPos_GMem = CurYPos;

		XPos_GMem += ThreadsPerBlock;

		XPos_GMem = max(XPos_GMem, 0);
		YPos_GMem = max(YPos_GMem, 0);

		XPos_GMem = min(XPos_GMem, ImageWidth - 1);
		YPos_GMem = min(YPos_GMem, ImageHigh - 1);


		ImageRegion[tid + ThreadsPerBlock] = d_iImage[YPos_GMem*ImageWidth + XPos_GMem];
	}

	__syncthreads();

	// filter part

	float FilteredData = 0;

#pragma unroll
	for (int cnt = 0; cnt < LineFilterSize; cnt++)
	{
		FilteredData += ImageRegion[tid + cnt] * d_LineFilter[cnt];
	}

	if ((CurXPos < ImageWidth) && (CurYPos < ImageHigh))
	{
		d_oImage[CurYPos*ImageWidth + CurXPos] = FilteredData;
	}

}


// seperable convolution vertical
void ImageConv_Seperable_V(unsigned short *d_iImage, unsigned short *d_oImage, float *d_LineFilter, int ImageWidth, int ImageHigh, cudaStream_t cstream)
{

	int ImageHigh_Thread = GetImageWidth_Thread(ImageHigh);

	int TotalThreadNum = ImageHigh_Thread * ImageWidth;

	int BlockDim = ThreadsPerBlock;
	int BlockNum = (TotalThreadNum + ThreadsPerBlock - 1) / ThreadsPerBlock;


	gpuImageConv_Seperable_V << < BlockNum, BlockDim, 0, cstream >> > (d_iImage, d_oImage, d_LineFilter, ImageWidth, ImageHigh);

}


// seperable convolution vertical
__global__ void gpuImageConv_Seperable_V(unsigned short *d_iImage, unsigned short *d_oImage, float *d_LineFilter, int ImageWidth, int ImageHigh)
{
	enum {
		FilterSize_Half = (int)(LineFilterSize / 2),

		SharedMemHigh = ThreadsPerBlock + 2 * FilterSize_Half,
	};

	// side pixels, valid pixels, side pixels
	__shared__ unsigned short ImageRegion[SharedMemHigh];


	int gid = threadIdx.x + blockDim.x*blockIdx.x;
	int tid = threadIdx.x;


	int ImageHigh_Thread = GetImageWidth_Thread(ImageHigh);


	int CurYPos = gid % ImageHigh_Thread;
	int CurXPos = gid / ImageHigh_Thread;

	// read global memory to shared memory, first part
	int XPos_GMem = CurXPos;
	int YPos_GMem = CurYPos - FilterSize_Half;

	XPos_GMem = max(XPos_GMem, 0);
	YPos_GMem = max(YPos_GMem, 0);

	XPos_GMem = min(XPos_GMem, ImageWidth - 1);
	YPos_GMem = min(YPos_GMem, ImageHigh - 1);

	int Addr_GMem = YPos_GMem*ImageWidth + XPos_GMem;

	ImageRegion[tid] = d_iImage[Addr_GMem];

	// read global memory to shared memory, second part
	if (tid < SharedMemHigh - ThreadsPerBlock)
	{
		XPos_GMem = CurXPos;
		YPos_GMem = CurYPos - FilterSize_Half;

		YPos_GMem += ThreadsPerBlock;

		XPos_GMem = max(XPos_GMem, 0);
		YPos_GMem = max(YPos_GMem, 0);

		XPos_GMem = min(XPos_GMem, ImageWidth - 1);
		YPos_GMem = min(YPos_GMem, ImageHigh - 1);

		Addr_GMem = YPos_GMem*ImageWidth + XPos_GMem;

		ImageRegion[tid + ThreadsPerBlock] = d_iImage[Addr_GMem];
	}

	__syncthreads();

	// filter part

	float FilteredData = 0;

#pragma unroll
	for (int cnt = 0; cnt < LineFilterSize; cnt++)
	{
		FilteredData += ImageRegion[tid + cnt] * d_LineFilter[cnt];
	}

	if ((CurXPos < ImageWidth) && (CurYPos < ImageHigh))
	{
		d_oImage[CurYPos*ImageWidth + CurXPos] = FilteredData;
	}

}

// ImageSubstract, d_ImageC = d_ImageA - d_ImageB
void ImageSubstract(unsigned short *d_ImageA, unsigned short *d_ImageB, unsigned short *d_ImageC, int ImageWidth, int ImageHigh, cudaStream_t cstream)
{

	int TotalThreadNum = ImageWidth * ImageHigh;

	int BlockDim = ThreadsPerBlock;
	int BlockNum = (TotalThreadNum + ThreadsPerBlock - 1) / ThreadsPerBlock;

	gpuImageSubstract << < BlockNum, BlockDim, 0, cstream >> > (d_ImageA, d_ImageB, d_ImageC, ImageWidth, ImageHigh);

}


// ImageSubstract, d_ImageC = d_ImageA - d_ImageB
__global__ void gpuImageSubstract(unsigned short *d_ImageA, unsigned short *d_ImageB, unsigned short *d_ImageC, int ImageWidth, int ImageHigh)
{
	int gid = threadIdx.x + blockDim.x*blockIdx.x;

	if (gid < ImageWidth * ImageHigh)
	{
		unsigned short DataA = d_ImageA[gid];
		unsigned short DataB = d_ImageB[gid];

		d_ImageC[gid] = DataA >= DataB ? (DataA - DataB) : 0;
	}
}



void ImageVarianceCalcOnSelectRegion(unsigned short *d_RawImg, float *d_MeanDataX, float *d_MeanDataX2, int ImageCenterX, int ImageCenterY, int StdImageSize, int ImageWidth, int ImageHigh, cudaStream_t cstream)
{
	int TotalThreadNum = StdImageSize * StdImageSize;

	int BlockDim = ThreadsPerBlock;
	int BlockNum = (TotalThreadNum + ThreadsPerBlock - 1) / ThreadsPerBlock;

	gpuImageVarianceCalcOnSelectRegion << < BlockNum, BlockDim, 0, cstream >> > (d_RawImg, d_MeanDataX, d_MeanDataX2, ImageCenterX, ImageCenterY, StdImageSize, ImageWidth, ImageHigh);

}

__global__ void gpuImageVarianceCalcOnSelectRegion(unsigned short *d_RawImg, float *d_MeanDataX, float *d_MeanDataX2, int ImageCenterX, int ImageCenterY, int StdImageSize, int ImageWidth, int ImageHigh)
{
	int gid = threadIdx.x + blockDim.x*blockIdx.x;

	float TotalDataNum = StdImageSize * StdImageSize;

	if (gid < StdImageSize* StdImageSize)
	{
		int CurXPos = gid % StdImageSize;
		int CurYPos = gid / StdImageSize;

		// read data from global memory
		int XPos_GMem = ImageCenterX + CurXPos;
		int YPos_GMem = ImageCenterY + CurYPos;

		unsigned short data1 = d_RawImg[YPos_GMem*ImageWidth + XPos_GMem];
		unsigned short data2 = d_RawImg[YPos_GMem*ImageWidth + (XPos_GMem - StdImageSize)];
		unsigned short data3 = d_RawImg[(YPos_GMem - StdImageSize)*ImageWidth + XPos_GMem];
		unsigned short data4 = d_RawImg[(YPos_GMem - StdImageSize)*ImageWidth + (XPos_GMem - StdImageSize)];

		data1 = min(data1, data2);
		data3 = min(data3, data4);

		data1 = min(data1, data3);

		atomicAdd(d_MeanDataX, (float)data1 / TotalDataNum);
		atomicAdd(d_MeanDataX2, (float)(data1*data1) / TotalDataNum);

	}

}

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

// line seperable Gaussian filter kernel
__constant__ float d_LineFilterH_Signal[LineFilterSize];


void LDROIExtractData_TypeDef::FilterInit(int ROISize, int LocType, int MultiEmitterFitEn)
{

	// filter sigma
	float Sigma7 = 0.5f;
	float Sigma11 = 1.0f;

	float SigmaS = Sigma7 + (ROISize - 7)*(Sigma11 - Sigma7) / (11 - 7);

	if (SigmaS < 0.5f)SigmaS = 0.5f;

	if (LocType == LocType_AS3D)
	{
		SigmaS *= 1.5f;
	}
	else if ((MultiEmitterFitEn == 0) && (LocType == LocType_GS2D))
	{
		SigmaS *= 1.2f;
	}

	//
	SigmaS = -1 / (2 * SigmaS*SigmaS);

	float LineFilter_Center = (int)(LineFilterSize / 2);

	for (int cnt = 0; cnt < LineFilterSize; cnt++)
	{
		h_LineFilterH_Signal[cnt] = expf((cnt - LineFilter_Center)*(cnt - LineFilter_Center)*SigmaS);
	}

	// normalize filter
	float SumData = 0;
	for (int i = 0; i < LineFilterSize; i++)
	{
		for (int j = 0; j < LineFilterSize; j++)
		{
			SumData += h_LineFilterH_Signal[i] * h_LineFilterH_Signal[j];
		}
	}

	SumData = sqrtf(1.0f / SumData);

	for (int cnt = 0; cnt < LineFilterSize; cnt++)
	{
		h_LineFilterH_Signal[cnt] = h_LineFilterH_Signal[cnt] * SumData;
	}


	cudaStream_t loc_stream1;
	cudaStreamCreate(&loc_stream1);


	// copy host filter data to constant memory

	cudaMemcpyToSymbolAsync(d_LineFilterH_Signal, h_LineFilterH_Signal, LineFilterSize * sizeof(float), 0, cudaMemcpyHostToDevice, loc_stream1);


	cudaStreamSynchronize(loc_stream1);
	cudaStreamDestroy(loc_stream1);
}



void LDROIExtractData_TypeDef::ImageFiltering(int LocType, int MultiEmitterFitEn, int ImageWidth, int ImageHigh, int BatchedImageNum, cudaStream_t cstream)
{
	const int BatchedImageHigh = BatchedImageNum*ImageHigh;


	BackgroundIntensityRemove(d_RawImg, d_RawImg_Smoothed, ImageWidth, BatchedImageHigh, cstream);


	ImageVarianceCalc(d_RawImg_Smoothed, ImageWidth, BatchedImageHigh, cstream);


	if ((LocType == LocType_GS2D) && (MultiEmitterFitEn > 0))
	{
		// smooth d_RawImg
		ImageConv_Seperable_H(d_RawImg_Smoothed, d_LineFilterImage_t, ImageWidth, BatchedImageHigh, cstream);
		ImageConv_Seperable_V(d_LineFilterImage_t, d_LineFilterImage_t1, ImageWidth, BatchedImageHigh, cstream);

		// sharpen filtered image
		ImageSharpen(d_LineFilterImage_t1, d_RawImg_Smoothed, ImageWidth, BatchedImageHigh, cstream);

	}
	else
	{
		// smooth d_RawImg
		ImageConv_Seperable_H(d_RawImg_Smoothed, d_LineFilterImage_t, ImageWidth, BatchedImageHigh, cstream);
		ImageConv_Seperable_V(d_LineFilterImage_t, d_RawImg_Smoothed, ImageWidth, BatchedImageHigh, cstream);

	}

	// debug
	//	cudaMemcpyAsync(h_RawImg, d_RawImg_Smoothed, ImageWidth*BatchedImageHigh * sizeof(short), cudaMemcpyDeviceToHost, cstream);
	//	printf("var:%f\n", ImageVariance);


	cudaStreamSynchronize(cstream);

}


// calculate background intensity
__global__ void gpuRemoveBackgroundIntensity(unsigned short *d_RawImg, unsigned short *d_BackgroundImage, int ImageWidth, int ImageHigh)
{
	enum {
		SharedMemHigh = BlockFilter_BatchLineNum + 2 * BkgFilterSize_Half,
		SharedMemWidth = ThreadsPerBlock + 2 * BkgFilterSize_Half,

		BkgCmpArrayLen = BkgFilterSize - 1,
	};

	// side pixels, valid pixels, side pixels
	__shared__ unsigned short ImageROI[SharedMemHigh][SharedMemWidth];

	__shared__ float BkgVal[ThreadsPerBlock];


	int gid = threadIdx.x + blockDim.x*blockIdx.x;
	int tid = threadIdx.x;


	int ImageWidth_Thread = GetImageWidth_Thread(ImageWidth);


	// x,y pos of current thread
	int CurXPos = gid % ImageWidth_Thread;
	int CurYPos = (gid / ImageWidth_Thread)*BlockFilter_BatchLineNum;


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

		ImageROI[cnt][tid] = d_RawImg[Addr_GMem];
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

			ImageROI[cnt][tid + ThreadsPerBlock] = d_RawImg[Addr_GMem];
		}
	}

	__syncthreads();

	// filter part

	// margin pixel values, up, right, bottom, left
	float PixelNumL1[BkgCmpArrayLen];
	float PixelNumL2[BkgCmpArrayLen];
	float PixelNumL3[BkgCmpArrayLen];
	float PixelNumL4[BkgCmpArrayLen];


	for (int cnt = 0; cnt < BlockFilter_BatchLineNum; cnt++)
	{
		int XCenterPos_SMem = tid + BkgFilterSize_Half;
		int YCenterPos_SMem = cnt + BkgFilterSize_Half;

#pragma unroll
		for (int i = 0; i < BkgCmpArrayLen; i++)
		{
			PixelNumL1[i] = ImageROI[YCenterPos_SMem - BkgFilterSize_Half][XCenterPos_SMem - BkgFilterSize_Half + i];
			PixelNumL3[i] = ImageROI[YCenterPos_SMem + BkgFilterSize_Half][XCenterPos_SMem - BkgFilterSize_Half + i + 1];
			PixelNumL2[i] = ImageROI[YCenterPos_SMem - BkgFilterSize_Half + i][XCenterPos_SMem + BkgFilterSize_Half];
			PixelNumL4[i] = ImageROI[YCenterPos_SMem - BkgFilterSize_Half + i + 1][XCenterPos_SMem - BkgFilterSize_Half];
		}

#pragma unroll
		for (int i = 0; i < BkgCmpArrayLen; i++)
		{
			PixelNumL1[i] = min(PixelNumL1[i], PixelNumL2[i]);
			PixelNumL3[i] = min(PixelNumL3[i], PixelNumL4[i]);
		}

		float CurBkgVal = 0;

#pragma unroll
		for (int i = 0; i < BkgCmpArrayLen; i++)
		{
			CurBkgVal += min(PixelNumL1[i], PixelNumL3[BkgCmpArrayLen - 1 - i]);
		}

		CurBkgVal /= BkgCmpArrayLen;

		// average in threads
		BkgVal[tid] = CurBkgVal;

		__syncthreads();

		int Sel0 = tid - 1;
		int Sel1 = tid ;
		int Sel2 = tid + 1;
		Sel0 = max(Sel0, 0);
		Sel2 = min(Sel2, ThreadsPerBlock - 1);

		float MeanBkg = (BkgVal[Sel0] + BkgVal[Sel1] + BkgVal[Sel2]) / 3.0f;

		XPos_GMem = CurXPos;
		YPos_GMem = CurYPos + cnt;

		/**/
		float CenterPixelData = ImageROI[YCenterPos_SMem][XCenterPos_SMem];
		CenterPixelData = CenterPixelData - MeanBkg;

		CenterPixelData = max(CenterPixelData, 0.0f);
		

		if ((XPos_GMem < ImageWidth) && (YPos_GMem < ImageHigh))
		{
			d_BackgroundImage[YPos_GMem*ImageWidth + XPos_GMem] = CenterPixelData;
		}
	}
}


// calculate background intensity
void BackgroundIntensityRemove(unsigned short *d_RawImg, unsigned short *d_BackgroundImage, int ImageWidth, int ImageHigh, cudaStream_t cstream)
{
	int ImageWidth_Thread = GetImageWidth_Thread(ImageWidth);

	int TotalThreadNum = ((ImageHigh + BlockFilter_BatchLineNum - 1) / BlockFilter_BatchLineNum) * ImageWidth_Thread;

	int BlockDim = ThreadsPerBlock;
	int BlockNum = (TotalThreadNum + ThreadsPerBlock - 1) / ThreadsPerBlock;

	gpuRemoveBackgroundIntensity << < BlockNum, BlockDim, 0, cstream >> > (d_RawImg, d_BackgroundImage, ImageWidth, ImageHigh);

}


// seperable convolution horizontal
__global__ void gpuImageConv_Seperable_H(unsigned short *d_iImage, unsigned short *d_oImage, int ImageWidth, int ImageHigh)
{
	enum {
		FilterSize_Half = (int)(LineFilterSize / 2),

		SharedMemWidth = ThreadsPerBlock + 2 * FilterSize_Half,
	};

	// side pixels, valid pixels, side pixels
	__shared__ unsigned short ImageROI[SharedMemWidth];


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


	ImageROI[tid] = d_iImage[YPos_GMem*ImageWidth + XPos_GMem];


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


		ImageROI[tid + ThreadsPerBlock] = d_iImage[YPos_GMem*ImageWidth + XPos_GMem];
	}

	__syncthreads();

	// filter part

	float FilteredData = 0;

#pragma unroll
	for (int cnt = 0; cnt < LineFilterSize; cnt++)
	{
		FilteredData += ImageROI[tid + cnt] * d_LineFilterH_Signal[cnt];
	}

	if ((CurXPos < ImageWidth) && (CurYPos < ImageHigh))
	{
		d_oImage[CurYPos*ImageWidth + CurXPos] = FilteredData;
	}

}



// seperable convolution horizontal
void ImageConv_Seperable_H(unsigned short *d_iImage, unsigned short *d_oImage, int ImageWidth, int ImageHigh, cudaStream_t cstream)
{

	int ImageWidth_Thread = GetImageWidth_Thread(ImageWidth);

	int TotalThreadNum = ImageHigh * ImageWidth_Thread;

	int BlockDim = ThreadsPerBlock;
	int BlockNum = (TotalThreadNum + ThreadsPerBlock - 1) / ThreadsPerBlock;

	gpuImageConv_Seperable_H << < BlockNum, BlockDim, 0, cstream >> > (d_iImage, d_oImage, ImageWidth, ImageHigh);


}



// seperable convolution vertical
__global__ void gpuImageConv_Seperable_V(unsigned short *d_iImage, unsigned short *d_oImage, int ImageWidth, int ImageHigh)
{
	enum {
		FilterSize_Half = (int)(LineFilterSize / 2),

		SharedMemHigh = ThreadsPerBlock + 2 * FilterSize_Half,
	};

	// side pixels, valid pixels, side pixels
	__shared__ unsigned short ImageROI[SharedMemHigh];


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

	ImageROI[tid] = d_iImage[Addr_GMem];

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

		ImageROI[tid + ThreadsPerBlock] = d_iImage[Addr_GMem];
	}

	__syncthreads();

	// filter part

	float FilteredData = 0;

#pragma unroll
	for (int cnt = 0; cnt < LineFilterSize; cnt++)
	{
		FilteredData += ImageROI[tid + cnt] * d_LineFilterH_Signal[cnt];
	}

	if ((CurXPos < ImageWidth) && (CurYPos < ImageHigh))
	{
		d_oImage[CurYPos*ImageWidth + CurXPos] = FilteredData;
	}

}

// seperable convolution vertical
void ImageConv_Seperable_V(unsigned short *d_iImage, unsigned short *d_oImage, int ImageWidth, int ImageHigh, cudaStream_t cstream)
{

	int ImageHigh_Thread = GetImageWidth_Thread(ImageHigh);

	int TotalThreadNum = ImageHigh_Thread * ImageWidth;

	int BlockDim = ThreadsPerBlock;
	int BlockNum = (TotalThreadNum + ThreadsPerBlock - 1) / ThreadsPerBlock;


	gpuImageConv_Seperable_V << < BlockNum, BlockDim, 0, cstream >> > (d_iImage, d_oImage, ImageWidth, ImageHigh);

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


// ImageSubstract, d_ImageC = d_ImageA - d_ImageB
void ImageSubstract(unsigned short *d_ImageA, unsigned short *d_ImageB, unsigned short *d_ImageC, int ImageWidth, int ImageHigh, cudaStream_t cstream)
{

	int TotalThreadNum = ImageWidth * ImageHigh;

	int BlockDim = ThreadsPerBlock;
	int BlockNum = (TotalThreadNum + ThreadsPerBlock - 1) / ThreadsPerBlock;

	gpuImageSubstract << < BlockNum, BlockDim, 0, cstream >> > (d_ImageA, d_ImageB, d_ImageC, ImageWidth, ImageHigh);

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

void ImageVarianceCalcOnSelectRegion(unsigned short *d_RawImg, float *d_MeanDataX, float *d_MeanDataX2, int ImageCenterX, int ImageCenterY, int StdImageSize, int ImageWidth, int ImageHigh, cudaStream_t cstream)
{
	int TotalThreadNum = StdImageSize * StdImageSize;

	int BlockDim = ThreadsPerBlock;
	int BlockNum = (TotalThreadNum + ThreadsPerBlock - 1) / ThreadsPerBlock;

	gpuImageVarianceCalcOnSelectRegion << < BlockNum, BlockDim, 0, cstream >> > (d_RawImg, d_MeanDataX, d_MeanDataX2, ImageCenterX, ImageCenterY, StdImageSize, ImageWidth, ImageHigh);

}

// calculate background intensity
__global__ void gpuImageSharpen(unsigned short *d_RawImg, unsigned short *d_FilteredImg, int ImageWidth, int ImageHigh)
{
	enum {
		SharedMemWidth = ThreadsPerBlock + 2 * SharpFilterSize_Half,

		SharedMemHigh = BlockFilter_BatchLineNum + 2 * SharpFilterSize_Half,

	};

	// side pixels, valid pixels, side pixels
	__shared__ unsigned short ImageROI[SharedMemHigh][SharedMemWidth];


	// init image filter
	float FilterKernel[9]; // 3x3

	for (int cnt = 0; cnt < 9; cnt++)
	{
		FilterKernel[cnt] = -0.125f;
	}
	FilterKernel[4] = 2;

	float(*pFilter)[SharpFilterSize] = (float(*)[SharpFilterSize])FilterKernel;


	int gid = threadIdx.x + blockDim.x*blockIdx.x;
	int tid = threadIdx.x;


	int ImageWidth_Thread = GetImageWidth_Thread(ImageWidth);


	// x,y pos of current thread
	int CurXPos = gid % ImageWidth_Thread;
	int CurYPos = (gid / ImageWidth_Thread)*BlockFilter_BatchLineNum;


	int XPos_GMem = 0;
	int YPos_GMem = 0;

	// load global memory to shared memory, first part
#pragma unroll
	for (int cnt = 0; cnt < SharedMemHigh; cnt++)
	{
		XPos_GMem = CurXPos - SharpFilterSize_Half;
		YPos_GMem = CurYPos - SharpFilterSize_Half + cnt;

		XPos_GMem = max(XPos_GMem, 0);
		YPos_GMem = max(YPos_GMem, 0);

		XPos_GMem = min(XPos_GMem, ImageWidth - 1);
		YPos_GMem = min(YPos_GMem, ImageHigh - 1);

		int Addr_GMem = YPos_GMem*ImageWidth + XPos_GMem;

		ImageROI[cnt][tid] = d_RawImg[Addr_GMem];
	}

	// load global memory to shared memory, second part
#pragma unroll
	for (int cnt = 0; cnt < SharedMemHigh; cnt++)
	{
		if (tid < SharedMemWidth - ThreadsPerBlock)
		{
			XPos_GMem = CurXPos - SharpFilterSize_Half;
			YPos_GMem = CurYPos - SharpFilterSize_Half + cnt;

			// right part
			XPos_GMem = XPos_GMem + ThreadsPerBlock;

			XPos_GMem = max(XPos_GMem, 0);
			YPos_GMem = max(YPos_GMem, 0);

			XPos_GMem = min(XPos_GMem, ImageWidth - 1);
			YPos_GMem = min(YPos_GMem, ImageHigh - 1);

			int Addr_GMem = YPos_GMem*ImageWidth + XPos_GMem;

			ImageROI[cnt][tid + ThreadsPerBlock] = d_RawImg[Addr_GMem];
		}
	}

	__syncthreads();


	// filter part

	for (int cnt = 0; cnt < BlockFilter_BatchLineNum; cnt++)
	{
		int XCenterPos_SMem = tid + SharpFilterSize_Half;
		int YCenterPos_SMem = cnt + SharpFilterSize_Half;

		float FilteredData = ImageROI[YCenterPos_SMem][XCenterPos_SMem];


#pragma unroll
		for (int j = 0; j < SharpFilterSize; j++)
		{
#pragma unroll
			for (int i = 0; i < SharpFilterSize; i++)
			{
				FilteredData += pFilter[j][i] * ImageROI[YCenterPos_SMem - SharpFilterSize_Half + j][XCenterPos_SMem - SharpFilterSize_Half + i];
			}
		}

		FilteredData = max(FilteredData, 0.0f);

		// write filtered data
		XPos_GMem = CurXPos;
		YPos_GMem = CurYPos + cnt;

		if ((XPos_GMem < ImageWidth) && (YPos_GMem < ImageHigh))
		{
			d_FilteredImg[YPos_GMem*ImageWidth + XPos_GMem] = FilteredData;
		}
	}
}


void ImageSharpen(unsigned short *d_RawImg, unsigned short *d_FilteredImg, int ImageWidth, int ImageHigh, cudaStream_t cstream)
{
	int ImageWidth_Thread = GetImageWidth_Thread(ImageWidth);

	int TotalThreadNum = ((ImageHigh + BlockFilter_BatchLineNum - 1) / BlockFilter_BatchLineNum) * ImageWidth_Thread;

	int BlockDim = ThreadsPerBlock;
	int BlockNum = (TotalThreadNum + ThreadsPerBlock - 1) / ThreadsPerBlock;

	gpuImageSharpen << < BlockNum, BlockDim, 0, cstream >> > (d_RawImg, d_RawImg, ImageWidth, ImageHigh);
}

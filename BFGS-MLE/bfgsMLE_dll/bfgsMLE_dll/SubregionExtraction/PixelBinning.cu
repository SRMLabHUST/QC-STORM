#include "PixelBinning.h"



__global__ void gpu_PixelBinning(unsigned short *d_Image, unsigned short *d_oImage, float CameraOffset, int XBin, int YBin, int ImageWidth, int oImageWidth, int oImageHigh);


void PixelBinning_TypeDef::GetPixelBinnedImageForCPU(unsigned short *h_iImage, float CameraOffset, int XBin, int YBin, cudaStream_t cstream)
{

	GetPixelBinnedImageForGPU(h_iImage, CameraOffset, XBin, YBin, cstream);


	cudaMemcpyAsync(h_oImage, d_oImage, oImageWidth*oImageHigh*sizeof(unsigned short), cudaMemcpyDeviceToHost, cstream);
	cudaStreamSynchronize(cstream);

}

void PixelBinning_TypeDef::GetPixelBinnedImageForGPU(unsigned short *h_iImage, float CameraOffset, int XBin, int YBin, cudaStream_t cstream)
{
	oImageWidth = ImageWidth / XBin;
	oImageHigh = ImageHigh / YBin;


	cudaMemcpyAsync(d_Image, h_iImage, ImageWidth*ImageHigh*sizeof(unsigned short), cudaMemcpyHostToDevice, cstream);

	int TotalThreadNum = oImageWidth*oImageHigh;

	int BlockDim = ThreadsPerBlock;
	int BlockNum = ((TotalThreadNum + ThreadsPerBlock - 1) / ThreadsPerBlock);

	gpu_PixelBinning << <BlockNum, BlockDim, 0, cstream >> >(d_Image, d_oImage, CameraOffset, XBin, YBin, ImageWidth, oImageWidth, oImageHigh);


	cudaStreamSynchronize(cstream);

}	


void PixelBinning_TypeDef::Init(int ImageWidth, int ImageHigh)
{
	this->ImageWidth = ImageWidth;
	this->ImageHigh = ImageHigh;


//	cudaMallocHost((void **)&h_Image, ImageWidth*ImageHigh*sizeof(unsigned short));
	cudaMalloc((void **)&d_Image, ImageWidth*ImageHigh*sizeof(unsigned short));

	cudaMallocHost((void **)&h_oImage, ImageWidth*ImageHigh*sizeof(unsigned short));
	cudaMalloc((void **)&d_oImage, ImageWidth*ImageHigh*sizeof(unsigned short));


}

void PixelBinning_TypeDef::Deinit()
{
//	cudaFreeHost(h_Image);
	cudaFree(d_Image);

	cudaFreeHost(h_oImage);
	cudaFree(d_oImage);



}


__global__ void gpu_PixelBinning(unsigned short *d_Image, unsigned short *d_oImage, float CameraOffset, int XBin, int YBin, int ImageWidth, int oImageWidth, int oImageHigh)
{
	int gid = threadIdx.x + blockDim.x*blockIdx.x;

	int TotalThreadNum = oImageWidth*oImageHigh;

	// position in d_oImage
	int XPos = gid % oImageWidth;
	int YPos = gid / oImageWidth;

	int Raw_XPos, Raw_YPos;

	float TotalDat = 0;

	if (gid < TotalThreadNum)
	{
		for (int ycnt = 0; ycnt < YBin; ycnt++)
		{
			Raw_YPos = YPos*YBin + ycnt;

			for (int xcnt = 0; xcnt < XBin; xcnt++)
			{
				Raw_XPos = XPos*XBin + xcnt;

				TotalDat += d_Image[Raw_YPos*ImageWidth + Raw_XPos];
			}

		}

		TotalDat -= (XBin*YBin - 1)*CameraOffset;

		d_oImage[YPos*oImageWidth + XPos] = TotalDat;

	}
}



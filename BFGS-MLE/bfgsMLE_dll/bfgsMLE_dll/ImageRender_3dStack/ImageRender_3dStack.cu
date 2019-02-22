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

#include "ImageRender_3dStack.h"



#define RendROIWidth			7
#define RendROIWidth_half		((int)RendROIWidth/2)

__device__ static void CreateRenderPSF(float oPSFArry[RendROIWidth*RendROIWidth], float RendXPos, float RendYPos, float RendSigma)
{
	float(*RenderPSF)[RendROIWidth] = (float(*)[RendROIWidth])oPSFArry;
	int xcnt, ycnt;

	RendSigma = -0.5f / (RendSigma*RendSigma);

	RendXPos = RendXPos - (int)RendXPos; // shift to 0-0.99
	RendYPos = RendYPos - (int)RendYPos; // shift to 0-0.99
	RendXPos = RendXPos + RendROIWidth_half; // 0 is center
	RendYPos = RendYPos + RendROIWidth_half; // 0 is center

#pragma unroll
	for (ycnt = 0; ycnt < RendROIWidth; ycnt++)
	{
#pragma unroll
		for (xcnt = 0; xcnt < RendROIWidth; xcnt++)
		{
			RenderPSF[ycnt][xcnt] = __expf(((xcnt - RendXPos)*(xcnt - RendXPos) + (ycnt - RendYPos)*(ycnt - RendYPos)) * RendSigma);

		}
	}
}


__global__ void PointFillRegion_Stack3D(float *d_LocArry, float *d_SRIntensityImg, float RenderZDepth, int RGBImageEncodeMode, float FixedlocPrec_x, float SNR_th, float PixelSize, float PixelZoom, int SRImageWidth, int SRImageHigh, int FluoNum)
{

	//    int tid = threadIdx.x;
	int gid = threadIdx.x + blockDim.x*blockIdx.x;

	// only for 2d
	float PSFArry[RendROIWidth*RendROIWidth];
	float(*RenderPSF)[RendROIWidth] = (float(*)[RendROIWidth])PSFArry;

	// parameters for 2d and 3d
	float(*pLocArry)[OutParaNumGS2D];
	pLocArry = (float(*)[OutParaNumGS2D])d_LocArry;


	float PeakPhoton;
	float RendXPos;
	float RendYPos;

	float CurZPos; // z pos in nm

	float CurSNR;

	float LocPrec_x;
	float LocPrec_z;
	float RendSigma_x;
	float RendSigma_z;


	float amp1;

	if (gid < FluoNum)
	{
		PeakPhoton = pLocArry[gid][Pos_PPho]; // intensity as weight 
		RendXPos = pLocArry[gid][Pos_XPos] * PixelZoom;
		RendYPos = pLocArry[gid][Pos_YPos] * PixelZoom;

		CurZPos = pLocArry[gid][Pos_ZPos];

		CurSNR = pLocArry[gid][Pos_PSNR];
		LocPrec_x = pLocArry[gid][Pos_CrbX]; // unit is nm
											 //		LocPrec_z = pLocArry[gid][Pos_CrbZ]; // unit is nm
		LocPrec_z = LocPrec_x * 2;

		// valid fluos
		if ((CurSNR >= SNR_th) && (PeakPhoton > 1.0f) && (RendXPos > RendROIWidth) && (RendXPos < SRImageWidth - RendROIWidth) && (RendYPos > RendROIWidth) && (RendYPos < SRImageHigh - RendROIWidth))
		{
			if ((RGBImageEncodeMode == 1) || (RGBImageEncodeMode == 2))
			{
				// use the same weight but not the peak photon
				PeakPhoton = 1.0f;
			}
			if (RGBImageEncodeMode == 2)
			{
				LocPrec_x = FixedlocPrec_x;
			}

			// rendered psf don't be too small
			LocPrec_x = max(LocPrec_x, PixelSize / PixelZoom / 2);

			LocPrec_z = LocPrec_x * 2;

			// calculate rendering PSF, 
			RendSigma_x = LocPrec_x / (PixelSize / PixelZoom); //unit: sr image pixel
			RendSigma_z = LocPrec_z; // unit: nm


			float ZPosDiff = abs(CurZPos - RenderZDepth);
			float ZAmpDecay = __expf(-ZPosDiff*ZPosDiff / (2 * RendSigma_z*RendSigma_z));

			PeakPhoton = PeakPhoton*ZAmpDecay;

			CreateRenderPSF(PSFArry, RendXPos, RendYPos, RendSigma_x);


#pragma unroll
			for (int rcnt = 0; rcnt < RendROIWidth; rcnt++)
			{
				int Offset = SRImageWidth*((int)RendYPos - RendROIWidth_half + rcnt) + (int)RendXPos - RendROIWidth_half;

#pragma unroll
				for (int ccnt = 0; ccnt < RendROIWidth; ccnt++)
				{
					atomicAdd(&d_SRIntensityImg[Offset + ccnt], PeakPhoton*RenderPSF[rcnt][ccnt]);

				}

			}
		}
	}
}


///////////////
// cuda wrapper


void ImageRender3DStackData_TypeDef::GetSaveImgTop(float *ih_LocArry, LocalizationPara & LocPara, float RenderZDepth, int RGBImageEncodeMode, float FixedlocPrec_x, int TotalFluoNum, cudaStream_t cstream)
{

	cudaMemcpyAsync(d_LocArry, ih_LocArry, TotalFluoNum*OutParaNumGS2D * sizeof(float), cudaMemcpyHostToDevice, cstream);

	cudaMemsetAsync(d_SRIntensityImg, 0, LocPara.SRImageWidth*LocPara.SRImageHigh * sizeof(float), cstream);

	int BlockDim = ThreadsPerBlock;
	int BlockNum = ((TotalFluoNum + ThreadsPerBlock - 1) / ThreadsPerBlock);


	PointFillRegion_Stack3D << <BlockNum, BlockDim, 0, cstream >> > (d_LocArry, d_SRIntensityImg, RenderZDepth, RGBImageEncodeMode, FixedlocPrec_x, LocPara.SNR_th, LocPara.PixelSize, LocPara.PixelZoom, LocPara.SRImageWidth, LocPara.SRImageHigh, TotalFluoNum);

	cudaMemcpyAsync(h_SRIntensityImg, d_SRIntensityImg, LocPara.SRImageWidth*LocPara.SRImageHigh * sizeof(float), cudaMemcpyDeviceToHost, cstream);
	cudaStreamSynchronize(cstream);

}



void ImageRender3DStackData_TypeDef::Init(int SRImageWidth, int SRImageHigh, int TotalFluoNum)
{
	// host and gpu
	cudaError_t err;

	err = cudaMallocHost((void **)&h_LocArry, TotalFluoNum*OutParaNumGS2D * sizeof(float));
	err = cudaMalloc((void **)&d_LocArry, TotalFluoNum*OutParaNumGS2D * sizeof(float));

	err = cudaMallocHost((void **)&h_SRIntensityImg, SRImageWidth*SRImageHigh * sizeof(float));
	err = cudaMalloc((void **)&d_SRIntensityImg, SRImageWidth*SRImageHigh * sizeof(float));


}

void ImageRender3DStackData_TypeDef::Deinit()
{

	cudaFreeHost(h_LocArry);
	cudaFree(d_LocArry);

	cudaFreeHost(h_SRIntensityImg);
	cudaFree(d_SRIntensityImg);

}


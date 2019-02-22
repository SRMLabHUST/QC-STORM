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

#include "SubregionExtraction_ld.h"



#define LineBlockNum			(ThreadsPerBlock)
#define LineBlockNum_7x7		(LineBlockNum + 7) // LineBlockNum +6, here to avoid bank conflict
#define ExtPixelNum_7x7			6

#define LineBlockNum_3x3		(LineBlockNum + 3) // LineBlockNum +6, here to avoid bank conflict
#define ExtPixelNum_3x3			2


#define BatchLineProcNum		ThreadsPerBlock




//ImageSource = 0: pRawImg is a cpu pointer. 1:pRawImg is a GPU pointer
void LDROIExtractData_TypeDef::ExtractMolecules(unsigned short *pImgData, int ImageSource, LocalizationPara & LocPara, int StartFrame, int BatchFrameNum, cudaStream_t cstream)
{
	//ImageSource = 0: pRawImg is a cpu pointer. 1:pRawImg is a GPU pointer

	int GroupImgHigh = BatchFrameNum*LocPara.ImageHigh;
	int GroupImgSize = BatchFrameNum*LocPara.ImageWidth*LocPara.ImageHigh*sizeof(short);

	const int EachGroupFluoNum = MaxPointNum / BatchFrameNum;

	int ROIDataLen = LocPara.ROISize*(LocPara.ROISize + 1);

	int RegionAddrOffset = (h_TotalRegionCount)*ROIDataLen;

	cudaMemsetAsync(d_RegionNum, 0, BatchFrameNum*sizeof(int), cstream);

//	printf("batch num:%d %d %d\n", BatchFrameNum, ImageWidth, ImageHigh);
	// get raw image

	if (ImageSource == ImageSource_CPU_Pinned)
	{
		// image from CPU
		cudaMemcpyAsync(d_RawImg, pImgData, GroupImgSize, cudaMemcpyHostToDevice, cstream); // h_RawImg
	}
	else if (ImageSource == ImageSource_CPU_Normal)
	{
		// image from CPU
		cudaMemcpy(d_RawImg, pImgData, GroupImgSize, cudaMemcpyHostToDevice); // h_RawImg

	}
	else if (ImageSource == ImageSource_GPU)
	{
		// image from GPU
		d_RawImg = pImgData;
	}
	else
	{
		// image data error
		return;
	}


	// image filtering 
	// filtering algorithm could be improved futher, this algorithm performance rely on the image hight

	int ImgWidthBlockNum = (LocPara.ImageWidth + ValidProcCore_7x7 - 1) / ValidProcCore_7x7;
	int ImgWidthThreadNum = ImgWidthBlockNum * 32;
	int ImgHighBlockNum = (GroupImgHigh + BatchLineNum - 1) / BatchLineNum;

	int TotalThreadNum = ImgHighBlockNum*ImgWidthThreadNum;

	int BlockDim = ThreadsPerBlock;
	int BlockNum = TotalThreadNum / ThreadsPerBlock;

	gpuImgFilleringLD1 << < BlockNum, BlockDim, 0, cstream >> >(d_RawImg, d_GaussImg, d_StdImg, LocPara.ImageWidth, GroupImgHigh);
	gpuImgFilleringLD2 << < BlockNum, BlockDim, 0, cstream >> >(d_GaussImg, d_AnnuImg, LocPara.ImageWidth, GroupImgHigh);

	ImgWidthBlockNum = (LocPara.ImageWidth + ValidProcCore_3x3 - 1) / ValidProcCore_3x3;
	ImgWidthThreadNum = ImgWidthBlockNum * 32;

	TotalThreadNum = ImgHighBlockNum*ImgWidthThreadNum;

	BlockNum = TotalThreadNum / ThreadsPerBlock;
	gpuSubregionFindingLD << < BlockNum, BlockDim, 0, cstream >> >(d_AnnuImg, d_StdImg, d_RegionPosMem, d_RegionNum, LocPara.ImageWidth, GroupImgHigh, LocPara.ImageHigh);
	
	cudaMemcpyAsync(h_RegionNum, d_RegionNum, BatchFrameNum * sizeof(int), cudaMemcpyDeviceToHost, cstream);
//	cudaMemcpyAsync(h_RegionPosMem, d_RegionPosMem, MaxPointNum * RegionPosInfNum * sizeof(int), cudaMemcpyDeviceToHost, cstream);
	cudaStreamSynchronize(cstream);

	int CurRegionNum = 0;
	int TotalRegionNum = 0;

	int cnt = 0;

	int pcnt = 0;


	TotalRegionNum = 0;

	for (cnt = 0; cnt < BatchFrameNum; cnt++)
	{
		// since whole memory is seperated for batch frames, and there are zero gaps, see gpuSubregionFindingLD
		CurRegionNum = h_RegionNum[cnt];


		BlockNum = ((CurRegionNum + ThreadsPerBlock - 1) / ThreadsPerBlock);

		// region extraction with point position
		gpuSubregionExtractionLD << < BlockNum, BlockDim, 0, cstream >> >(d_RawImg, &d_RegionMem[TotalRegionNum*ROIDataLen], &d_RegionPosMem[cnt*EachGroupFluoNum*RegionPosInfNum], CurRegionNum, LocPara.ROISize, LocPara.ImageWidth, GroupImgHigh, LocPara.ImageHigh, StartFrame);
		
		TotalRegionNum += CurRegionNum;
	}


	cudaMemcpyAsync(&h_RegionMem[RegionAddrOffset], d_RegionMem, TotalRegionNum * ROIDataLen * sizeof(short), cudaMemcpyDeviceToHost, cstream);
	cudaStreamSynchronize(cstream); // wait task of this stream finish


	h_TotalRegionCount = h_TotalRegionCount + TotalRegionNum;


//	printf("total fluo:%d \n", h_TotalRegionCount);
}



void LDROIExtractData_TypeDef::Init(LocalizationPara & LocPara)
{
	//		cudaError_t err;
	cudaError_t err;

	int ROIDataLen = LocPara.ROISize*(LocPara.ROISize + 1);

	int MaxBatchImgNum = (2048 * 2048 / LocPara.ImageWidth / LocPara.ImageHigh) + 2;
	if (MaxBatchImgNum < 1)MaxBatchImgNum = 1;

	// host and gpu
	err = cudaMallocHost((void **)&h_RawImg, MaxBatchedImageSize);
	HandleErr(err, "cudaMalloc h_RawImg");
	err = cudaMalloc((void **)&d_RawImg, MaxBatchedImageSize);
	HandleErr(err, "cudaMalloc d_RawImg");

	cudaMalloc((void **)&d_GaussImg, MaxBatchedImageSize);
	cudaMalloc((void **)&d_AnnuImg, MaxBatchedImageSize);
	cudaMalloc((void **)&d_StdImg, MaxBatchedImageSize);

	cudaMallocHost((void **)&h_RegionMem, MaxPointNum * ROIDataLen * sizeof(unsigned short));

	cudaMalloc((void **)&d_RegionMem, MaxPointNum * ROIDataLen * sizeof(unsigned short));

	cudaMallocHost((void **)&h_RegionNum, MaxBatchImgNum *sizeof(int));

	cudaMalloc((void **)&d_RegionNum, MaxBatchImgNum *sizeof(int));
	cudaMalloc((void **)&d_RegionPosMem, MaxPointNum * RegionPosInfNum * sizeof(int));

	
	h_TotalRegionCount = 0;

}


void LDROIExtractData_TypeDef::Deinit()
{
	// host and gpu
	cudaFreeHost(h_RawImg);
	cudaFree(d_RawImg);

	cudaFree(d_GaussImg);
	cudaFree(d_AnnuImg);
	cudaFree(d_StdImg);

	cudaFreeHost(h_RegionMem);

	cudaFree(d_RegionMem);

	cudaFreeHost(h_RegionNum);
//	cudaFreeHost(h_RegionPosMem);

	cudaFree(d_RegionNum);
	cudaFree(d_RegionPosMem);

}

int LDROIExtractData_TypeDef::GetRegionNum()
{
	return h_TotalRegionCount;
}

void LDROIExtractData_TypeDef::ResetRegionNum()
{
	h_TotalRegionCount = 0;
}

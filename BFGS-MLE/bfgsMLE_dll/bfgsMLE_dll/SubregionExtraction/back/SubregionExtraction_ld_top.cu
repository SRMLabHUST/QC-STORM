#include "SubregionExtraction_ld.h"



#define LineBlockNum			(ThreadsPerBlock)
#define LineBlockNum_7x7		(LineBlockNum + 7) // LineBlockNum +6, here to avoid bank conflict
#define ExtPixelNum_7x7			6

#define LineBlockNum_3x3		(LineBlockNum + 3) // LineBlockNum +6, here to avoid bank conflict
#define ExtPixelNum_3x3			2


#define BatchLineProcNum		ThreadsPerBlock

#define _AFXDLL
#include <afx.h>



//ImageSource = 0: pRawImg is a cpu pointer. 1:pRawImg is a GPU pointer
void LDROIExtractData_TypeDef::ExtractMolecules(unsigned short *pImgData, int ImageSource, LocalizationPara & LocPara, int StartFrame, int BatchFrameNum, cudaStream_t cstream)
{
	//ImageSource = 0: pRawImg is a cpu pointer. 1:pRawImg is a GPU pointer

	int GroupImgHigh = BatchFrameNum*LocPara.ImageHigh;
	int GroupImgSize = BatchFrameNum*LocPara.ImageWidth*LocPara.ImageHigh*sizeof(unsigned short);

	const int EachGroupFluoNum = MaxPointNum / BatchFrameNum;

	int SubRegionWidth = LocPara.ROISize*(LocPara.ROISize + 1);

	int RegionAddrOffset = (h_TotalRegionCount)*SubRegionWidth;

	cudaMemsetAsync(d_RegionNum, 0, BatchFrameNum*sizeof(int), cstream);

//	printf("batch num:%d %d %d\n", BatchFrameNum, ImageWidth, ImageHigh);
	// get raw image

	if (ImageSource == ImageSource_CPU)
	{
		// image from CPU
		cudaMemcpyAsync(d_RawImg, pImgData, GroupImgSize, cudaMemcpyHostToDevice, cstream); // h_RawImg
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

	int BlockSize = TotalThreadNum / ThreadsPerBlock;
	gpuImgFilleringLD1 << < BlockSize, ThreadsPerBlock, 0, cstream >> >(d_RawImg, d_GaussImg, d_StdImg, LocPara.ImageWidth, GroupImgHigh);
	gpuImgFilleringLD2 << < BlockSize, ThreadsPerBlock, 0, cstream >> >(d_GaussImg, d_AnnuImg, LocPara.ImageWidth, GroupImgHigh);

	ImgWidthBlockNum = (LocPara.ImageWidth + ValidProcCore_3x3 - 1) / ValidProcCore_3x3;
	ImgWidthThreadNum = ImgWidthBlockNum * 32;

	TotalThreadNum = ImgHighBlockNum*ImgWidthThreadNum;

	BlockSize = TotalThreadNum / ThreadsPerBlock;
	gpuSubregionFindingLD << < BlockSize, ThreadsPerBlock, 0, cstream >> >(d_AnnuImg, d_StdImg, d_RegionPosMem, d_RegionNum, LocPara.ImageWidth, GroupImgHigh, LocPara.ImageHigh);
	
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


		BlockSize = ((CurRegionNum + ThreadsPerBlock - 1) / ThreadsPerBlock);

		// region extraction with point position
		gpuSubregionExtractionLD << < BlockSize, ThreadsPerBlock, 0, cstream >> >(d_RawImg, &d_RegionMem[TotalRegionNum*SubRegionWidth], &d_RegionPosMem[cnt*EachGroupFluoNum*RegionPosInfNum], CurRegionNum, LocPara.ROISize, LocPara.ImageWidth, GroupImgHigh, LocPara.ImageHigh, StartFrame);
		
		TotalRegionNum += CurRegionNum;

	}


	cudaMemcpyAsync(&h_RegionMem[RegionAddrOffset], d_RegionMem, TotalRegionNum * SubRegionWidth * sizeof(unsigned short), cudaMemcpyDeviceToHost, cstream);
	cudaStreamSynchronize(cstream); // wait task of this stream finish


	h_TotalRegionCount = h_TotalRegionCount + TotalRegionNum;


//	printf("total fluo:%d \n", h_TotalRegionCount);
}



void LDROIExtractData_TypeDef::Init(LocalizationPara & LocPara)
{
	//		cudaError_t err;
	cudaError_t err;

//	int ImgSize = ImageWidth*ImageHigh*sizeof(unsigned short);
	int SubRegionWidth = LocPara.ROISize*(LocPara.ROISize + 1);

	int BatchImgNum = (2048 * 2048 / LocPara.ImageWidth / LocPara.ImageHigh) + 2;
	if (BatchImgNum < 1)BatchImgNum = 1;


	// malloc resources for the pointers in the structure

	// host and gpu
	err = cudaMallocHost((void **)&h_RawImg, MaxGroupImgSize);
	HandleErr(err, "cudaMalloc d_SRIntensityImg");
	err = cudaMalloc((void **)&d_RawImg, MaxGroupImgSize);
	HandleErr(err, "cudaMalloc d_SRIntensityImg");

	cudaMalloc((void **)&d_GaussImg, MaxGroupImgSize);
	cudaMalloc((void **)&d_AnnuImg, MaxGroupImgSize);
	cudaMalloc((void **)&d_StdImg, MaxGroupImgSize);

	cudaMallocHost((void **)&h_RegionMem, MaxPointNum * SubRegionWidth * sizeof(unsigned short));

	cudaMalloc((void **)&d_RegionMem, MaxPointNum * SubRegionWidth * sizeof(unsigned short));

	cudaMallocHost((void **)&h_RegionNum, BatchImgNum*sizeof(int));
//	cudaMallocHost((void **)&h_RegionPosMem, MaxPointNum * RegionPosInfNum * sizeof(int));

	cudaMalloc((void **)&d_RegionNum, BatchImgNum*sizeof(int));
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

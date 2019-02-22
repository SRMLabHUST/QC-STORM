#include "SubregionExtraction_hd.h"


void GetClusterPos(unsigned short PosArry[49], int *xoffset, int *yoffset);
char IsValidClusterPos(unsigned short PosArry[49]);
void cpuPointCluster(HDROIExtractData_TypeDef * HDPara, int ImageWidth, int ImageHigh);


///////////////////////////

void HDROIExtract_ExtractMolecules(HDROIExtractData_TypeDef *h_HDROIExtractData, int LargeImgWidth, int LargeImgHigh, int ImgHighAdj, int RegionSize, int FrameOffset, cudaStream_t cstream)
{
	cudaError_t err;

	int ImgSize = LargeImgWidth*LargeImgHigh*sizeof(unsigned short);
	int WholeImageWidth = RegionSize*(RegionSize + 1);

	err=cudaMemcpyAsync(h_HDROIExtractData->d_RawImg, h_HDROIExtractData->h_RawImg, ImgSize, cudaMemcpyHostToDevice, cstream);

	// image filtering
	int GridSize = ((LargeImgHigh + ThreadsPerBlock - 1) / ThreadsPerBlock);

	gpuImgFiller_gauss_std << < GridSize, ThreadsPerBlock, 0, cstream >> >(h_HDROIExtractData->d_RawImg, h_HDROIExtractData->d_GaussImg, h_HDROIExtractData->d_StdImg1, LargeImgWidth, LargeImgHigh);
	gpuImgFiller_annu_stdannu << < GridSize, ThreadsPerBlock, 0, cstream >> >(h_HDROIExtractData->d_GaussImg, h_HDROIExtractData->d_StdImg1, h_HDROIExtractData->d_AnnuImg, h_HDROIExtractData->d_StdImg2, LargeImgWidth, LargeImgHigh);
	
	cudaStreamQuery(cstream);

	// debug
//	cudaMemcpyAsync(h_HDROIExtractData->h_PointImg, h_HDROIExtractData->d_StdImg2, ImgSize, cudaMemcpyDeviceToHost, cstream);
//	cudaStreamSynchronize(cstream); // wait task of this stream finish


	// sub-region finding
	cudaMemsetAsync(h_HDROIExtractData->d_PointImg, 0, LargeImgWidth*LargeImgHigh*sizeof(unsigned short), cstream);
	cudaMemsetAsync(h_HDROIExtractData->d_RegionNum, 0, sizeof(int), cstream);

	cudaStreamSynchronize(cstream); // wait task of this stream finish

	gpuSubregionFinding << < GridSize, ThreadsPerBlock, 0, cstream >> >(h_HDROIExtractData->d_AnnuImg, h_HDROIExtractData->d_StdImg2, h_HDROIExtractData->d_PointImg, h_HDROIExtractData->d_RegionNum, h_HDROIExtractData->d_RegionMem, LargeImgWidth, LargeImgHigh);


	cudaMemcpyAsync(h_HDROIExtractData->h_PointImg, h_HDROIExtractData->d_PointImg, ImgSize, cudaMemcpyDeviceToHost, cstream);
	cudaMemcpyAsync(h_HDROIExtractData->h_iPosMem, h_HDROIExtractData->d_RegionMem, MaxPointNum * 3 * sizeof(int), cudaMemcpyDeviceToHost, cstream);
	cudaMemcpyAsync(h_HDROIExtractData->h_iPosCount, h_HDROIExtractData->d_RegionNum, sizeof(int), cudaMemcpyDeviceToHost, cstream);


	cudaStreamSynchronize(cstream); // wait task of this stream finish

	// cluster the subregion pos by CPU
	cpuPointCluster(h_HDROIExtractData->h_PointImg, h_HDROIExtractData->h_iPosCount, h_HDROIExtractData->h_iPosMem, h_HDROIExtractData->h_cPointImg, h_HDROIExtractData->h_oPosCount, h_HDROIExtractData->h_oPosMem, LargeImgWidth, LargeImgHigh);

//	printf("region count:%d %d\n", *h_HDROIExtractData->h_iPosCount,*h_HDROIExtractData->h_oPosCount);

	cudaMemcpyAsync(h_HDROIExtractData->d_RegionMem, h_HDROIExtractData->h_oPosMem, MaxPointNum * 3 * sizeof(int), cudaMemcpyHostToDevice, cstream);
	cudaMemcpyAsync(h_HDROIExtractData->d_RegionNum, h_HDROIExtractData->h_oPosCount, sizeof(int), cudaMemcpyHostToDevice, cstream);

	GridSize = ((*h_HDROIExtractData->h_oPosCount + ThreadsPerBlock - 1) / ThreadsPerBlock);

	gpuSubregionExtraction << < GridSize, ThreadsPerBlock, 0, cstream >> >(h_HDROIExtractData->d_RawImg, h_HDROIExtractData->d_RegionNum, h_HDROIExtractData->d_RegionMem, h_HDROIExtractData->d_RegionMem, RegionSize, LargeImgWidth, LargeImgHigh, ImgHighAdj, FrameOffset);

	cudaMemcpyAsync(h_HDROIExtractData->h_RegionMem + ((*h_HDROIExtractData->h_TotalRegionNum)*WholeImageWidth), h_HDROIExtractData->d_RegionMem, (*h_HDROIExtractData->h_oPosCount) * WholeImageWidth*sizeof(unsigned short), cudaMemcpyDeviceToHost, cstream);
	*h_HDROIExtractData->h_TotalRegionNum = *h_HDROIExtractData->h_TotalRegionNum + *h_HDROIExtractData->h_oPosCount;

	cudaStreamSynchronize(cstream); // wait task of this stream finish
}


void HDROIExtract_Init(HDROIExtractData_TypeDef **ha_HDROIExtractData, int LargeImgWidth, int LargeImgHigh, int RegionSize)
{

	int ImgSize = LargeImgWidth*LargeImgHigh*sizeof(unsigned short);
	int WholeImageWidth = RegionSize*(RegionSize + 1);

	cudaMallocHost((void **)ha_HDROIExtractData, sizeof(HDROIExtractData_TypeDef)); // raw image in CPU


	HDROIExtractData_TypeDef *h_HDROIExtractData = *ha_HDROIExtractData;

	// malloc resources for the pointers in the structure
	// host
	cudaMallocHost((void **)&h_HDROIExtractData->h_RawImg, ImgSize); // raw image in CPU

	// gpu
	cudaMalloc((void **)&h_HDROIExtractData->d_RawImg, ImgSize);  //  raw image in GPU


	cudaMalloc((void **)&h_HDROIExtractData->d_GaussImg, ImgSize);  // dog filtered image in GPU
	cudaMalloc((void **)&h_HDROIExtractData->d_AnnuImg, ImgSize);  // dog filtered image in GPU
	cudaMalloc((void **)&h_HDROIExtractData->d_StdImg1, ImgSize); // std image
	cudaMalloc((void **)&h_HDROIExtractData->d_StdImg2, ImgSize); // annular smooth filtered std image

	cudaMallocHost((void **)&h_HDROIExtractData->h_iPosMem, MaxPointNum * 3 * sizeof(int)); // memory for raw pos
	cudaMallocHost((void **)&h_HDROIExtractData->h_iPosCount, sizeof(int));  				// memory for raw pos number
	cudaMallocHost((void **)&h_HDROIExtractData->h_oPosMem, MaxPointNum * 3 * sizeof(int)); // memory for clustered pos by CPU
	cudaMallocHost((void **)&h_HDROIExtractData->h_oPosCount, sizeof(int)); 				// memory for clustered pos number by CPU

	cudaMalloc((void **)&h_HDROIExtractData->d_RegionMem, MaxPointNum * 3 * sizeof(int));		// pos memory in GPU
	cudaMalloc((void **)&h_HDROIExtractData->d_RegionNum, sizeof(int));


	// point images for host and device useage
	cudaMallocHost((void **)&h_HDROIExtractData->h_PointImg, ImgSize);  // point marked image in CPU
	cudaMallocHost((void **)&h_HDROIExtractData->h_cPointImg, ImgSize); // clustered point  marked image in CPU
	cudaMalloc((void **)&h_HDROIExtractData->d_PointImg, ImgSize);


	cudaMallocHost((void **)&h_HDROIExtractData->h_RegionMem, MaxPointNum * WholeImageWidth * sizeof(unsigned short)); // memory for sub region in CPU
	cudaMallocHost((void **)&h_HDROIExtractData->h_TotalRegionNum, sizeof(int)); // sub region number

	// memory for sub region in GPU
	cudaMalloc((void **)&h_HDROIExtractData->d_RegionMem, MaxPointNum * WholeImageWidth * sizeof(unsigned short));


	HDROIExtract_ResetRegionNum(h_HDROIExtractData);

}


void HDROIExtract_Deinit(HDROIExtractData_TypeDef *h_HDROIExtractData)
{
	cudaFreeHost(h_HDROIExtractData->h_RawImg); // raw image in CPU

	cudaFree(h_HDROIExtractData->d_RawImg);  //  raw image in GPU

	cudaFree(h_HDROIExtractData->d_GaussImg);  // dog filtered image in GPU
	cudaFree(h_HDROIExtractData->d_AnnuImg);  // dog filtered image in GPU
	cudaFree(h_HDROIExtractData->d_StdImg1); // std image
	cudaFree(h_HDROIExtractData->d_StdImg2); // annular smooth filtered std image

	cudaFreeHost(h_HDROIExtractData->h_iPosMem);   // memory for raw pos
	cudaFreeHost(h_HDROIExtractData->h_iPosCount); // memory for raw pos number
	cudaFreeHost(h_HDROIExtractData->h_oPosMem);   // memory for clustered pos by CPU
	cudaFreeHost(h_HDROIExtractData->h_oPosCount); // memory for clustered pos number by CPU

	cudaFree(h_HDROIExtractData->d_RegionMem); 		// pos memory in GPU
	cudaFree(h_HDROIExtractData->d_RegionNum);

	cudaFreeHost(h_HDROIExtractData->h_PointImg); 
	cudaFreeHost(h_HDROIExtractData->h_cPointImg);
	cudaFree(h_HDROIExtractData->d_PointImg);

	cudaFreeHost(h_HDROIExtractData->h_RegionMem);
	cudaFreeHost(h_HDROIExtractData->h_TotalRegionNum); 

	cudaFree(h_HDROIExtractData->d_RegionMem); 

	// free
	cudaFreeHost(h_HDROIExtractData);

}

int HDROIExtract_GetRegionNum(HDROIExtractData_TypeDef * h_HDROIExtractData)
{
	return *h_HDROIExtractData->h_TotalRegionNum;
}

void HDROIExtract_ResetRegionNum(HDROIExtractData_TypeDef * h_HDROIExtractData)
{
	*h_HDROIExtractData->h_TotalRegionNum = 0;
}


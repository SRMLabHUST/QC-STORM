#include "bfgs_hd_top.h"

#include <time.h>


int HDLoc_BFGS_MLELocalization(unsigned short * h_SubRegion, HDLocData_TypeDef* h_HDLocData, int FluoNum, float Offset, float kadc, float PsfWidth, int isFPGAProc, cudaStream_t cstream)
{
	// , HDLocData_TypeDef* d_HDLocPara
	cudaError_t err;
	// FluoNum must be the integer multiples of 32
	int suboffset = 0;
	//	suboffset = GetSubRegionOffset_hd(h_HDLocData->h_SubRegion, ROISize);

	// #if bfgs_debug
	// 	printf("sub Offset:%d\n", suboffset);
	// #endif

	if (PsfWidth < 0.2f)PsfWidth = 1.10f;

	PsfWidth = (1.0f / (2.0f*PsfWidth*PsfWidth));

	//	 ROISize = 9;
	const int RegionDataSize = (9*(9 + 1)); // 56 = 7*8

	int Point1Num; // number of single point region
	int Point2Num; // number of two    point region
	int Point3Num; // number of three  point region


	cudaMemsetAsync(h_HDLocData->d_Point1Num, 0, sizeof(int), cstream);
	cudaMemsetAsync(h_HDLocData->d_Point2Num, 0, sizeof(int), cstream);
	cudaMemsetAsync(h_HDLocData->d_Point3Num, 0, sizeof(int), cstream);

	if (suboffset>0)
	{
		FluoNum--;
	}
	cudaMemcpyAsync(h_HDLocData->d_SubRegion, &h_SubRegion[suboffset], FluoNum*RegionDataSize*sizeof(short), cudaMemcpyHostToDevice, cstream);

	int BlockSize = ((FluoNum + ThreadsPerBlock - 1) / ThreadsPerBlock);
	// clasify the subregion into different category:1,2,3 point in a region
	GPUSubregionSeperate << <BlockSize, ThreadsPerBlock, 0, cstream >> >(h_HDLocData->d_SubRegion, FluoNum, h_HDLocData->d_Point1Num, h_HDLocData->d_Point1PosArry, h_HDLocData->d_Point2Num, h_HDLocData->d_Point2PosArry, h_HDLocData->d_Point3Num, h_HDLocData->d_Point3PosArry);


	cudaMemcpyAsync(&Point1Num, h_HDLocData->d_Point1Num, sizeof(int), cudaMemcpyDeviceToHost, cstream);
	cudaMemcpyAsync(&Point2Num, h_HDLocData->d_Point2Num, sizeof(int), cudaMemcpyDeviceToHost, cstream);
	cudaMemcpyAsync(&Point3Num, h_HDLocData->d_Point3Num, sizeof(int), cudaMemcpyDeviceToHost, cstream);

	cudaStreamSynchronize(cstream);

	// total point num
	*h_HDLocData->h_PointNum = Point1Num + Point2Num*2 + Point3Num*3;

//	printf("point num: 1,2,3:%d + %d + %d = %d\n", Point1Num, Point2Num, Point3Num, FluoNum);

	int addrOffset=0;

//	int time1, time2;
//	float p1time, p2time, p3time;

//	time1 = clock();

	// 1 point loc
	int GridSize1 = ((Point1Num + ThreadsPerBlock - 1) / ThreadsPerBlock);
	MLEROILocTop_1p << <GridSize1, ThreadsPerBlock, 0, cstream >> >(h_HDLocData->d_SubRegion, h_HDLocData->d_LocArry1, h_HDLocData->d_Point1PosArry, Offset, kadc, Point1Num, isFPGAProc);

	addrOffset=0;
	cudaMemcpyAsync(h_HDLocData->h_LocArry + addrOffset, h_HDLocData->d_LocArry1, Point1Num*OutParaNum1*sizeof(float), cudaMemcpyDeviceToHost, cstream);

//	cudaStreamSynchronize(cstream);
//	time2 = clock();
//	p1time = time2 - time1 + 0.1;
//	printf("p1 finished %f\n", p1time);

//	time1 = clock();
	// 2 point loc
	int GridSize2 = ((Point2Num + ThreadsPerBlock - 1) / ThreadsPerBlock);
	MLEROILocTop_2p << <GridSize2, ThreadsPerBlock, 0, cstream >> >(h_HDLocData->d_SubRegion, h_HDLocData->d_LocArry2, h_HDLocData->d_Point2PosArry, Offset, kadc, PsfWidth, Point2Num, isFPGAProc);

	addrOffset = PointParaNum*Point1Num;
	cudaMemcpyAsync(h_HDLocData->h_LocArry + addrOffset, h_HDLocData->d_LocArry2, Point2Num*OutParaNum2*sizeof(float), cudaMemcpyDeviceToHost, cstream);

//	cudaStreamSynchronize(cstream);
//	time2 = clock();
//	p2time = time2 - time1 + 0.1;
//	printf("p2 finished %f\n", p2time);

//	time1 = clock();
	// 3 point loc
	int GridSize3 = ((Point3Num + ThreadsPerBlock - 1) / ThreadsPerBlock);
	MLEROILocTop_3p << <GridSize3, ThreadsPerBlock, 0, cstream >> >(h_HDLocData->d_SubRegion, h_HDLocData->d_LocArry3, h_HDLocData->d_Point3PosArry, Offset, kadc, PsfWidth, Point3Num, isFPGAProc);
	
	addrOffset = PointParaNum*Point1Num + PointParaNum*Point2Num*2;
	cudaMemcpyAsync(h_HDLocData->h_LocArry + addrOffset, h_HDLocData->d_LocArry3, Point3Num*OutParaNum3*sizeof(float), cudaMemcpyDeviceToHost, cstream);
	
	cudaStreamQuery(cstream);

//	cudaStreamSynchronize(cstream);
//	time2 = clock();
//	p3time = time2 - time1 + 0.1;
//	printf("p3 finished %f\n", p3time);

//	printf("loc speed 1-3: %d %d %d\n", (int)(Point1Num * 1000 / p1time), (int)(Point2Num * 1000 / p2time), (int)(Point3Num * 1000 / p3time));
	return 0;
}


// get tha data Offset of received USB3.0 data, since the data Offset maybe unknown
int HDLoc_GetSubRegionOffset(unsigned short * FluoPoint, int ROISize)
{
	int SubRegionOffset = 0;
	int cnt = 0;

	int SensorID;
	SensorID = FluoPoint[54] & 0x03;

	if ((FluoPoint[49] == 0) && (FluoPoint[55] == 0))
	{
		if ((SensorID == 1) || (SensorID == 2))
		{
			SubRegionOffset = 0;
			return SubRegionOffset;
		}
	}
	// if the Offset is not 0, then find it
	for (cnt = 0; cnt < 56; cnt++)
	{
		if ((FluoPoint[cnt] == 0) && (FluoPoint[cnt + 6] == 0))
		{
			SensorID = FluoPoint[cnt + 5] & 0x03;
			if ((SensorID == 1) || (SensorID == 2))
			{
				SubRegionOffset = cnt + 7;
				break;
			}
		}
	}
	return SubRegionOffset;
}

void HDLoc_GetImgSizeFromRegion(unsigned short * FluoPoint, int ROISize, int *ImageWidth, int *ImageHigh)
{
	int SubRegionOffset = 0;
	SubRegionOffset = HDLoc_GetSubRegionOffset(FluoPoint, ROISize);

	int InfOffset = SubRegionOffset + ROISize*ROISize;
	// Offset lead by region size
	int RegionOffset = (ROISize - 5) / 2;


	*ImageHigh = (FluoPoint[InfOffset + RegionOffset + 3]) * 2;	// image high
	*ImageWidth = (FluoPoint[InfOffset + RegionOffset + 4] / 4) * 5;	// image width

}

int HDLoc_GetFirstFrameFromRegion(unsigned short * FluoPoint, int ROISize)
{
	int SubRegionOffset = 0;
	SubRegionOffset = HDLoc_GetSubRegionOffset(FluoPoint, ROISize);

	int InfOffset = SubRegionOffset + ROISize*ROISize;
	// Offset lead by region size
	int RegionOffset = (ROISize - 5) / 2;

	unsigned int CurFrameNum;

	CurFrameNum = FluoPoint[InfOffset + RegionOffset + 2];			// cur image frame

	return CurFrameNum;
}


void HDLoc_Init(HDLocData_TypeDef** ha_HDLocData, int ROISize)
{
	//, HDLocData_TypeDef** da_LocPara
	int WholeImageWidth = ROISize*(ROISize + 1);

	// allocate memory for structure
	cudaMallocHost((void **)ha_HDLocData, sizeof(HDLocData_TypeDef)); // raw image in CPU

//	cudaMalloc((void **)da_LocPara, sizeof(HDLocData_TypeDef));  //  raw image in GPU

	HDLocData_TypeDef *h_HDLocData=*ha_HDLocData;


	// bfgs mle high density
	cudaMallocHost((void **)&h_HDLocData->h_SubRegion, PointNumTh*WholeImageWidth*sizeof(short));

	cudaMallocHost((void **)&h_HDLocData->h_LocArry, MaxPointNum*OutParaNum3*sizeof(float)); // maximum nubmer of possible point number
	cudaMallocHost((void **)&h_HDLocData->h_PointNum, sizeof(int)); // maximum nubmer of possible point number
	
	cudaMalloc((void **)&h_HDLocData->d_SubRegion, MaxPointNum*WholeImageWidth*sizeof(short));

	cudaMalloc((void **)&h_HDLocData->d_LocArry1, MaxPointNum*OutParaNum1*sizeof(float)); // maximum number of single point region
	cudaMalloc((void **)&h_HDLocData->d_LocArry2, MaxPointNum*OutParaNum2*sizeof(float)); // maximum number of two    point region
	cudaMalloc((void **)&h_HDLocData->d_LocArry3, MaxPointNum*OutParaNum3*sizeof(float)); // maximum number of three  point region

	cudaMalloc((void **)&h_HDLocData->d_Point1Num, sizeof(int));
	cudaMalloc((void **)&h_HDLocData->d_Point2Num, sizeof(int));
	cudaMalloc((void **)&h_HDLocData->d_Point3Num, sizeof(int));

	cudaMalloc((void **)&h_HDLocData->d_Point1PosArry, MaxPointNum*sizeof(int));
	cudaMalloc((void **)&h_HDLocData->d_Point2PosArry, MaxPointNum*sizeof(int));
	cudaMalloc((void **)&h_HDLocData->d_Point3PosArry, MaxPointNum*sizeof(int));


}

void HDLoc_Deinit(HDLocData_TypeDef* h_HDLocData)
{
	//, HDLocData_TypeDef* d_HDLocPara
	cudaFreeHost(h_HDLocData->h_SubRegion); // maximum nubmer of possible point number
	cudaFreeHost(h_HDLocData->h_LocArry); // maximum nubmer of possible point number
	cudaFreeHost(h_HDLocData->h_PointNum); // maximum nubmer of possible point number

	cudaFree(h_HDLocData->d_SubRegion);

	cudaFree(h_HDLocData->d_LocArry1); // maximum number of single point region
	cudaFree(h_HDLocData->d_LocArry2); // maximum number of two    point region
	cudaFree(h_HDLocData->d_LocArry3); // maximum number of three  point region

	cudaFree(h_HDLocData->d_Point1Num);
	cudaFree(h_HDLocData->d_Point2Num);
	cudaFree(h_HDLocData->d_Point3Num);

	cudaFree(h_HDLocData->d_Point1PosArry);
	cudaFree(h_HDLocData->d_Point2PosArry);
	cudaFree(h_HDLocData->d_Point3PosArry);

	// free
	cudaFreeHost(h_HDLocData);
//	cudaFree(d_HDLocPara);
}



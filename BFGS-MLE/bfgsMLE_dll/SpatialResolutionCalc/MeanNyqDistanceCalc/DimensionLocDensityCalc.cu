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


#include "DimensionLocDensityCalc.h"




float NyqDimensionDensityCalc_TypeDef::GetMeanMinNeighborDistance_AFrameGap(float PixelSize_nm, int OffsetFrame, int CurFrame, bool Is3DImaging, cudaStream_t cstream)

{
	// from CalcStartFrame:CurFrame, and negnect frame before OffsetFrame
	// CalcStartFrame must not less than OffsetFrame
	//	in a group, after consecutive filter, the first several frames are not filtered completely, so they should be rejected

	int CalcStartPos = GetFrameStartPos(OffsetFrame);
	int CurAccuFluoNum = GetAccumulatedFluoNum(OffsetFrame, CurFrame);

	float DistanceTh_pixel;
	float MeanDistance;

	// compute mean value of some molecule's distance compared to all others
	// 10000 is enough, but after consecutive, valid fluo number will significantly decrease, thus select more molecules
	int SelDatLen = min(NEIGHBOR_DISTANCE_CALC_DATA_SELECT_NUMBER, CurAccuFluoNum);

	CalcMinNeighborDistance(d_LocArry, PixelSize_nm, CalcStartPos, SelDatLen, Is3DImaging, CurAccuFluoNum, cstream);


	// normal mean distance calculation
	DistanceTh_pixel = 200;
	MeanDistance = CalcMeanDistance(DistanceTh_pixel, SelDatLen, cstream);;


	// filter background noise
	DistanceTh_pixel = MeanDistance * 5.0f;
	MeanDistance = CalcMeanDistance(DistanceTh_pixel, SelDatLen, cstream);


	return MeanDistance;
}


void NyqDimensionDensityCalc_TypeDef::GetDimensionLocDensity_AGroup(int OffsetFrame, int StartFrame, int EndFrame, float PixelSize_nm, bool Is3DImaging, cudaStream_t cstream)
{
	OffsetFrame = max(OffsetFrame, 2);
	OffsetFrame = min(OffsetFrame, 6);

	StartFrame = max(StartFrame, OffsetFrame);
	StartFrame = min(StartFrame, ImagesPerGroup_Valid / 2);

	EndFrame = min(EndFrame, ImagesPerGroup_Valid);

	ValidAccuFrameVaryDataNum = EndFrame - StartFrame + 1;


	AccumulatedFrameNumVary.clear();
	AccumulatedMeanDistanceVary.clear();


	for (int fcnt = 0; fcnt < ValidAccuFrameVaryDataNum; fcnt++)
	{

		int CurFrame = StartFrame + fcnt;

		float MeanDistance_Pixel = GetMeanMinNeighborDistance_AFrameGap(PixelSize_nm, OffsetFrame, CurFrame, Is3DImaging, cstream);

		// valid frame number and corresponding distance
		AccumulatedFrameNumVary.push_back(CurFrame - OffsetFrame + 1);
		AccumulatedMeanDistanceVary.push_back(MeanDistance_Pixel * PixelSize_nm / 1000.0f);

	}

	GetDimensionDensity();

}


void NyqDimensionDensityCalc_TypeDef::GetDimensionDensity()
{
	bool ExecptionOccor = false;

	float *ix = new float[ValidAccuFrameVaryDataNum];
	float *iy = new float[ValidAccuFrameVaryDataNum];

	for (int cnt = 0; cnt < ValidAccuFrameVaryDataNum; cnt++)
	{
		if (AccumulatedMeanDistanceVary[cnt] <= 0.0000001f)
		{
			ExecptionOccor = true;
			return;
		}

		ix[cnt] = logf(AccumulatedFrameNumVary[cnt]);
		iy[cnt] = logf(AccumulatedMeanDistanceVary[cnt]);
	}


	// para type, ParaNum, ItNum, ItNum_bs; 
	BFGSOptimizer< float, 2, 5, 11> LineOrder1CurveFit(DimFit_LineOrder1_PreFitting, DimFit_LineOrder1_TargerF);

	LineOrder1CurveFit.BFGSOptimize(ix, iy, ValidAccuFrameVaryDataNum);


	// y=a*x+b
	float a = LineOrder1CurveFit.FitPara[0];
	float b = LineOrder1CurveFit.FitPara[1];
	
	if ((a < -2) || (a > -0.25))
	{
		ExecptionOccor = true;
		return;

	}

	float n = -1 / a;
	float d = exp((logf(0.5f) - b)*n);

	DimensionFD = n;
	LocDensityFD = d;

	//	printf("dimension: %f %f\n", n, d);


	delete[] ix;
	delete[] iy;
}


int NyqDimensionDensityCalc_TypeDef::AddLocArray_FewFrames(float * h_iLocArry, int FluoNum, float ConsecFilter_DistanceTh_pixel, bool IsEndCalc, cudaStream_t cstream)
{

	memcpy(&h_LocArry[AccumulatedFluoNum*OutParaNumGS2D], h_iLocArry, FluoNum*OutParaNumGS2D*sizeof(float));

	AccumulatedFluoNum += FluoNum;

	int firstFrame = GetFirstFrame(h_LocArry, AccumulatedFluoNum);
	int lastFrame = GetLastFrame(h_LocArry, AccumulatedFluoNum);
	int totalFrame = lastFrame - firstFrame + 1;

	// 2500 is commonly enough, while in real localization, there are some invalid molecules, such as background noise, consecutive emitting
	const int LeastFluoNum = 4000; // 2500/0.7

	int IsEnough = 0;

	// enough molecules and frames
	if ((AccumulatedFluoNum >= LeastFluoNum) && (totalFrame >= ImagesPerGroup_Ideal))
	{
		IsEnough = 1;
	}

	if ((AccumulatedFluoNum >= 40000) && (totalFrame >= 35))
	{
		IsEnough = 1;
	}

	// too many molecules
	if (AccumulatedFluoNum >= (MAX_FLUO_NUM_PER_GROUP - 1.1*PointNumTh))
	{
		IsEnough = 1;
	}

	// too many frames
	if (totalFrame >= MAX_FRAME_NUM_PER_GROUP*0.8f)
	{
		IsEnough = 1;
	}

	// at the end of calculation, if the frame is roughly
	if (IsEndCalc && (totalFrame >= 30))
	{
		IsEnough = 1;
	}

	if (IsEnough)
	{

		AddLocArray_AGroup(h_LocArry, AccumulatedFluoNum, ConsecFilter_DistanceTh_pixel, cstream);
	}

	return IsEnough;
}


void NyqDimensionDensityCalc_TypeDef::AddLocArray_AGroup(float * h_iLocArry, int FluoNum, float ConsecFilter_DistanceTh_pixel, cudaStream_t cstream)
{

	FindAllFramePos(h_iLocArry, FluoNum);

	cudaMemcpyAsync(d_LocArry, h_iLocArry, FluoNum*OutParaNumGS2D*sizeof(float), cudaMemcpyHostToDevice, cstream);

	// filter consecutive molecules
	ConsecutiveFilter.FilterConsecutiveFluo(d_LocArry, FluoNum, ConsecFilter_DistanceTh_pixel, cstream);


	cudaStreamSynchronize(cstream);

}


void NyqDimensionDensityCalc_TypeDef::ResetAccumulatedData()
{
	AccumulatedFluoNum = 0;

	ImagesPerGroup_Valid = 0;
}

void NyqDimensionDensityCalc_TypeDef::Init(int MaxFluoNumPerGroup, int ImagesPerGroup)
{
	cudaError_t err;

	this->ImagesPerGroup_Ideal = ImagesPerGroup;

	err = cudaMallocHost((void **)&h_LocArry, MaxFluoNumPerGroup*OutParaNumGS2D*sizeof(float));
	HandleErr(err, "cudaMallocHost MeanNyqDistanceCalc h_LocArry");

	err = cudaMalloc((void **)&d_LocArry, MaxFluoNumPerGroup*OutParaNumGS2D*sizeof(float));
	HandleErr(err, "cudaMalloc MeanNyqDistanceCalc d_LocArry");

	// mean and min distance
	cudaMallocHost((void **)&h_ValidNum, sizeof(float));
	cudaMallocHost((void **)&h_TotalValue, sizeof(float));

	cudaMalloc((void **)&d_ValidNum, sizeof(float));
	cudaMalloc((void **)&d_TotalValue, sizeof(float));

	cudaMalloc((void **)&d_MinDistance, NEIGHBOR_DISTANCE_CALC_DATA_SELECT_NUMBER * 2 * sizeof(float));

	//
	ConsecutiveFilter.Init(MaxFluoNumPerGroup);

	ResetAccumulatedData();

	DimensionFD = 0;
	LocDensityFD = 0;

}



void NyqDimensionDensityCalc_TypeDef::DeInit()
{
	cudaFreeHost(h_LocArry);
	cudaFree(d_LocArry);


	cudaFreeHost(h_ValidNum);
	cudaFreeHost(h_TotalValue);

	cudaFree(d_ValidNum);
	cudaFree(d_TotalValue);

	cudaFree(d_MinDistance);


	ConsecutiveFilter.DeInit();

}

void NyqDimensionDensityCalc_TypeDef::CalcMinNeighborDistance(float *d_LocArry, float PixelSize_nm, int CalcStartPos, int SelDatLen, bool Is3DImaging, int FluoNum, cudaStream_t cstream)
{

	int BlockDim = ThreadsPerBlock;

	int BlockNum = (SelDatLen + ThreadsPerBlock - 1) / ThreadsPerBlock;


	cudaMemsetAsync(d_MinDistance, 0, SelDatLen * sizeof(float), cstream);

	if (!Is3DImaging)
	{
		gpuMinDistanceCalc_2D << <BlockNum, BlockDim, 0, cstream >> >(d_LocArry, d_MinDistance, PixelSize_nm, CalcStartPos, SelDatLen, FluoNum);
	}
	else
	{
		gpuMinDistanceCalc_3D << <BlockNum, BlockDim, 0, cstream >> >(d_LocArry, d_MinDistance, PixelSize_nm, CalcStartPos, SelDatLen, FluoNum);
	}

	cudaStreamSynchronize(cstream);

}



float NyqDimensionDensityCalc_TypeDef::CalcMeanDistance(float DistanceTh_pixel, int SelDatLen, cudaStream_t cstream)
{

	int BlockDim = ThreadsPerBlock;

	int BlockNum = (SelDatLen + ThreadsPerBlock - 1) / ThreadsPerBlock;


	cudaMemsetAsync(d_ValidNum, 0, sizeof(float), cstream);
	cudaMemsetAsync(d_TotalValue, 0, sizeof(float), cstream);

	gpuMeanDistanceCalc << <BlockNum, BlockDim, 0, cstream >> >(d_MinDistance, d_ValidNum, d_TotalValue, DistanceTh_pixel, SelDatLen);


	//	cudaMemcpyAsync(h_MinDistance, d_MinDistance, TotalFluoNum*sizeof(float), cudaMemcpyDeviceToHost, cstream);
	cudaMemcpyAsync(h_ValidNum, d_ValidNum, sizeof(float), cudaMemcpyDeviceToHost, cstream);
	cudaMemcpyAsync(h_TotalValue, d_TotalValue, sizeof(float), cudaMemcpyDeviceToHost, cstream);

	cudaStreamSynchronize(cstream);

	float MeanDistance = (*h_TotalValue) / (*h_ValidNum);

	return MeanDistance;

}


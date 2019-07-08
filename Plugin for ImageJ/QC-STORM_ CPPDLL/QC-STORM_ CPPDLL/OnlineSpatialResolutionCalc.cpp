#include "stdafx.h"
#include "OnlineSpatialResolutionCalc.h"

#include <time.h>

tbb::concurrent_queue<qLocArray> LocArray_Resolution_Queue;


UINT th_OnlineSpatialResolutionCalc(LPVOID params)
{
	cudaSetDevice(GPUID_2Best);


	int CurDevice = 0;

	cudaGetDevice(&CurDevice);
	printf("Resolution dev: %d\n", CurDevice);


	LocArray_Resolution_Queue.clear();

	// spatial resolution calculation
	DimensionDensityCalc.ResetAccumulatedData();
	DimensionDensityCalc.FrameNumberPerGroupCalc(LocPara_Global.ImageWidth, LocPara_Global.ImageHigh);



	SpatialResolutionCalc.ResetData();
	SpatialResolutionCalc.SetStructureSize(LocPara_Global.StrucuteSize_2D);
	

	qLocArray LocArray_Rec;

	
	bool IsBreak = false;

	unsigned int ResolutionTime = 0;
	unsigned int time1 = 0;

	float *WriteLocArry = NULL;
	int WriteLocNum = 0;


	while (1)
	{
		// get localizations send from OnlineLocalizationLD
		if (LocArray_Resolution_Queue.try_pop(LocArray_Rec))
		{
			time1 = clock();

			WriteLocArry = LocArray_Rec.h_LocArray;
			WriteLocNum = LocArray_Rec.FluoNum;
			IsBreak = LocArray_Rec.IsEndCalc;

			// dimension and density calculation, only calculate after enough frames are accumulated

			int IsEnough = DimensionDensityCalc.AddLocArray_FewFrames(WriteLocArry, WriteLocNum, LocPara_Global.ConsecFit_DistanceTh_nm / LocPara_Global.PixelSize, IsBreak, Resolution_stream1);
		
			delete[] WriteLocArry;

			if (IsEnough)
			{
				//	bool Is3DImaging = LocPara_Global.LocType == LocType_AS3D;
				bool Is3DImaging = false;

				DimensionDensityCalc.GetDimensionLocDensity_AGroup(DimensionDensityCalc.ImagesPerGroup_Valid * 1 / 10, DimensionDensityCalc.ImagesPerGroup_Valid * 2 / 10, DimensionDensityCalc.ImagesPerGroup_Valid, LocPara_Global.PixelSize, Is3DImaging, Resolution_stream1);

				//					printf("Dim data: %.2f %.2f %d\n", DimensionDensityCalc.DimensionFD, DimensionDensityCalc.LocDensityFD, DimensionDensityCalc.ImagesPerGroup_Valid);

				// add current group's dimension and density to calculate spatial resolution
				SpatialResolutionCalc.AddDimensionLocDensity_AGroup(DimensionDensityCalc.DimensionFD, DimensionDensityCalc.LocDensityFD, DimensionDensityCalc.ImagesPerGroup_Valid);


				// get spatial resolution vary data

				float Mean_LocPrecisionXY = FluoStatData.TimeVaryMean_LocPrecisionXY;;
				float Mean_LocPrecisionZ = Mean_LocPrecisionXY / 1.414f * 2.0f;

				SpatialResolutionCalc.GetSpatialResolutionVary(Is3DImaging, LocPara_Global.IsHollowTube, Mean_LocPrecisionXY, Mean_LocPrecisionZ, NYQUIST_RESOLUTION_OVERSAMPLING);


				DimensionDensityCalc.ResetAccumulatedData();
			}


			ResolutionTime += (clock() - time1);

		}

		if (IsBreak)break;
	}


	printf("Resolution calc time : %d ms\n", ResolutionTime);


	cudaGetDevice(&CurDevice);
	printf("Resolution dev: %d\n", CurDevice);

	return 0;
}



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

#include "stdafx.h"
#include "SpatialResolutionCalc.h"

#include <math.h>

#define max(a,b)    (((a) > (b)) ? (a) : (b))
#define min(a,b)    (((a) < (b)) ? (a) : (b))


template<typename DataType>
DataType GetAccumulatedVectorValue(vector<DataType>&iData_Vary)
{
	DataType TotalValue = 0;

	for (int i = 0; i < iData_Vary.size(); i++)
	{
		TotalValue += iData_Vary[i];
	}

	return TotalValue;
}

template<typename DataType>
void PrintData_core(vector<DataType> &iDataVary, const char*pstr)
{
	printf("%s :\n", pstr);

	for (int cnt = 0; cnt < iDataVary.size(); cnt++)
	{
		printf("%.2f ", (float)iDataVary[cnt]);
	}
	printf("\n");

}


void SpatialResolutionCalc_TypeDef::AddDimensionLocDensity_AGroup(float DimensionFD, float LocDensityFD, int FrameNum_CurGroup)
{
	//
	FrameNum_Vary_Group.push_back(FrameNum_CurGroup);

	AccumulatedFrameNum_Vary_Group.push_back(GetAccumulatedValue(FrameNum_Vary_Group));
	
	if ((DimensionFD > 0) && (DimensionFD < 3.5))
	{
		// avoid nan
	}
	else
	{
		DimensionFD = 5;
		LocDensityFD = 0;
	}

	// process abnormal value
	if (DimensionFD > 3.5f)
	{
		int VarySize = Dimension_Vary_Group_FD.size();
		if (VarySize >= 2)
		{
			float *DimensionFDVaryData = Dimension_Vary_Group_FD.data();
			float *LocDensityFDVaryData = LocDensity_Vary_Group_FD.data();

			DimensionFD = (DimensionFDVaryData[VarySize - 2] + DimensionFDVaryData[VarySize - 1]) / 2;
			LocDensityFD = (LocDensityFDVaryData[VarySize - 2] + LocDensityFDVaryData[VarySize - 1]) / 2;
		}
	}

	//
	Dimension_Vary_Group_FD.push_back(DimensionFD);
	LocDensity_Vary_Group_FD.push_back(LocDensityFD);


	// convert localization density in FD to 2d 3d
	LocDensity_Vary_Group_2D.push_back(expf(2 / DimensionFD*logf(LocDensityFD)));
	// convert localization density in FD to 2d 3d
	LocDensity_Vary_Group_3D.push_back(expf(3 / DimensionFD*logf(LocDensityFD)));


	// accumulated localization density
	AccumulatedLocDensity_Vary_FD.push_back(GetAccumulatedDensity(LocDensity_Vary_Group_FD, FrameNum_Vary_Group));
	AccumulatedLocDensity_Vary_2D.push_back(GetAccumulatedDensity(LocDensity_Vary_Group_2D, FrameNum_Vary_Group));
	AccumulatedLocDensity_Vary_3D.push_back(GetAccumulatedDensity(LocDensity_Vary_Group_3D, FrameNum_Vary_Group));

}


void SpatialResolutionCalc_TypeDef::GetSpatialResolutionVary(bool Is3DImaging, bool IsHollowTube, float MeanLocPrecXY, float MeanLocPrecZ, float OverSamplingRatio)
{

	// MeanLocPrecXY is sqrt(loc prec X^2 + loc prec Y^2)
	

	float DimensionConversionSize_2D = 0;
	float DimensionConversionSize_3D = 0;

	GetDimensionConversionWidth(&DimensionConversionSize_2D, &DimensionConversionSize_3D, Is3DImaging, IsHollowTube, MeanLocPrecXY, MeanLocPrecZ);


	// each data point is seperated by 10 frames
	NyqResolutionVary_FD_10f.clear(); // mean Nyquist resolution
	NyqResolutionVary_2D_10f.clear(); // mean Nyquist resolution
	NyqResolutionVary_3D_10f.clear(); // mean Nyquist resolution

	NyquistResolutionVary_10f.clear(); // mean Nyquist resolution
	SpatialResolutionVary_10f.clear();


	// valid display data point number
	int AccumulatedFrameNum = GetAccumulatedFrameNum();
	int TotalGroupNum_LocDensity = GetLocDensityGroupNum();

	int ResolutionVaryNum_10f = AccumulatedFrameNum / SPATIAL_RESOLUTION_FRAME_GAP;


	if (TotalGroupNum_LocDensity < 1)return;


	for (int cnt = 0; cnt < ResolutionVaryNum_10f; cnt++)
	{
		int CurFrame = (cnt + 1)*SPATIAL_RESOLUTION_FRAME_GAP;
		
		int WholeGroupId = GetGroupId_CurFrame(CurFrame);
		int CurGroupId = min(WholeGroupId + 1, TotalGroupNum_LocDensity - 1);
		
		int ResiFrame = 0;

		float AccuDensityOffset_FD = 0;
		float AccuDensityOffset_2D = 0;
		float AccuDensityOffset_3D = 0;

		if (WholeGroupId < 0)
		{
			AccuDensityOffset_FD = 0;
			AccuDensityOffset_2D = 0;
			AccuDensityOffset_3D = 0;
			ResiFrame = CurFrame;
		}
		else
		{
			AccuDensityOffset_FD = AccumulatedLocDensity_Vary_FD[WholeGroupId];
			AccuDensityOffset_2D = AccumulatedLocDensity_Vary_2D[WholeGroupId];
			AccuDensityOffset_3D = AccumulatedLocDensity_Vary_3D[WholeGroupId];
			ResiFrame = CurFrame - AccumulatedFrameNum_Vary_Group[WholeGroupId];
		}


		float MeanDimensionFD = GetMeanDimensionFD(CurGroupId);

		float CurAccuDensity_FD = AccuDensityOffset_FD + LocDensity_Vary_Group_FD[CurGroupId] * ResiFrame;
		float CurAccuDensity_2D = AccuDensityOffset_2D + LocDensity_Vary_Group_2D[CurGroupId] * ResiFrame;
		float CurAccuDensity_3D = AccuDensityOffset_3D + LocDensity_Vary_Group_3D[CurGroupId] * ResiFrame;


		// nm
		NyqResolutionVary_FD_10f.push_back(2 * 1000 * 0.5*powf(CurAccuDensity_FD, -1.0f / MeanDimensionFD));
		NyqResolutionVary_2D_10f.push_back(2 * 1000 * 0.5*powf(CurAccuDensity_2D, -1.0f / 2));
		NyqResolutionVary_3D_10f.push_back(2 * 1000 * 0.5*powf(CurAccuDensity_3D, -1.0f / 3));


		// 2D conversion size is first met
		if ((NyqResolutionVary_FD_10f[cnt] <= DimensionConversionSize_2D) && (MetStructureSize2D_Pos < 0))
		{
			MetStructureSize2D_Pos = cnt;
			MetStructureSize2D_Value = NyqResolutionVary_FD_10f[cnt];
		}

			
		float CorrectedNyquistResolution = NyqResolutionVary_FD_10f[cnt];

		// apply 2D correction by normalization
		if (NyqResolutionVary_FD_10f[cnt] <= DimensionConversionSize_2D)
		{
			if ((!Is3DImaging) || (Is3DImaging && (MeanDimensionFD < 2)))
			{
				CorrectedNyquistResolution = MetStructureSize2D_Value * NyqResolutionVary_2D_10f[cnt] / NyqResolutionVary_2D_10f[MetStructureSize2D_Pos];
			}
		}

		// 3D conversion size is first met
		if ((CorrectedNyquistResolution <= DimensionConversionSize_3D) && (MetStructureSize3D_Pos < 0))
		{
			MetStructureSize3D_Pos = cnt;
			MetStructureSize3D_Value = CorrectedNyquistResolution;
		}

		// apply 3D correction by normalization
		if (CorrectedNyquistResolution <= DimensionConversionSize_3D)
		{
			if (Is3DImaging)
			{
				CorrectedNyquistResolution = MetStructureSize3D_Value * NyqResolutionVary_3D_10f[cnt] / NyqResolutionVary_3D_10f[MetStructureSize3D_Pos];
			}
		}

		float LocPrec_FWHM = MeanLocPrecXY*2.35f;
		float CurValidNyqResolution = sqrtf(OverSamplingRatio)*CorrectedNyquistResolution;
		float CurSpatialResolution = sqrtf(LocPrec_FWHM*LocPrec_FWHM + CurValidNyqResolution*CurValidNyqResolution);

		NyquistResolutionVary_10f.push_back(CurValidNyqResolution);

		SpatialResolutionVary_10f.push_back(CurSpatialResolution);

	}
}


float SpatialResolutionCalc_TypeDef::GetCurSpatialResolution()
{
	int SpatialResolutionDataNum = SpatialResolutionVary_10f.size();

	float CurResolution = 1000;

	if (SpatialResolutionDataNum > 0)
	{
		CurResolution = SpatialResolutionVary_10f[SpatialResolutionDataNum-1];
	}

	return CurResolution;
}

float SpatialResolutionCalc_TypeDef::GetCurNyquistResolution()
{

	int NyqResolutionDataNum = NyquistResolutionVary_10f.size();

	float CurResolution = 1000;

	if (NyqResolutionDataNum > 0)
	{
		CurResolution = NyquistResolutionVary_10f[NyqResolutionDataNum - 1];
	}

	return CurResolution;
}


int SpatialResolutionCalc_TypeDef::GetEstimatedAcquisitionFrame(float RSCResolution_th, float DimensionFD, float LocDensityFD, bool Is3DImaging, bool IsHollowTube, float MeanLocPrecXY, float MeanLocPrecZ, float OverSamplingRatio)
{

#define MAX_ACQ_FRAME_NUMBER			50000
#define ACQ_ESTIMATION_GROUP_FRAME		50

	// MeanLocPrecXY is sqrt(loc prec X^2 + loc prec Y^2)

	if ((DimensionFD <= 0) || (DimensionFD > 3.2) || (LocDensityFD <= 0) || (RSCResolution_th < 1))
	{
		return 0;
	}

	//

	int MetStructureSize2D_Pos_Estimate = -1;
	int MetStructureSize3D_Pos_Estimate = -1;

	float MetStructureSize2D_Value_Estimate = 0;
	float MetStructureSize3D_Value_Estimate = 0;

	float DimensionConversionSize_2D = 0;
	float DimensionConversionSize_3D = 0;

	GetDimensionConversionWidth(&DimensionConversionSize_2D, &DimensionConversionSize_3D, Is3DImaging, IsHollowTube, MeanLocPrecXY, MeanLocPrecZ);


	int MaxAcqFrameNum = MAX_ACQ_FRAME_NUMBER;

	int ResolutionVaryGroupNum_Estimate = MAX_ACQ_FRAME_NUMBER / ACQ_ESTIMATION_GROUP_FRAME;

	// each data point is seperated by 100 frames
	vector<float> NyqResolutionVary_FD_Estimate; // mean Nyquist resolution
	vector<float> NyqResolutionVary_2D_Estimate; // mean Nyquist resolution
	vector<float> NyqResolutionVary_3D_Estimate; // mean Nyquist resolution


	float LocDensity2D = expf(2 / DimensionFD*logf(LocDensityFD));
	float LocDensity3D = expf(3 / DimensionFD*logf(LocDensityFD));


	for (int cnt = 0; cnt < ResolutionVaryGroupNum_Estimate; cnt++)
	{
		int CurFrame = (cnt + 1)*ACQ_ESTIMATION_GROUP_FRAME;

		float CurAccuDensity_FD = CurFrame*LocDensityFD;
		float CurAccuDensity_2D = CurFrame*LocDensity2D;
		float CurAccuDensity_3D = CurFrame*LocDensity3D;

		// nm
		NyqResolutionVary_FD_Estimate.push_back(2 * 1000 * 0.5*powf(CurAccuDensity_FD, -1.0f / DimensionFD));
		NyqResolutionVary_2D_Estimate.push_back(2 * 1000 * 0.5*powf(CurAccuDensity_2D, -1.0f / 2));
		NyqResolutionVary_3D_Estimate.push_back(2 * 1000 * 0.5*powf(CurAccuDensity_3D, -1.0f / 3));


		// 2D conversion size is first met
		if ((NyqResolutionVary_FD_Estimate[cnt] <= DimensionConversionSize_2D) && (MetStructureSize2D_Pos_Estimate < 0))
		{
			MetStructureSize2D_Pos_Estimate = cnt;
			MetStructureSize2D_Value_Estimate = NyqResolutionVary_FD_Estimate[cnt];
		}


		float CorrectedNyquistResolution = NyqResolutionVary_FD_Estimate[cnt];

		// apply 2D correction by normalization
		if (NyqResolutionVary_FD_Estimate[cnt] <= DimensionConversionSize_2D)
		{
			if ((!Is3DImaging) || (Is3DImaging && (DimensionFD < 2)))
			{
				CorrectedNyquistResolution = MetStructureSize2D_Value_Estimate * NyqResolutionVary_2D_Estimate[cnt] / NyqResolutionVary_2D_Estimate[MetStructureSize2D_Pos_Estimate];
			}
		}
		// 3D conversion size is first met
		if ((CorrectedNyquistResolution <= DimensionConversionSize_3D) && (MetStructureSize3D_Pos_Estimate < 0))
		{
			MetStructureSize3D_Pos_Estimate = cnt;
			MetStructureSize3D_Value_Estimate = CorrectedNyquistResolution;
		}

		// apply 3D correction by normalization
		if (CorrectedNyquistResolution <= DimensionConversionSize_3D)
		{
			if (Is3DImaging)
			{
				CorrectedNyquistResolution = MetStructureSize3D_Value_Estimate * NyqResolutionVary_3D_Estimate[cnt] / NyqResolutionVary_3D_Estimate[MetStructureSize3D_Pos_Estimate];
			}
		}

		float LocPrec_FWHM = MeanLocPrecXY*2.35f;
		float CurValidNyqResolution = sqrtf(OverSamplingRatio)*CorrectedNyquistResolution;
		float CurSpatialResolution = sqrtf(LocPrec_FWHM*LocPrec_FWHM + CurValidNyqResolution*CurValidNyqResolution);


		if (CurSpatialResolution <= RSCResolution_th)
		{
			MaxAcqFrameNum = CurFrame;
			break;
		}
	}

	return MaxAcqFrameNum;
}

void SpatialResolutionCalc_TypeDef::GetDimensionConversionWidth(float *oWidth2D, float *oWidth3D, bool Is3DImaging, bool IsHollowTube, float MeanLocPrecXY, float MeanLocPrecZ)
{
	float DimensionConversionSize_2D = 0;
	float DimensionConversionSize_3D = 0;

	float MeanLocPrec;

	if (!Is3DImaging)
	{
		MeanLocPrec = MeanLocPrecXY / 1.414f; // only in x or y
	}
	else
	{
		MeanLocPrec = (MeanLocPrecXY / 1.414f * 2 + MeanLocPrecZ) / 3.0f;
	}

	if (!Is3DImaging)
	{
		DimensionConversionSize_2D = StructureSize_2D + MeanLocPrec*2.35f;
	}
	else
	{
		DimensionConversionSize_2D = 3.14f * StructureSize_2D + MeanLocPrec*2.35f;
	}

	if (Is3DImaging)
	{
		// IsHollowTube : StructureSize_2D significantly larget than MeanLocPrec*2.35
		if (IsHollowTube)
		{
			DimensionConversionSize_3D = sqrt(MeanLocPrecXY*MeanLocPrecXY + MeanLocPrecZ*MeanLocPrecZ)*2.35f;
		}
		else
		{
			DimensionConversionSize_3D = StructureSize_2D + MeanLocPrec*2.35f; // microtube
		}
	}

	*oWidth2D = DimensionConversionSize_2D;
	*oWidth3D = DimensionConversionSize_3D;

//	printf("get dim covert width:%.2f %.2f %.2f %.2f\n", DimensionConversionSize_2D, DimensionConversionSize_3D, MeanLocPrecXY, MeanLocPrecZ);
}

int SpatialResolutionCalc_TypeDef::GetAccumulatedFrameNum()
{
	int AccuFrameNum = 0;
	int GroupNum = AccumulatedFrameNum_Vary_Group.size();

	if (GroupNum > 0)
	{
		AccuFrameNum = AccumulatedFrameNum_Vary_Group[GroupNum - 1];
	}

	return AccuFrameNum;
}

int SpatialResolutionCalc_TypeDef::GetGroupId_CurFrame(int CurFrame)
{

	int GroupNum = GetLocDensityGroupNum();

	int cnt = 0;

	for (cnt = 0; cnt < GroupNum; cnt++)
	{
		if (AccumulatedFrameNum_Vary_Group[cnt] > CurFrame)break;
	}

	int WholeGroupId = cnt - 1;

	return WholeGroupId;
}

int SpatialResolutionCalc_TypeDef::GetLocDensityGroupNum()
{
	return FrameNum_Vary_Group.size();

}

float SpatialResolutionCalc_TypeDef::GetAccumulatedDensity(vector<float> &iLocDensity_Vary_Group, vector<int> &iFrameNum_Vary_Group)
{
	int GroupSize = min(iLocDensity_Vary_Group.size(), iFrameNum_Vary_Group.size());

	float AccumulatedDensity = 0;

	for (int cnt = 0; cnt < GroupSize; cnt++)
	{
		AccumulatedDensity += iLocDensity_Vary_Group[cnt] * iFrameNum_Vary_Group[cnt];
	}

	return AccumulatedDensity;
}


float SpatialResolutionCalc_TypeDef::GetMeanDimensionFD(int WholeGroupId)
{
	int TotalGroupNum_LocDensity = GetLocDensityGroupNum();

	WholeGroupId = min(WholeGroupId, TotalGroupNum_LocDensity - 1);

	float TotalDat = 0;
	int ValidNum = 0;
	float CurDat;

	if (TotalGroupNum_LocDensity > 0)
	{
		for (int cnt = 0; cnt <= WholeGroupId; cnt++)
		{
			CurDat = Dimension_Vary_Group_FD[cnt];
			if ((CurDat > 0) && (CurDat <= 3))
			{
				TotalDat += CurDat;
				ValidNum += 1;
			}
		}
	}

	if (ValidNum == 0)ValidNum = 1;

	return TotalDat / ValidNum;

}

int SpatialResolutionCalc_TypeDef::GetAccumulatedValue(vector<int> & iData_Vary)
{
	return GetAccumulatedVectorValue<int>(iData_Vary);
}

float SpatialResolutionCalc_TypeDef::GetAccumulatedValue(vector<float> & iData_Vary)
{
	return GetAccumulatedVectorValue<float>(iData_Vary);

}


void SpatialResolutionCalc_TypeDef::ResetData()
{

	// store dimension and density of each group data
	FrameNum_Vary_Group.clear();
	AccumulatedFrameNum_Vary_Group.clear();

	Dimension_Vary_Group_FD.clear();
	LocDensity_Vary_Group_FD.clear();
	LocDensity_Vary_Group_2D.clear();
	LocDensity_Vary_Group_3D.clear();



	// each data point is seperated by 10 frames
	NyqResolutionVary_FD_10f.clear();
	NyqResolutionVary_2D_10f.clear();
	NyqResolutionVary_3D_10f.clear();

	// each data point is seperated by 10 frames
	NyquistResolutionVary_10f.clear();
	SpatialResolutionVary_10f.clear();


	//
	MetStructureSize2D_Pos = -1;
	MetStructureSize3D_Pos = -1;

	MetStructureSize2D_Value = 0;
	MetStructureSize3D_Value = 0;


}

void SpatialResolutionCalc_TypeDef::Init()
{
	StructureSize_2D = 40; // nm


	ResetData();
}

void SpatialResolutionCalc_TypeDef::DeInit()
{


}

void SpatialResolutionCalc_TypeDef::SetStructureSize(float StructureSize_2D)
{
	this->StructureSize_2D = StructureSize_2D;
}

void SpatialResolutionCalc_TypeDef::PrintData()
{
	PrintData_core<int>(FrameNum_Vary_Group, "FrameNum_Vary_Group");
	PrintData_core<float>(Dimension_Vary_Group_FD, "Dimension_Vary_Group_FD");
	PrintData_core<float>(LocDensity_Vary_Group_FD, "LocDensity_Vary_Group_FD");
	PrintData_core<float>(LocDensity_Vary_Group_2D, "LocDensity_Vary_Group_2D");

	PrintData_core<int>(AccumulatedFrameNum_Vary_Group, "AccumulatedFrameNum_Vary_Group");
	PrintData_core<float>(AccumulatedLocDensity_Vary_FD, "AccumulatedLocDensity_Vary_FD");
	PrintData_core<float>(AccumulatedLocDensity_Vary_2D, "AccumulatedLocDensity_Vary_2D");

	PrintData_core<float>(NyqResolutionVary_FD_10f, "NyqResolutionVary_FD_10f");
	PrintData_core<float>(NyqResolutionVary_2D_10f, "NyqResolutionVary_2D_10f");

	PrintData_core<float>(NyquistResolutionVary_10f, "NyquistResolutionVary_10f");
	PrintData_core<float>(SpatialResolutionVary_10f, "SpatialResolutionVary_10f");
}

float SpatialResolutionCalc_TypeDef::GetMeanVectorValue(vector<float> &iData_Vary_Group)
{
	float SumData = 0;
	int DataNum = 0;

	int VectorSize = iData_Vary_Group.size();
	if (VectorSize > 1)
	{
		for (int cnt = 0; cnt < VectorSize-1; cnt++)
		{
			SumData += iData_Vary_Group[cnt];
			DataNum++;
		}
	}

	if (DataNum <= 0)DataNum = 1;

	return SumData / DataNum;
}



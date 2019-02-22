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

#include "statisticsInfo.h"
#include "cudaWrapper.h"


float FluoStatisticData_TypeDef::GetActivationDensity(float Ov1MoleculesRatio, float RadiusTh_um)
{
	float density = 0;
	// some times when image is very small and molecular is very little, the ratio can be 1, lead the density to be 0
	if (Ov1MoleculesRatio <= 0.001f)Ov1MoleculesRatio = 0.001f;
	if (Ov1MoleculesRatio > 1)Ov1MoleculesRatio = 1;

	float ratio1 = logf(Ov1MoleculesRatio);
	float pi = 3.14159f;

	density = -ratio1 / pi / RadiusTh_um / RadiusTh_um;

	// a compensation from practical to theoretical
	density = density*density* 0.3156f + density* 1.186f; // second order polynomial compensation

	return density;
}


float FluoStatisticData_TypeDef::GetHistogramMeanData(int *HistData, int DatLen, float PercentTh)
{
	//PercentTh=0.05; or 0.1
	int MaxDistribVal = 0;

	float CurValue = 0;

	float WeightedSum = 0;
	float OrdinarySum = 0;
	float MeanDat = 0;
	int cnt;

	MaxDistribVal = GetHistogramMaxData(HistData, DatLen);
	PercentTh = PercentTh*MaxDistribVal;

	// get mean by center of mass

	for (cnt = 0; cnt < DatLen; cnt++)
	{
		CurValue = HistData[cnt];

		if (CurValue >= PercentTh)
		{
			WeightedSum += CurValue * (cnt + 0.5f);
			OrdinarySum += CurValue;
		}
	}

	if (OrdinarySum == 0)OrdinarySum = 1;

	MeanDat = WeightedSum / OrdinarySum;

	return MeanDat;
}

int FluoStatisticData_TypeDef::GetHistogramMaxData(int *HistData, int DatLen)
{
	int cnt;
	int MaxDat = HistData[0];

	for (cnt = 0; cnt < DatLen; cnt++)
	{
		MaxDat = max(MaxDat, HistData[cnt]);
	}

	return MaxDat;
}

int FluoStatisticData_TypeDef::GetHistogramMaxDataPos(int *HistData, int DatLen)
{
	int cnt;
	int MaxPos = 0;
	float MaxDat = 0;

	for (cnt = 0; cnt < DatLen; cnt++)
	{
		if (MaxDat < HistData[cnt])
		{
			MaxDat = HistData[cnt];
			MaxPos = cnt;
		}
	}

	return MaxPos;
}

float FluoStatisticData_TypeDef::GetHistogramWidth(int *HistData, int MaxPos, int DatLen)
{
	// get gaussian distribution width (full width at half maximum)(2.2.3548f*sigma) by center of mass

	float WeightedSum = 0;
	float OrdinarySum = 0;

	float SigmaL;
	float SigmaR;
	int cnt = 0;

	if (MaxPos < 1)MaxPos = 1;
	if (MaxPos > DatLen - 1)MaxPos = DatLen - 1;



	WeightedSum = 0;
	OrdinarySum = 0;
	for (cnt = 1; cnt < MaxPos; cnt++)
	{
		WeightedSum += HistData[cnt] * (MaxPos - cnt);
		OrdinarySum += HistData[cnt];
	}
	SigmaL = WeightedSum / OrdinarySum;

	WeightedSum = 0;
	OrdinarySum = 0;
	for (cnt = MaxPos; cnt < DatLen - 1; cnt++)
	{
		WeightedSum += HistData[cnt] * (cnt - MaxPos);
		OrdinarySum += HistData[cnt];
	}
	SigmaR = WeightedSum / OrdinarySum;

	SigmaL = (SigmaL + SigmaR) / 2;

	return SigmaL;

}



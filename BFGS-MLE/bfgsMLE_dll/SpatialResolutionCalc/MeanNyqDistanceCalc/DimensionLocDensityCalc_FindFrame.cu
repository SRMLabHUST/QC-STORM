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

#pragma once

#include <iostream>
using namespace std;

#include "DimensionLocDensityCalc.h"



void NyqDimensionDensityCalc_TypeDef::FindAllFramePos(float * h_LocArry, int TotalFluoNum)
{
	int FirstFrame = GetFirstFrame(h_LocArry, TotalFluoNum);
	int EndFrame = GetLastFrame(h_LocArry, TotalFluoNum);

//	printf("find frame:%d %d %d\n", FirstFrame, EndFrame, TotalFluoNum);

	FrameStartPosArry.clear();
	FrameFluoNumArry.clear();


	int CurFramePos = 0;
	int CurFluoNum = 0;

	int cnt = 0;
	for (cnt = FirstFrame; cnt <= EndFrame; cnt++)
	{
		CurFluoNum = GetFluoNumPerFrame(h_LocArry, TotalFluoNum, cnt, CurFramePos);


		FrameStartPosArry.push_back(CurFramePos);
		FrameFluoNumArry.push_back(CurFluoNum);

		//
		CurFramePos += CurFluoNum;
	}

	ImagesPerGroup_Valid = EndFrame - FirstFrame + 1;
}


int NyqDimensionDensityCalc_TypeDef::GetAccumulatedFluoNum(int OffsetFrame, int CurFrame)
{
	int CalcStartPos = GetFrameStartPos(OffsetFrame);
	int CurAccuFluoNum = GetFrameEndPos(CurFrame) - CalcStartPos;

	return CurAccuFluoNum;

}

int NyqDimensionDensityCalc_TypeDef::GetFrameStartPos(int Frame)
{
	int pos = FrameStartPosArry[Frame - 1];

	return pos;
}


int NyqDimensionDensityCalc_TypeDef::GetFrameEndPos(int Frame)
{
	int pos = FrameStartPosArry[Frame - 1] + FrameFluoNumArry[Frame - 1];

	return pos;
}


int NyqDimensionDensityCalc_TypeDef::GetFirstFrame(float * h_LocArry, int FluoNum)
{
	float(*pLocArry)[OutParaNumGS2D]; // for parameter array

	pLocArry = (float(*)[OutParaNumGS2D])h_LocArry;

	int cnt = 0;
	int curFrame = 0;
	for (cnt = 0; cnt < FluoNum; cnt++)
	{
		curFrame = pLocArry[cnt][Pos_Frme];

		if (curFrame != 0)break;
	}


	return curFrame;
}


int NyqDimensionDensityCalc_TypeDef::GetLastFrame(float * h_LocArry, int FluoNum)
{
	float(*pLocArry)[OutParaNumGS2D]; // for parameter array

	pLocArry = (float(*)[OutParaNumGS2D])h_LocArry;

	int cnt = 0;
	int curFrame = 0;
	for (cnt = FluoNum - 1; cnt > 0; cnt--)
	{
		curFrame = pLocArry[cnt][Pos_Frme];

		if (curFrame != 0)break;
	}


	return curFrame;
}

// find end pos and molecule number from its start position
int NyqDimensionDensityCalc_TypeDef::GetFluoNumPerFrame(float * h_LocArry, int FluoNum, int CurFrame, int StartPos)
{
	float(*pLocArry)[OutParaNumGS2D]; // for parameter array

	pLocArry = (float(*)[OutParaNumGS2D])h_LocArry;


	int FluoNum_CurFrame = 0;

	for (int cnt = StartPos; cnt < FluoNum; cnt++)
	{
		int tFrame = pLocArry[cnt][Pos_Frme];

		if (tFrame > CurFrame)
		{
			break;
		}
		else
		{
			FluoNum_CurFrame++;
		}
	}

	return FluoNum_CurFrame;
}



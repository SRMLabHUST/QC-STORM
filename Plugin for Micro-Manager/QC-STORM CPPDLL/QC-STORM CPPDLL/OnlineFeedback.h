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
#include <vector>
#include <string>

using namespace std;

// average data among several data points
#define Feedback_SamplingNum_Average				4

// wait data to be stable after feedback
#define Feedback_SamplingNum_WaitStable			1




class FeedbackCtlData
{

#define CtlTypeNum			2 // auto and manual target value
#define CtlType_Auto		0
#define CtlType_Manual		1


private:
	vector<float>TimeVaryData;

	// feedback target value
	float TargetValue[CtlTypeNum];
	int TargetValid[CtlTypeNum];

	volatile int CtlMode; // 0: auto target value, 1:nanual target value
	volatile int FeedbackEnable;

public:

	FeedbackCtlData()
	{
		FeedbackEnable = 0;
		memset(TargetValue, 0, CtlTypeNum*sizeof(float));
		memset(TargetValid, 0, CtlTypeNum*sizeof(int));

		CtlMode = CtlType_Auto; // auto targer value


		ResetAutoTarget();
	}

	void ResetAutoTarget()
	{
		TimeVaryData.clear();

		TargetValue[CtlType_Auto] = 0;
		TargetValid[CtlType_Auto] = 0;

	}

	void SetCtlMode(int Mode)
	{
		CtlMode = Mode;
		TimeVaryData.clear();
	}

	void SetFeedbackEn(int en)
	{
		FeedbackEnable = (en > 0);
	}

	int IsFeedbackEnable()
	{
		return FeedbackEnable;
	}

	int IsTargetValid()
	{
		return TargetValid[CtlMode];
	}

	float GetTargetValue()
	{
		return TargetValue[CtlMode];
	}

	void SetManualTargetValue(float iManualTarget)
	{
		SetCtlMode(CtlType_Manual);

		TargetValue[CtlType_Manual] = iManualTarget;
		TargetValid[CtlType_Manual] = 1;

		printf("set manual target:%f \n", TargetValue[CtlType_Manual]);
	}

	int GetVaryDataGroupNum()
	{
		// average data among several data points
		int GroupNum = TimeVaryData.size() / Feedback_SamplingNum_Average;

		return GroupNum;
	}

	void UpdateVaryData(float idat)
	{

		TimeVaryData.push_back(idat);

		// goal is not set, then set auto goal
		if (TargetValid[CtlType_Auto] == 0)
		{
			if (GetVaryDataGroupNum() > 0)
			{
				TargetValue[CtlType_Auto] = GetGroupMeanData_FromEnd(1);
				TargetValid[CtlType_Auto] = 1;

				printf("get auto goal:%f\n", TargetValue[CtlType_Auto]);
			}
		}
	}

	// group id start from 1
	float GetGroupMeanData_FromEnd(int GroupId)
	{
		float MeanData = 0;

		if (GetVaryDataGroupNum() >= GroupId)
		{
			float SumData = 0;
			int ValidNum = 0;

			int StartID = TimeVaryData.size() - Feedback_SamplingNum_Average;
			int EndID = TimeVaryData.size();

			for (int cnt = StartID; cnt < EndID; cnt++)
			{
				SumData += TimeVaryData[cnt];
				ValidNum++;
			}

			if (ValidNum <= 0)ValidNum = 1;
			MeanData = SumData / ValidNum;
		}

		return MeanData;
	}

};


void StartSendCmdUARTRunning();


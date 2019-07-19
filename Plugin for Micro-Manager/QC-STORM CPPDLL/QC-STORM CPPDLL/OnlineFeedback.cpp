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

#include "OnlineLocalizationLD.h"
#include "PIDController.h"
#include "OnlineFeedback_rc.h"


void UpdateFeedbackDat();


// uart cmd
#define UARTCMDBufSize				64

// 10ms wait if no cmd
#define FeedbackThreadIdleTime		5


ThreadCmdProcessState FeedbackState_LocDensity;

volatile int IsFeedbackCmdUARTRunning = 0;


// by adjusting activation laser via MCU
UINT th_OnlineFeedback_LocalizationDensity(LPVOID params)
{

	// send to MCU
	char *UARTBuf = new char[UARTCMDBufSize];
	string cmdString;

	LocDensity_CtlData.ResetAutoTarget();


	while (OnlineFeedbackAlive)
	{
		if (FeedbackState_LocDensity.HaveProcessWait())
		{
			FeedbackState_LocDensity.ProcessFinish();

			// activation laser density
			if (LocDensity_CtlData.IsFeedbackEnable() && LocDensity_CtlData.IsTargetValid())
			{

				float LocDensityTarget = LocDensity_CtlData.GetTargetValue();

				LocDensity_PID.SetSetpoint(LocDensityTarget);


				float CurLocDensity = LocDensity_CtlData.GetGroupMeanData_FromEnd(1);


				float LocDensityCtl_Ajust = LocDensity_PID.Compute(CurLocDensity);


				printf("LocDensity feedback:%.2f %%.2f %%.2f\n", LocDensityTarget, CurLocDensity, LocDensityTarget - CurLocDensity);


				// sent to the MCU to drive the step motor to adjust laser power by rotating the ND filter	
				sprintf(UARTBuf, "@a%.4f#", LocDensityCtl_Ajust);
				cmdString = UARTBuf;

				UARTCmdQueue.push(cmdString);

			}
		}
		else
		{
			Sleep(FeedbackThreadIdleTime);
		}
	}

	delete[]UARTBuf;

	return 0;
}



// receive commands from multi threads from specific hardware control
// send commands to MCU by UART
UINT th_SendFeedbackCmdToMCU(LPVOID params)
{
	string cmdString;
	UARTCmdQueue.clear();

	IsFeedbackCmdUARTRunning = 1;

	while (OnlineFeedbackAlive)
	{
		// get cmd from parallel queue
		if (UARTCmdQueue.try_pop(cmdString))
		{
			printf("send cmd:%s\n", cmdString.c_str());


			if (FeedbackCmdUART.isOpen())
			{
				// send cmd
				FeedbackCmdUART.write(cmdString);

			}
		}

	}

	UARTCmdQueue.clear();

	/*
	// feedback finish
	if (FeedbackCmdUART.isOpen())
	{
		FeedbackCmdUART.close();
	}
	*/

	IsFeedbackCmdUARTRunning = 0;

	return 0;
}



void StartSendCmdUARTRunning()
{
	if (IsFeedbackCmdUARTRunning == 0)
	{
		UARTCmdQueue.clear();

		AfxBeginThread(th_SendFeedbackCmdToMCU, NULL);

	}
}

// feedback state machine
// feed some data and wait some inteval for parameter stabalization
UINT th_OnlineFeedback(LPVOID params)
{

	// send feedback cmd from localization stastical data to feedback state machine
	OnlineFeedbackState.Reset();

	// send feedback cmd from feedback state machine to specific feedback hardware control
	FeedbackState_LocDensity.Reset();


	// threads for specific feedback hardware control
	AfxBeginThread(th_OnlineFeedback_LocalizationDensity, NULL);


	StartSendCmdUARTRunning();


	// implemented a state machine to decide how to organize feedback data
#define FeedbackState_Wait			0
#define FeedbackState_Feed			1

	unsigned int FeedbackState = FeedbackState_Wait;


	unsigned int RecGroupCnt = 0;


	// same with online processing
	while (OnlineFeedbackAlive)
	{
		if (OnlineFeedbackState.HaveProcessWait())
		{
			OnlineFeedbackState.ProcessFinish();

			// 1 feedback for 10 acquisition, and wait 2 acquisition for data stable
			switch (FeedbackState)
			{
			// wait data stable after feedback
			case FeedbackState_Wait:

				RecGroupCnt++;
				if (RecGroupCnt >= Feedback_SamplingNum_WaitStable)
				{
					RecGroupCnt = 0;
					FeedbackState = FeedbackState_Feed;
				}
				else
				{
					FeedbackState = FeedbackState_Wait;
				}

				break;

			// feedback every Feedback_SamplingNum_Average measure
			case FeedbackState_Feed:

				UpdateFeedbackDat();

				RecGroupCnt++;
				if (RecGroupCnt >= Feedback_SamplingNum_Average)
				{
					RecGroupCnt = 0;
					FeedbackState = FeedbackState_Wait;


					// feedback data is prepared
					FeedbackState_LocDensity.MakeAProcess();

				}
				else
				{
					FeedbackState = FeedbackState_Feed;
				}

				break;
			default:
				RecGroupCnt = 0;
				FeedbackState = FeedbackState_Wait;

				break;
			}
		}
		else
		{
			Sleep(FeedbackThreadIdleTime);
		}
	}
	
	return 0;
}



void UpdateFeedbackDat()
{

	// feed localization density feedback control data
	int LocDensityArraySize = FluoStatData.TimeVary_LocDensity2D.size();

	if (LocDensityArraySize > 0)
	{
		float CurLocDensity = FluoStatData.MeanLocDensity2D;

		LocDensity_CtlData.UpdateVaryData(CurLocDensity);

	}




}


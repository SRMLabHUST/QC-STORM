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
#include "OnlineFeedback_rc.h"

#include <iostream>
using namespace std;




ThreadCmdProcessState OnlineFeedbackState;


// PID controler
PIDController LocDensity_PID;


// store data variation for feedback
FeedbackCtlData LocDensity_CtlData;



// activation density and z drift
serial::Serial FeedbackCmdUART;

// uart cmd queue, send cmd from different thread to the same device
tbb::concurrent_queue<string> UARTCmdQueue;

// translation stage
serial::Serial TranslationStageUART;


// close all serial ports after library close
CloseSerialPorts closeSerialPorts;



void ActivationLaserPowerSet(float PowerPercentage)
{

	// send to MCU
	char *UARTBuf = new char[64];
	string cmdString;


	// sent to the MCU to drive the step motor to adjust laser power by rotating the ND filter	
	sprintf(UARTBuf, "@A%.2f#", PowerPercentage);
	cmdString = UARTBuf;

	UARTCmdQueue.push(cmdString);


	delete[]UARTBuf;

}




void ResetFeedbackCtl()
{
	// send to MCU
	char *UARTBuf = new char[64];
	string cmdString;

	// reset feedback MCU para to initial state
	sprintf(UARTBuf, "@r#");
	cmdString = UARTBuf;

	UARTCmdQueue.push(cmdString);
	delete[]UARTBuf;

}
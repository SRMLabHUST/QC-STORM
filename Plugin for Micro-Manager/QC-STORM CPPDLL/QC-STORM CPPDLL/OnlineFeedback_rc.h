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

#include "LocResource.h"

// for multi-thread program useage
#include "tbb.h"
//using namespace tbb;

#include "MySerial.h"

#include "PIDController.h"

#include "OnlineFeedback.h"


extern ThreadCmdProcessState OnlineFeedbackState;


extern PIDController LocDensity_PID;


extern FeedbackCtlData LocDensity_CtlData;



// activation density and z drift
extern serial::Serial FeedbackCmdUART;

extern tbb::concurrent_queue<string> UARTCmdQueue;

// translation stage
extern serial::Serial TranslationStageUART;



void ResetFeedbackCtl();

void LocDensityTest_Set();
void LocDensityTest_Reset();


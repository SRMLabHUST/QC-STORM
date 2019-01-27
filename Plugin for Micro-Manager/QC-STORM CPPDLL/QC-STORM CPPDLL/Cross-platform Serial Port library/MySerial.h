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


#include <serial/serial.h>

//using namespace serial;

#include <iostream>
#include <vector>
using namespace std;


#include "OnlineFeedback_rc.h"


int GetFirstPortId();
int IsSerialPortValid(string Port);
int IsSerialPortValid(int PortId);


class CloseSerialPorts
{
public:

	~CloseSerialPorts()
	{

		if (FeedbackCmdUART.isOpen())
		{
			FeedbackCmdUART.close();
		}

		if (TranslationStageUART.isOpen())
		{
			TranslationStageUART.close();
		}

	}
};

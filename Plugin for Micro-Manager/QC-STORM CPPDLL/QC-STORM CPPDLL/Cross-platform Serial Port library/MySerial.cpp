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
#include "MySerial.h"

int GetFirstPortId()
{
	int FirstId = 0;

	vector<serial::PortInfo> devices_found = serial::list_ports();

	vector<serial::PortInfo>::iterator iter = devices_found.begin();

	while (iter != devices_found.end())
	{
		serial::PortInfo device = *iter++;
//		printf("(%s, %s, %s %d)\n", device.port.c_str(), device.description.c_str(), device.hardware_id.c_str(), id);

		//	
		const char *PortID = device.port.c_str();

		string str1 = device.description.c_str();

		int pos = str1.find("CP2");
		if (pos >= 0)
		{
			FirstId = atoi(PortID + 3);

			break;
		}
	}
	return FirstId;
}

int IsSerialPortValid(string Port)
{
	int Valid = 0;

	vector<serial::PortInfo> devices_found = serial::list_ports();

	vector<serial::PortInfo>::iterator iter = devices_found.begin();

	while (iter != devices_found.end())
	{
		serial::PortInfo device = *iter++;

		if (device.port == Port)
		{
			Valid = 1;
			break;
		}
	}
	return Valid;
}

int IsSerialPortValid(int PortId)
{
	char buf[10];
	sprintf(buf, "COM%d", PortId);


	return IsSerialPortValid(buf);

}

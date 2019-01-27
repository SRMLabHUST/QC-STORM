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

#include <math.h>

#define AntiWindup		1


class PIDController
{

public:
	// state data
	float Setpoint; // refrence value

	float Output; // output data

	// error
	float PTerm;
	float ITerm; // ki*total error
	float DTerm;

	float lastErr;

	// PID coefficient
	float kp, ki, kd;

	float SampleTimeRatio;
	float outMin, outMax;
	float IntMin, IntMax;


public:
	PIDController()
	{
		Setpoint = 0; // refrence value
		Output = 0; // output data

		kp = 0; // proportion coefficient
		ki = 0; // integral   coefficient
		kd = 0; // derivation coefficient

		PTerm = 0;
		ITerm = 0;
		DTerm = 0;

		outMin = -100000;
		IntMin = -100000;
		outMax = 100000;
		IntMax = 100000;

		lastErr = 0;
	}

	// set PID coefficient
	void SetPIDPara(float Kp, float Ki, float Kd = 0)
	{
		SampleTimeRatio = 1;

		kp = Kp;
		ki = Ki * SampleTimeRatio;
		kd = Kd / SampleTimeRatio;

		PTerm = 0;
		ITerm = 0;
		DTerm = 0;

		lastErr = 0;

	}

	void SetSetpoint(float sp)
	{
		Setpoint = sp;

		//		Output = 0; 
		//		ITerm = 0;
		//		lastErr = 0;
	}

	float Compute(float curdat)
	{
		float curErr = Setpoint - curdat;
		float dErr = (curErr - lastErr);// d err=-d input
		lastErr = curErr;


		PTerm = kp * curErr;
		ITerm += ki*curErr;
		DTerm = kd * dErr;

#if AntiWindup
		if (ITerm> IntMax) ITerm = IntMax;
		else if (ITerm< IntMin) ITerm = IntMin;
#endif

		Output = PTerm + ITerm + DTerm;

#if AntiWindup
		if (Output > outMax) Output = outMax;
		else if (Output < outMin) Output = outMin;
#endif

		return Output;

	}

	void SetOutputLimits(float oMin, float oMax, float iMin, float iMax)
	{
		outMin = oMin;
		outMax = oMax;
		IntMin = iMin;
		IntMax = iMax;
	}

	void SetIniOutput(float offset)
	{
		Output = offset;
	}

	float GetOutput()
	{
		return Output;
	}

};



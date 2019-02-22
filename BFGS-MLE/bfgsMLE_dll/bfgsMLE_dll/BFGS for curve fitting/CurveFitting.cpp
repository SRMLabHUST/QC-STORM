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

#include "CurveFitting.h"

/*
exp fitting 
y=a*exp(b*x)
FitPara=[a,b]

*/

void ExpFit_PreFitting(float *FitPara, float *ix, float *iy, int DataNum)
{
	FitPara[0] = iy[0];

	if ((iy[0] > 0.01) && (iy[1] > 0.01))
	{
		FitPara[1] = logf(iy[1] / iy[0]);
	}
	else
	{
		FitPara[1] = 0.0f;
	}

}

float ExpFit_TargerF(float *FitPara, float *ix, float *iy, int DataNum)
{

	float a = FitPara[0];
	float b = FitPara[1];

	float y0;

	int cnt = 0;
	float SquareError = 0;

	for (cnt = 0; cnt < DataNum; cnt++)
	{
		y0 = a*expf(b*cnt);

		SquareError += powf(y0 - iy[cnt], 2);
	}

	return SquareError;

}

/*
Gaussian Fitting 1 0
Y=a*exp(-(x-x0)^2/(2*sigma^2))
FitPara=[A,x0,sigma]
*/


void GausFit10_PreFitting(float *FitPara, float *ix, float *iy, int DataNum)
{
	float mdat = iy[0];
	float mpos = 0;

	int cnt = 0;
	for (cnt = 0; cnt < DataNum; cnt++)
	{
		if (mdat < iy[cnt])
		{
			mdat = iy[cnt];
			mpos = cnt;
		}
	}
	FitPara[0] = mdat;
	FitPara[1] = mpos;
	FitPara[2] = mpos / 2.5f;

}

float GausFit10_TargerF(float *FitPara, float *ix, float *iy, int DataNum)
{
	float a = FitPara[0];
	float x0 = FitPara[1];
	float sigma = FitPara[2];
	float y0;

	int cnt = 0;
	float SquareError = 0;

	for (cnt = 0; cnt < DataNum; cnt++)
	{
		y0 = a*expf(-(cnt - x0)*(cnt - x0) / (2 * sigma*sigma));

		SquareError += powf(y0 - iy[cnt], 2);
	}

	return SquareError;


}


/*
Gaussian Fitting 1 1
Y=a*exp(-(x-x0)^2/(2*sigma^2))+b
FitPara=[A,x0,sigma,b]
*/


void GausFit11_PreFitting(float *FitPara, float *ix, float *iy, int DataNum)
{
	float mdat = iy[0];
	float mpos = 0;

	int cnt = 0;

	for (cnt = 0; cnt < DataNum; cnt++)
	{
		if (mdat < iy[cnt])
		{
			mdat = iy[cnt];
			mpos = cnt;
		}
	}

	float b = Min(iy[0], iy[DataNum - 1]);

	FitPara[0] = mdat - b;
	FitPara[1] = mpos;
	FitPara[2] = mpos / 2.5f;
	FitPara[3] = b;

}

float GausFit11_TargerF(float *FitPara, float *ix, float *iy, int DataNum)
{
	float a = FitPara[0];
	float x0 = FitPara[1];
	float sigma = FitPara[2];
	float b = FitPara[3];

	float y0;

	int cnt = 0;
	float SquareError = 0;

	for (cnt = 0; cnt < DataNum; cnt++)
	{
		y0 = a*expf(-(cnt - x0)*(cnt - x0) / (2 * sigma*sigma)) + b;

		SquareError += powf(y0 - iy[cnt], 2);
	}

	return SquareError;


}


////
/*
Gaussian Fitting 1 0 int
Y=a*exp(-(x-x0)^2/(2*sigma^2))
FitPara=[A,x0,sigma]
*/


void GausFit10i_PreFitting(float *FitPara, int *ix, int *iy, int DataNum)
{
	int mdat = iy[0];
	int mpos = 0;

	int cnt = 0;
	for (cnt = 0; cnt < DataNum; cnt++)
	{
		if (mdat < iy[cnt])
		{
			mdat = iy[cnt];
			mpos = cnt;
		}
	}
	FitPara[0] = mdat;
	FitPara[1] = mpos;
	FitPara[2] = mpos / 2.5f;

}

float GausFit10i_TargerF(float *FitPara, int *ix, int *iy, int DataNum)
{
	float a = FitPara[0];
	float x0 = FitPara[1];
	float sigma = FitPara[2];
	float y0;

	int cnt = 0;
	float SquareError = 0;

	for (cnt = 0; cnt < DataNum; cnt++)
	{
		y0 = a*expf(-(cnt - x0)*(cnt - x0) / (2 * sigma*sigma));

		SquareError += powf(y0 - iy[cnt], 2);
	}

	return SquareError;


}

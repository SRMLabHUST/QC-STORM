/*
@ author: luchang li, luchangli1993@163.com
@ huazhong university of science and technology
@ 2017/12/07
@ free use for research
*/
#include "stdafx.h"
#include "CurveFitting.h"


void ZPlaneFit_PreFitting(float *FitPara, float *DataArrayList[3], int DataNum)
{
	float *ix = DataArrayList[0];
	float *iy = DataArrayList[1];
	float *iz = DataArrayList[2];


	FitPara[0] = 1; // a
	FitPara[1] = 1;
	FitPara[2] = -100;

	printf("prefit:%f %f %f\n", FitPara[0], FitPara[1], FitPara[2]);
}

float ZPlaneFit_TargerF(float *FitPara, float *DataArrayList[3], int DataNum)
{
	float *ix = DataArrayList[0];
	float *iy = DataArrayList[1];
	float *iz = DataArrayList[2];


	float kx = FitPara[0];
	float ky = FitPara[1];
	float b0 = FitPara[2];


	int cnt = 0;
	float SquareError = 0;

	for (cnt = 0; cnt < DataNum; cnt++)
	{
		float zdiff = iz[cnt] - (kx * ix[cnt] + ky * iy[cnt] + b0);

		SquareError += powf(zdiff, 2);
	}

	return SquareError;


}

/*
Line order 1 fitting
y=a*x+b

*/

void LineOrder1_PreFitting(float *FitPara, float *DataArrayList[2], int DataNum)
{
	float *ix = DataArrayList[0];
	float *iy = DataArrayList[1];


	FitPara[0] = (iy[DataNum - 1] - iy[0]) / (ix[DataNum - 1] - ix[0]); // a
	FitPara[1] = iy[0] - ix[0] * FitPara[0];

	printf("prefit:%f %f\n", FitPara[0], FitPara[1]);
}

float LineOrder1_TargerF(float *FitPara, float *DataArrayList[2], int DataNum)
{
	float *ix = DataArrayList[0];
	float *iy = DataArrayList[1];


	float a = FitPara[0];
	float b = FitPara[1];

	float y0;

	int cnt = 0;
	float SquareError = 0;

	for (cnt = 0; cnt < DataNum; cnt++)
	{
		y0 = a*ix[cnt] + b;

		SquareError += powf(y0 - iy[cnt], 2);
	}

	return SquareError;


}

/*
int main()
{
#define DataLen	10
float x[DataLen] = { 0,1,2,3,4,5,6,7,8,9 };
float y[DataLen] = { 0.601981941401637,1.26297128454014,2.65407909847678,3.68921450314001,4.74815159282371,5.45054159850250,6.08382137799693,7.22897696871682,8.91333736150167,9.15237801896922 };

#define DataArrayNum		2

float *DataArray[DataArrayNum] = { x,y };


BFGSOptimizer_TypeDef< float, DataArrayNum, 2, 5, 11>BFGSOptimizer(LineOrder1_PreFitting, LineOrder1_TargerF);

BFGSOptimizer.BFGSOptimize(DataArray, DataLen);

BFGSOptimizer.PrintfFitPara();

return 0;
}


*/

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

#include "cuda_runtime.h"
#include "device_launch_parameters.h"

#include <stdio.h>



template <int InputDataLen, int OutputDataLen>
__host__ __device__ void InterplotData(float *x_i, float *x_o, float *y_o, float *An, float *Bn, float *Cn, float *Dn, int(*FindPositionID)(float *x_i, float xn, int InputDataLen, float UserPara), float UserPara);


template <int InputDataLen>
__host__ __device__ void ConstractSplineLinearEquations(float *A, float *b, float *x_i, float *y_i);


template <int DataLen>
__host__ __device__ void SolveSplineEquations(float *A, float *b, float *x);


// the FindPositionID is the function to find the index i where xn >= x_i[i]
// the UserPara the para to send to FindPositionID
// the input of xi must be sorted
template <int InputDataLen, int OutputDataLen>
__host__ __device__ void SplineInterpolation(float *x_i, float *y_i, float *x_o, float *y_o, int(*FindPositionID)(float *x_i, float xn, int InputDataLen, float UserPara), float UserPara)
{
	float An[InputDataLen];
	float Bn[InputDataLen];
	float Cn[InputDataLen];
	float Dn[InputDataLen];

	// initilize
	for (int i = 0; i < InputDataLen; i++)
	{
		An[i] = y_i[i];
	}

	// sove linear equations
	float A[InputDataLen*InputDataLen];
	float b[InputDataLen];

	ConstractSplineLinearEquations<InputDataLen>(A, b, x_i, y_i);

	// solve Cn
	SolveSplineEquations<InputDataLen>(A, b, Cn);

	// solve Bn and Dn
	for (int i = 0; i<InputDataLen - 1; i++)
	{
		float hi = x_i[i + 1] - x_i[i];
		Bn[i] = 1 / hi*(An[i + 1] - An[i]) - hi / 3 * (2 * Cn[i] + Cn[i + 1]);
		Dn[i] = (Cn[i + 1] - Cn[i]) / (3 * hi);
	}

	// yn(i) = An(ParaSet) + Bn(ParaSet)*(curx-xi) + Cn(ParaSet)*(curx-xi)^2 + Dn(ParaSet)*(curx-xi)^3;
	InterplotData <InputDataLen, OutputDataLen>(x_i, x_o, y_o, An, Bn, Cn, Dn, FindPositionID, UserPara);

}




template <int InputDataLen, int OutputDataLen>
__host__ __device__ void InterplotData(float *x_i, float *x_o, float *y_o, float *An, float *Bn, float *Cn, float *Dn, int(*FindPositionID)(float *x_i, float xn, int InputDataLen, float UserPara), float UserPara)
{
	int ParaSet = 0;

	for (int i = 0; i < OutputDataLen; i++)
	{
		float xn = x_o[i];

		ParaSet = FindPositionID(x_i, xn, InputDataLen, UserPara);
		if (ParaSet < 0)ParaSet = 0;
		if (ParaSet >= InputDataLen - 1)ParaSet = InputDataLen - 2;


//		printf("ParaSet:%d\n", ParaSet);
		float xi = x_i[ParaSet];

//		printf("xn xi:%f %f\n", xn, xi);

		y_o[i] = An[ParaSet] + Bn[ParaSet] * (xn - xi) + Cn[ParaSet] * (xn - xi)*(xn - xi) + Dn[ParaSet] * (xn - xi)*(xn - xi)*(xn - xi);

	}
}

template <int InputDataLen>
__host__ __device__ void ConstractSplineLinearEquations(float *A, float *b, float *x_i, float *y_i)
{
	float(*pA)[InputDataLen] = (float(*)[InputDataLen])A; // for parameter array

														  // 
	for (int r = 0; r<InputDataLen; r++)
	{
		for (int c = 0; c<InputDataLen; c++)
		{
			pA[r][c] = 0;
		}
	}
	// construct A

	pA[0][0] = 1;
	pA[InputDataLen - 1][InputDataLen - 1] = 1;


	for (int i = 1; i < InputDataLen - 1; i++)
	{
		float h0 = x_i[i] - x_i[i - 1];
		float h1 = x_i[i + 1] - x_i[i];
		pA[i][i - 1] = h0;
		pA[i][i] = 2 * (h0 + h1);
		pA[i][i + 1] = h1;
	}

	// construct b
	b[0] = 0;
	b[InputDataLen - 1] = 0;

	for (int i = 1; i<InputDataLen - 1; i++)
	{
		float h0 = x_i[i] - x_i[i - 1];
		float h1 = x_i[i + 1] - x_i[i];

		float a0 = y_i[i - 1];
		float a1 = y_i[i];
		float a2 = y_i[i + 1];

		b[i] = 3 / h1*(a2 - a1) - 3 / h0*(a1 - a0);
	}
}




template <int DataLen>
__host__ __device__ void SolveSplineEquations(float *A, float *b, float *x)
{
	float y[DataLen];
	float Beta[DataLen];

	float(*pA)[DataLen] = (float(*)[DataLen])A; // for parameter array


	for (int i = 0; i<DataLen - 1; i++)
	{
		float ci = pA[i][i + 1];
		float bi = pA[i][i];

		if (i == 0)
		{
			Beta[i] = ci / bi;
		}
		else
		{
			float ai = pA[i][i - 1];
			Beta[i] = ci / (bi - ai*Beta[i - 1]);
		}
	}

	for (int i = 0; i<DataLen; i++)
	{
		float bi = pA[i][i];

		if (i == 0)
		{
			y[i] = b[i] / bi;
		}
		else
		{
			float ai = pA[i][i - 1];
			y[i] = (b[i] - ai*y[i - 1]) / (bi - ai*Beta[i - 1]);
		}
	}

	for (int i = DataLen - 1; i >= 0; i--)
	{
		if (i == DataLen - 1)
		{
			x[i] = y[i];
		}
		else
		{
			x[i] = y[i] - Beta[i] * x[i + 1];
		}
	}
}





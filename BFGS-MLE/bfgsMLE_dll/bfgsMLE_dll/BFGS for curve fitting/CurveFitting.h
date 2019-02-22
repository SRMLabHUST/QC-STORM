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



#include "BFGSOptimizer.h"

// initial guess function for optimization/ curve fitting
// void xx_PreFitting(float *FitPara, float *ix, float *iy, int DataNum);
// optimization target function
// float xx_TargerF(float *FitPara, float *ix, float *iy, int DataNum);




/*
exp fitting
y=a*exp(b*x)
FitPara=[a,b]

*/

void ExpFit_PreFitting(float *FitPara, float *ix, float *iy, int DataNum);
float ExpFit_TargerF(float *FitPara, float *ix, float *iy, int DataNum);

/*
Gaussian Fitting 1 0
Y=a*exp(-(x-x0)^2/(2*sigma^2))
FitPara=[A,x0,sigma]
*/

void GausFit10_PreFitting(float *FitPara, float *ix, float *iy, int DataNum);
float GausFit10_TargerF(float *FitPara, float *ix, float *iy, int DataNum);

/*
Gaussian Fitting 1 1
Y=a*exp(-(x-x0)^2/(2*sigma^2))+b
FitPara=[A,x0,sigma,b]
*/
void GausFit11_PreFitting(float *FitPara, float *ix, float *iy, int DataNum);
float GausFit11_TargerF(float *FitPara, float *ix, float *iy, int DataNum);


////
/*
Gaussian Fitting 1 0 int
Y=a*exp(-(x-x0)^2/(2*sigma^2))
FitPara=[A,x0,sigma]
*/


void GausFit10i_PreFitting(float *FitPara, int *ix, int *iy, int DataNum);
float GausFit10i_TargerF(float *FitPara, int *ix, int *iy, int DataNum);

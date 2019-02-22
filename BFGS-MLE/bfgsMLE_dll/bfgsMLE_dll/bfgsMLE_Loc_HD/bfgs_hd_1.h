#include "cuda_runtime.h"
#include "device_launch_parameters.h"

#include "device_functions.h"
#include "math_functions.h"

#include <stdio.h>

#define ThreadsPerBlock		32 //Threads Per Block


#define SystemPSF			1.3f
#define ParaSigma			(1.0f/(2.0f*SystemPSF*SystemPSF))



__global__ void MLEROILocTop_1p(unsigned short *d_SubRegion, float *d_LocArry,int *d_Point1PosArry, float Offset, float kadc, int FluoNum, int isFPGAProc);

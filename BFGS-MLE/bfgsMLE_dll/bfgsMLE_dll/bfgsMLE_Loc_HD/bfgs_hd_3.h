#include "cuda_runtime.h"
#include "device_launch_parameters.h"

#include "device_functions.h"
#include "math_functions.h"

#include <stdio.h>

#define ThreadsPerBlock		32 //Threads Per Block



__global__ void MLEROILocTop_3p(unsigned short *d_SubRegion, float *d_LocArry,int *d_Point3PosArry, float Offset, float kadc, float PsfWidth, int FluoNum, int isFPGAProc);

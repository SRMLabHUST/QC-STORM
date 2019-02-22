#include "cuda_runtime.h"
#include "device_launch_parameters.h"

#include <stdio.h>

#define ThreadsPerBlock		32 //Threads Per Block

__global__ void GPUSubregionSeperate(unsigned short *d_SubRegion, int FluoNum, int *d_Point1Num, int *d_Point1PosArry, int *d_Point2Num, int *d_Point2PosArry, int *d_Point3Num, int *d_Point3PosArry);


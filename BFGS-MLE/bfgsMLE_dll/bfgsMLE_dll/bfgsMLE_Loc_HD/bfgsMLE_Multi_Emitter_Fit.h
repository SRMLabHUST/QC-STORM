#include "cuda_runtime.h"
#include "device_launch_parameters.h"

#include "bfgs_base.h"


void HDLoc_BFGS_MLELocalization_2D_2Emitter(float * d_LocArry, unsigned short * d_ImageROI, float *d_WLEPara, int MultiFitFluoNum, int * d_MultiFitFluoPos, LocalizationPara& LocPara, cudaStream_t cstream);

void HDLoc_BFGS_MLELocalization_2D_3Emitter(float * d_LocArry, unsigned short * d_ImageROI, float *d_WLEPara, int MultiFitFluoNum, int * d_MultiFitFluoPos, LocalizationPara& LocPara, cudaStream_t cstream);

void HDLoc_BFGS_MLELocalization_AS3D_2Emitter(float * d_LocArry, unsigned short * d_ImageROI, float *d_WLEPara, int MultiFitFluoNum, int * d_MultiFitFluoPos, LocalizationPara& LocPara, cudaStream_t cstream);

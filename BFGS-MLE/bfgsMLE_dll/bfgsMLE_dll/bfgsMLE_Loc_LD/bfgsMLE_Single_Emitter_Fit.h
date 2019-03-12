#include "cuda_runtime.h"
#include "device_launch_parameters.h"

#include "bfgs_base.h"




// bfgs 2D and as3D loc function
void LDLoc_BFGS_MLELocalizationGS2D(float * d_LocArry, unsigned short * d_ImageROI, float *d_WLEPara, int * d_MultiFitFluoNum, int * d_MultiFitFluoPos, LocalizationPara& LocPara, int FluoNum, cudaStream_t cstream);

void LDLoc_BFGS_MLELocalizationAS3D(float * d_LocArry, unsigned short * d_ImageROI, float *d_WLEPara, int * d_MultiFitFluoNum, int * d_MultiFitFluoPos, LocalizationPara& LocPara, int FluoNum, cudaStream_t cstream);

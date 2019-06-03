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

#include "cuda_runtime.h"
#include "device_launch_parameters.h"


#include <stdio.h>

#include "bfgs_CommonPara.h"



// wrapper for cuda function
BFGS_MLE_API char AllocHostMemory(void **ptr, long size);
BFGS_MLE_API char FreeHostMemory(void * ptr);
BFGS_MLE_API char AllocGPUMemory(void **ptr, long size);
BFGS_MLE_API char FreeGPUMemory(void * ptr);
BFGS_MLE_API char WaitGPUStream(cudaStream_t cstream);

BFGS_MLE_API char CreatStream(cudaStream_t*pstream);
BFGS_MLE_API char CreatStreamWithPriority(cudaStream_t *pstream, int prio);

BFGS_MLE_API char FreeStream(cudaStream_t cstream);

BFGS_MLE_API void HandleErr(cudaError_t err, const char * str );
BFGS_MLE_API int CudaDeviceNum();
BFGS_MLE_API char CudaSetDevice(int id);



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

#include "cudaWrapper.h"


// cuda wrapper

char AllocHostMemory(void **ptr, long size)
{
	cudaError_t err;
	err = cudaMallocHost(ptr, size);

#if bfgs_debug
	HandleErr(err, "malloc host memory");
#endif

	if (err == cudaSuccess)return 0;
	else return 1;
}

char FreeHostMemory(void * ptr)
{
	cudaError_t err;
	err = cudaFreeHost(ptr);
#if bfgs_debug
	HandleErr(err, "free host memory");
#endif

	if (err == cudaSuccess)return 0;
	else return 1;

}

char AllocGPUMemory(void **ptr, long size)
{
	cudaError_t err;
	err = cudaMalloc(ptr, size);

#if bfgs_debug
	HandleErr(err, "malloc gpu memory");
#endif

	if (err == cudaSuccess)return 0;
	else return 1;
}

char FreeGPUMemory(void * ptr)
{
	cudaError_t err;
	err = cudaFree(ptr);
#if bfgs_debug
	HandleErr(err, "free gpu memory");
#endif

	if (err == cudaSuccess)return 0;
	else return 1;

}

char CreatStream(cudaStream_t *pstream)
{
	cudaError_t err;
	err = cudaStreamCreate(pstream);

#if bfgs_debug
	HandleErr(err, "create stream");
#endif

	if (err == cudaSuccess)return 0;
	else return 1;
}

char CreatStreamWithPriority(cudaStream_t *pstream, int prio)
{
	cudaError_t err;
	err = cudaStreamCreateWithPriority(pstream, cudaStreamNonBlocking, prio);

#if bfgs_debug
	HandleErr(err, "create stream priority");
#endif

	if (err == cudaSuccess)return 0;
	else return 1;
}


char FreeStream(cudaStream_t cstream)
{
	cudaError_t err;
	err = cudaStreamDestroy(cstream);

#if bfgs_debug
	HandleErr(err, "free stream");
#endif

	if (err == cudaSuccess)return 0;
	else return 1;
}

// wait current stream to finish the task
char WaitGPUStream(cudaStream_t cstream)
{
	cudaError_t err;

	err = cudaStreamSynchronize(cstream);

#if bfgs_debug
	HandleErr(err, "wait stream");
#endif

	if (err == cudaSuccess)return 0;
	else return 1;

}

void HandleErr(cudaError_t err, const char * str)
{
	if (err != cudaSuccess)
	{
		printf("cuda err:%s, %s\n", str, cudaGetErrorString(err));
	}
	else
	{
		printf("cuda suc:%s\n", str);
	}
}

int CudaDeviceNum()
{
	int count;
	cudaGetDeviceCount(&count);
	printf("cuda device count:%d\n", count);


	if (count<1 || count>200)
	{
		printf("error: no cuda device\n");
		return 0;
	}

	int cnt;
	for (cnt = 0; cnt<count; cnt++)
	{
		cudaDeviceProp prop;

		cudaGetDeviceProperties(&prop, cnt);

		printf("device name:%s\n", prop.name);

	}
	return count;
}

char CudaSetDevice(int id)
{
	cudaError_t err;

	err = cudaSetDevice(id);;

#if bfgs_debug
	HandleErr(err, "set device num");
#endif

	if (err == cudaSuccess)return 0;
	else return 1;
}

#pragma once

#include "OnlineLocalizationLD.h"




struct qLocArray
{
	float *h_LocArray;
	int FluoNum;

	bool IsEndCalc;
};


extern tbb::concurrent_queue<qLocArray> LocArray_Resolution_Queue;



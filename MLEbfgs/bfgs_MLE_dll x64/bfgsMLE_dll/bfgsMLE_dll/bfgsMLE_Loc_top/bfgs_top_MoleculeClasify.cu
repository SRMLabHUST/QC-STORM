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

#include "bfgs_top.h"

#include "bfgsMLE_core.h"


#include "LDROIExtraction.h"




__global__ void gpuMoleculeFitClasify(unsigned short * d_ImageROI, int ROISize, float *d_WLEPara, int * d_SingleFitFluoNum, int * d_SingleFitFluoPos, int * d_MultiFitFluoNum_2E, int * d_MultiFitFluoPos_2E, int MultiEmitterFitEn, int FluoNum)
{
	int gid = threadIdx.x + blockDim.x*blockIdx.x;

	float(*pWLEPara)[WLE_ParaNumber] = (float(*)[WLE_ParaNumber])d_WLEPara;

	if (gid < FluoNum)
	{

#if(WLE_ENABLE == 1)
		int MoleculeType = pWLEPara[gid][WLE_Para_FluoType];

#else

		int ROIWholeSize = ROISize*(ROISize + 1);

		int RegionAddrOffset = gid*ROIWholeSize;
		int AddrOffset = ROISize*ROISize;

		int FitType = d_ImageROI[RegionAddrOffset + AddrOffset + 4];

		int MoleculeType = MoleculeType_MLEFit;
		
		if (FitType != ROIType_Fit_Single)
		{
			MoleculeType = MoleculeType_MultiFit;
		}

		
#endif // WLE_ENABLE

		if (MultiEmitterFitEn)
		{

#if(!SPEED_TEST)

			if (MoleculeType <= MoleculeType_WLEFit)
			{
				int StorePos = atomicAdd(d_SingleFitFluoNum, 1);
				d_SingleFitFluoPos[StorePos] = gid;
			}
			else
			{
				int StorePos = atomicAdd(d_MultiFitFluoNum_2E, 1);
				d_MultiFitFluoPos_2E[StorePos] = gid;
			}

#else  // only for speed testing


			int StorePos = 0;
			StorePos = atomicAdd(d_SingleFitFluoNum, 1);
			d_SingleFitFluoPos[StorePos] = gid;

			StorePos = atomicAdd(d_MultiFitFluoNum_2E, 1);
			d_MultiFitFluoPos_2E[StorePos] = gid;

#endif // SPEED_TEST

		}
		else
		{
			// all molecules are single molecule fitting
			d_SingleFitFluoPos[gid] = gid;
		}
	}
}



void LDLocData_TypeDef::MoleculePreFitClasify(int ROISize, int MultiEmitterFitEn, int FluoNum, cudaStream_t cstream)
{

	if (MultiEmitterFitEn)
	{
		int BlockDim = ThreadsPerBlock;
		int BlockNum = ((FluoNum + ThreadsPerBlock - 1) / ThreadsPerBlock);

		gpuMoleculeFitClasify << <BlockNum, BlockDim, 0, cstream >> > (d_ImageROI, ROISize, d_WLEPara, h_FitPosInf->d_SingleFitFluoNum, h_FitPosInf->d_SingleFitFluoPos, h_FitPosInf->d_MultiFitFluoNum_2E, h_FitPosInf->d_MultiFitFluoPos_2E, MultiEmitterFitEn, FluoNum);
		
		cudaMemcpyAsync(h_FitPosInf->h_SingleFitFluoNum, h_FitPosInf->d_SingleFitFluoNum, sizeof(int), cudaMemcpyDeviceToHost, cstream);
		
		cudaMemcpyAsync(h_FitPosInf->h_MultiFitFluoNum_2E, h_FitPosInf->d_MultiFitFluoNum_2E, sizeof(int), cudaMemcpyDeviceToHost, cstream);

		cudaStreamSynchronize(cstream);

	}
	else
	{
		*(h_FitPosInf->h_SingleFitFluoNum) = FluoNum;
	}


}



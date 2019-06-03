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

#ifndef __STACK_RENDER_3D_TOP_H
#define __STACK_RENDER_3D_TOP_H


#include "cuda_runtime.h"
#include "device_launch_parameters.h"


#include <stdio.h>

#include "bfgs_CommonPara.h"

#include "bfgs_base.h"

#include "ImageRender_base.h"



/*
RGBImageEncodeMode 0 : peak photon as weight for each molecule, rendered by localization precision calculated by CRLB
RGBImageEncodeMode 1 : 1 as weight for each molecule, rendered by localization precision calculated by CRLB
RGBImageEncodeMode 2 : 1 as weight for each molecule, rendered by fixed localization precision
*/



class ImageRender3DStackData_TypeDef
{
public:

	float *h_LocArry;
	float *d_LocArry;

	float *h_SRIntensityImg;
	float *d_SRIntensityImg;


public:


	void Init(int SRImageWidth, int SRImageHigh, int TotalFluoNum);
	void Deinit();

	// put all molecules data once for each z depth
	// FixedlocPrec_z is double of FixedlocPrec_x
	void GetSaveImgTop(float *ih_LocArry, LocalizationPara & LocPara, float RenderZDepth, int RGBImageEncodeMode, float FixedlocPrec_x, int TotalFluoNum, cudaStream_t cstream);

};




#endif //__STACK_RENDER_3D_TOP_H


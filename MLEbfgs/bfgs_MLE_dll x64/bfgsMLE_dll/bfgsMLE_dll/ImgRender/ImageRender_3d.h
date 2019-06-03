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

//#ifndef __LOCLIZATION_RENDER3D_H
//#define __LOCLIZATION_RENDER3D_H

#include "cuda_runtime.h"
#include "device_launch_parameters.h"


#include "bfgs_base.h"
#include "ImageRender_2d.h"
#include "ImageRender_top.h"




#define ImageRender_Debug	0



// image render 3d methods

void ImgRender_FluoRender_3D(float *h_LocArry, ImageRenderData_TypeDef *h_RendData, LocalizationPara & LocPara, int RGBImageEncodeMode, float RendlocPrec, int FluoNum, cudaStream_t cstream);


void ImgRender_GetDispImg_3D(ImageRenderData_TypeDef *h_RendData, float MinZDepth, float MaxZDepth, int ColorMode_3D, float BrightRatio, int oImgWidth, int oImgHigh, int cposX, int cposY, float DispZoom, float PixelZoom, int SRImageWidth, int SRImageHigh, cudaStream_t cstream);

void ImgRender_GetSaveImg_3D(ImageRenderData_TypeDef *h_RendData, LocalizationPara & LocPara, float BrightRatio, int RGBImageEncodeMode, cudaStream_t cstream);

void ImgRender_ResetFillImg_3D(ImageRenderData_TypeDef *h_RendData, int SRImageWidth, int SRImageHigh);


float GetSRImageHistgram(float *d_SRIntensityImg, float *d_SRImageHist, float *h_SRImageHist, float MaxHistVal, int SRImageWidth, int SRImageHigh, cudaStream_t cstream);


//#endif // __LOCLIZATION_RENDER3D_H

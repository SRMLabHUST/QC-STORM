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

//#ifndef __LOCLIZATION_RENDER2D_H
//#define __LOCLIZATION_RENDER2D_H


#include "cuda_runtime.h"
#include "device_launch_parameters.h"

#include "bfgs_base.h"
#include "ImageRender_top.h"


// image render 2d methods
// some are shared with 3d rendering

void ImgRender_FluoRenderWithPrec_2D(float *h_LocArry, ImageRenderData_TypeDef *h_RendData, LocalizationPara & LocPara, int RGBImageEncodeMode, float RendlocPrec, int FluoNum, cudaStream_t cstream);


void ImgRender_GetDispImg_2D(ImageRenderData_TypeDef *h_RendData, float BrightRatio, int oImgWidth, int oImgHigh, int cposX, int cposY, float DispZoom, float PixelZoom, int SRImageWidth, int SRImageHigh, cudaStream_t cstream);
void ImgRender_GetSaveImg_2D(ImageRenderData_TypeDef *h_RendData, float BrightRatio, int RGBImageEncodeMode, int SRImageWidth, int SRImageHigh, cudaStream_t cstream);
void ImgRender_ResetFillImg_2D(ImageRenderData_TypeDef *h_RendData, int SRImageWidth, int SRImageHigh);


// appendix methods

void ImgRender_GetFillMaxVal(int *h_MaxImageVal, int *d_MaxImageVal, cudaStream_t cstream);
void ImgRender_GetHistDat(int *h_HistMaxDat, int *d_HistMaxDat, cudaStream_t cstream);
float GetLightRatio(int *h_HistMaxDat, int *d_HistMaxDat);


//#endif // __LOCLIZATION_RENDER2D_H

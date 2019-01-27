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

#include "LocResource.h"
#include "OnlineLocalizationLD.h"

#include "SuperResolutionImageDisplay.h"

// online rend
extern float *h_RendFloatImage2D;

extern ThreadCmdProcessState RenderingState;


extern int CurDispPosX;
extern int CurDispPosY;
extern float CurDispZoom;


void RendDispImage();
void RendSaveImage();

void ConvertBMPToCImg(unsigned char*poImgData, unsigned char*piImgData, int ImageWidth, int ImageHigh);

void ConvertCImgToImageJ(unsigned char*poImgData, unsigned char*piImgData, int ImageWidth, int ImageHigh);



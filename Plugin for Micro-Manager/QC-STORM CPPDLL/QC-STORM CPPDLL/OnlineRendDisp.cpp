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

#include "stdafx.h"
#include "OnlineRendDisp.h"


int CurDispPosX = 0;
int CurDispPosY = 0;
float CurDispZoom = 1.0f;


// image for display and save in ImageJ
float *h_RendFloatImage2D = NULL;

// image zoom for CImg based library
void UpdateSRImageDispPara(int key, int clkX, int clkY, int WinWidth, int WinHigh);



ThreadCmdProcessState RenderingState;



UINT th_OnlineRendDispLD(LPVOID params)
{
	cudaSetDevice(GPUID_1Best);

	IsRendRunning = true;

	int clkbtn = 0;
	int WinWidth, WinHigh;


	ResetDisplay();

	RenderingState.Reset();

	while (OnlineRendAlive)
	{
		if (RenderingState.HaveProcessWait())
		{
			RenderingState.ProcessFinish();

			RendDispImage(); //

			UpdateStatInfDisplay();


			DisplayStatInfImageBySelection(StatInfDispSel);

		}

		CImgDisp_SRImage.wait(1);
		if (CImgDisp_SRImage.button() && (CImgDisp_SRImage.mouse_x() > 0) && (CImgDisp_SRImage.mouse_y() > 0))
		{
			WinWidth = CImgDisp_SRImage.width();
			WinHigh = CImgDisp_SRImage.height();
			// adjust display ratio
			UpdateSRImageDispPara(CImgDisp_SRImage.button(), CImgDisp_SRImage.mouse_x(), CImgDisp_SRImage.mouse_y(), WinWidth, WinHigh);
			
			RenderingState.MakeAProcess();

			CImgDisp_SRImage._button = 0;
		}
	}

	IsRendRunning = false;

	return 0;
}



void RendDispImage()
{

	RendData.GetDispImgTop(LocPara_Global, 0.3f, SR_IMAGE_DISPLAY_WIDTH, SR_IMAGE_DISPLAY_HIGH, CurDispPosX, CurDispPosY, CurDispZoom, render_stream1);
	WaitGPUStream(render_stream1);


	unsigned char* pDispImg = CImg_SRImage.data();
	ConvertBMPToCImg(pDispImg, (unsigned char*)(RendData.h_DispRendImg), SR_IMAGE_DISPLAY_WIDTH, SR_IMAGE_DISPLAY_HIGH);

}

void RendSaveImage()
{
	if (RendType_Is2D(LocPara_Global.LocType))
	{
		//2d image 
		cudaMemcpyAsync(h_RendFloatImage2D, RendData.d_SRIntensityImg, LocPara_Global.SRImageWidth* LocPara_Global.SRImageHigh*sizeof(float), cudaMemcpyDeviceToHost, render_stream1);
		WaitGPUStream(render_stream1);
	}
	else
	{
		// astigmatism 3d image
		// double-helix 3d
		RendData.GetSaveImgTop(LocPara_Global, 0.7f, RGBImage_EncodeMode_4B_BRGA, render_stream1);
		WaitGPUStream(render_stream1);
	}
}


void UpdateSRImageDispPara(int key, int clkX, int clkY, int WinWidth, int WinHigh)
{
	int PixelStep;
	int movDx = 0;
	int movDy = 0;

	PixelStep = 1.0f / CurDispZoom;


	movDx = (WinWidth / 2 - clkX) * PixelStep / LocPara_Global.PixelZoom;
	movDy = (WinHigh / 2 - clkY) * PixelStep / LocPara_Global.PixelZoom;

	//	printf("move :%d,%d,-,%d,%d,\n",curPicSizeX,curPicSizeY,movDx,movDy);

	CurDispPosX = CurDispPosX - movDx;
	CurDispPosY = CurDispPosY - movDy;

	if (key & 1)
	{
		// left click

		if (PixelStep>1)PixelStep = PixelStep - 1;
		else PixelStep = 1;
		CurDispZoom = 1.0f / PixelStep;

	}
	if (key & 2)
	{
		// left click
		PixelStep = PixelStep + 1;
		if (PixelStep>20)PixelStep = 20;
		CurDispZoom = 1.0f / PixelStep;

	}
	if (key & 4)
	{
		// middle click
		CurDispPosX = LocPara_Global.ImageWidth / 2;
		CurDispPosY = LocPara_Global.ImageHigh / 2;
		CurDispZoom = min(1.0f *WinWidth / (LocPara_Global.SRImageWidth), 1.0f *WinHigh / (LocPara_Global.SRImageHigh));
	}
}


void ConvertBMPToCImg(unsigned char*poImgData, unsigned char*piImgData, int ImageWidth, int ImageHigh)
{
	// convert it to bmp raw data for ImageJ
	// format is {B,G,R,0XFF} for each pixel

	int PixelNum = ImageWidth*ImageHigh;
	int rOffsetr = 0;
	int rOffsetg = ImageWidth*ImageHigh;
	int rOffsetb = ImageWidth*ImageHigh * 2;

	int cnt = 0;
	int cpos = 0;

	for (cnt = 0; cnt < PixelNum; cnt++)
	{
		cpos = cnt * 3;

		poImgData[rOffsetb + cnt] = piImgData[cpos + 0]; // b
		poImgData[rOffsetg + cnt] = piImgData[cpos + 1]; // g
		poImgData[rOffsetr + cnt] = piImgData[cpos + 2]; // r

	}
}




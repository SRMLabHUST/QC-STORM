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
#include "OnlineLocalizationLD.h"


ThreadCmdProcessState RenderingState;



void RendDispImage()
{
	
	if (RendType_Is2D(LocPara_Global.LocType))
	{
		// 2d image 
		cudaMemcpyAsync(h_RendFloatImage2D, RendData.d_SRIntensityImg, LocPara_Global.SRImageWidth* LocPara_Global.SRImageHigh*sizeof(float), cudaMemcpyDeviceToHost, render_stream1);
		WaitGPUStream(render_stream1);

	}
	else
	{
		int EncodeMode = RGBImage_EncodeMode_4B_BRGA;

		if(IsBatchLocRunning)EncodeMode = RGBImage_EncodeMode_3B_RGB;
		else EncodeMode = RGBImage_EncodeMode_4B_BRGA;

		// astigmatism 3d image
		RendData.GetSaveImgTop(LocPara_Global, 1.0, EncodeMode, render_stream1);
		WaitGPUStream(render_stream1);

	}
}


UINT th_OnlineRendDispLD(LPVOID params)
{
	cudaSetDevice(GPUID_1Best);

	IsRendRunning = 1;

	while (OnlineRendAlive)
	{
		if (RenderingState.HaveProcessWait())
		{
			RenderingState.ProcessFinish();

			RendDispImage(); //
		}
		else
		{
			Sleep(2);
		}
	}

	IsRendRunning = 0;

	return 0;
}

// imageJ RGB image is BGRA

void ConvertCImgToImageJ(unsigned char*poImgData, unsigned char*piImgData, int ImageWidth, int ImageHigh)
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
		cpos = cnt * 4;

		poImgData[cpos + 0] = piImgData[rOffsetb + cnt]; // b
		poImgData[cpos + 1] = piImgData[rOffsetg + cnt]; // g
		poImgData[cpos + 2] = piImgData[rOffsetr + cnt]; // r
		poImgData[cpos + 3] = 0xff; // alpha

	}
}


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

#pragma once

#include <iostream>
#include <vector>
#include <string>
#include <algorithm>

using namespace std;

#include "CImg.h"
using namespace cimg_library;

#include "CImgCurveDisplayParameters.h"


/*
DISP_TYPE
0 : curve, linear
1 : curve, logf
2 : histgram, linear
3 : histgram, logf

*/

#define DISP_TYPE_CURVE					0
#define DISP_TYPE_HIST					1

#define LOG_AXIS_OFF					0
#define LOG_AXIS_ON						1

#define FIXED_DATA_LEN_OFF				0
#define FIXED_DATA_LEN_ON				1

//
#define LABEL_DISPLAY_OFF			0
#define LABEL_DISPLAY_ON			1

// LABEL_TYPE 
#define LABEL_TYPE_FLOAT				0
#define LABEL_TYPE_INT					1


const unsigned char CImg_Red[] = { 255, 0, 0 };
const unsigned char CImg_Green[] = { 0, 255, 0 };
const unsigned char CImg_Blue[] = { 0, 0, 255 };
const unsigned char CImg_Orange[] = { 248, 147, 29 }; 
const unsigned char CImg_White[] = { 255, 255, 255 }; 
const unsigned char CImg_Black[] = { 0, 0, 0 }; 

// log y is just log, log x need re-sampling, if the input is linear

template<typename YDataType, int DISP_TYPE, bool IS_LOG_AXIS_Y, bool FIXED_DATA_LEN_EN, int FIXED_DATA_LEN, bool X_LABEL_NUMBER_EN, bool Y_LABEL_NUMBER_EN, int LABEL_TYPE, unsigned int CHANNEL_NUMBER>
class CImgCurveDisplay
{
private:

	vector<float> y_data[CHANNEL_NUMBER];


	string FigureName;

	float XLabelData[LABEL_NUMBER];
	float YLabelData[LABEL_NUMBER];

	bool StartPosValueEn;
	float StartPosValue0;
	float StartPosValue1;

	bool AutoSetYLabel;


	CImg<unsigned char> *CImg_Grid;
	CImg<unsigned char> *CImg_Axis;
	CImgDisplay *CImg_Display;



public:
	CImgCurveDisplay(string FigureName, bool StartPosValueEn, float StartPosValue0, float StartPosValue1, bool AutoSetYLabel)
	{

		this->FigureName = FigureName;

		this->StartPosValueEn = StartPosValueEn;
		this->StartPosValue0 = StartPosValue0;
		this->StartPosValue1 = StartPosValue1;

		this->AutoSetYLabel = AutoSetYLabel;

		CImg_Grid = new CImg<unsigned char> (GridImgWidth, GridImgHigh, 1, 3, 0);
		CImg_Axis = new CImg<unsigned char> (AxesImgWidth, AxesImgHigh, 1, 3, 0);

		CImg_Display = new CImgDisplay();
		
		
		memset(XLabelData, 0, LABEL_NUMBER * sizeof(float));
		memset(YLabelData, 0, LABEL_NUMBER * sizeof(float));


		if (FIXED_DATA_LEN_EN)
		{
			for (int i = 0; i < FIXED_DATA_LEN; i++)
			{
				for (int chn = 0; chn < CHANNEL_NUMBER; chn++)
				{
					y_data[chn].push_back(0);

				}
			}
		}
	}

	void ClearAllData(unsigned int Channel)
	{
		if (Channel >= CHANNEL_NUMBER)return;

		if (FIXED_DATA_LEN_EN)
		{
			for (int i = 0; i < FIXED_DATA_LEN; i++)
			{
				y_data[Channel][i] = 0;
			}
		}
		else
		{
			y_data[Channel].clear();
		}
	}

	void SetAllData(YDataType *yData, int DataLen, unsigned int Channel = 0)
	{
		if (DataLen <= 0)return;
		if (Channel >= CHANNEL_NUMBER)return;

		ClearAllData(Channel);

		if (FIXED_DATA_LEN_EN)
		{
			DataLen = min(DataLen, FIXED_DATA_LEN);
			for (int i = 0; i < DataLen; i++)
			{
				float curData = yData[i];

				if (IS_LOG_AXIS_Y)
				{
					if (curData < 1)curData = 1;

					curData = logf(curData);
				}

				y_data[Channel][i] = curData;
			}
		}
		else
		{
			y_data[Channel].clear();

			for (int i = 0; i < DataLen; i++)
			{
				float curYData = yData[i];

				if (IS_LOG_AXIS_Y)
				{
					if (curYData < 1)curYData = 1;
					curYData = logf(curYData);
				}

				y_data[Channel].push_back(curYData);

			}
		}

		InsertStartValue(Channel);
	}

	void AddAData(YDataType yData, int yId = 0, unsigned int Channel = 0)
	{
		if (Channel >= CHANNEL_NUMBER)return;

		if (IS_LOG_AXIS_Y)
		{
			if (yData < 1)yData = 1;
			yData = logf(yData);
		}

		if (FIXED_DATA_LEN_EN)
		{
			if (yId < y_data[Channel].size())
			{
				y_data[Channel][yId] = yData;
			}
		}
		else
		{
			y_data[Channel].push_back(yData);

		}

		InsertStartValue(Channel);
	}

	void CreateFigure()
	{
		if (AutoSetYLabel)
		{
			SetYLabel();
		}

		InitImage();

		PlotCurve();

	}

	void UpdateDisplay(bool isOpen)
	{

		if (isOpen)
		{
			CImg_Display->resize().display(*CImg_Axis);
			CImg_Display->set_title(FigureName.c_str());
			if (CImg_Display->is_closed())
			{
				CImg_Display->show();
			}
		}
		else
		{
			CImg_Display->close();
		}
	}

	auto GetAxisImage()
	{
		return CImg_Axis;
	}
	auto GetAxisImageData()
	{
		return CImg_Axis->data();
	}

	void SetLabelValue(float MinData, float MaxData, int LabelId)
	{
		for (int cnt = 0; cnt < LABEL_NUMBER; cnt++)
		{
			float CurData = cnt*(MaxData - MinData) / LABEL_SPACE_NUMBER + MinData;

			if (LabelId == 0)
			{
				XLabelData[cnt] = CurData;
			}
			else
			{
				if (IS_LOG_AXIS_Y)CurData = expf(CurData);

				YLabelData[cnt] = CurData;
			}
		}
	}

private:
	
	void PlotCurve()
	{

//		printf("DatLen:%d\n", DatLen);

		for (int chn = 0; chn < CHANNEL_NUMBER; chn++)
		{

			float* DispData = y_data[chn].data();
			int DatLen = y_data[chn].size();

			const unsigned char* CurCurveColor = CImg_Green;

			if (chn == 0)CurCurveColor = CImg_Green;
			if (chn == 1)CurCurveColor = CImg_Blue;
			if (chn == 2)CurCurveColor = CImg_Orange;
			if (chn == 3)CurCurveColor = CImg_White;
			if (chn == 4)CurCurveColor = CImg_Red;


			if (DatLen > 0)
			{
				CImg<float>CImg_DispDat(DatLen, 1, 1, 1);

				CImg_DispDat.assign(DispData, DatLen, 1, 1, 1);

				if (DISP_TYPE == DISP_TYPE_HIST)
				{
					CImg_Grid->draw_graph(CImg_DispDat, CurCurveColor, 1, 3, 0, 0, 0);
				}
				else
				{
					CImg_Grid->draw_graph(CImg_DispDat, CurCurveColor, 1, 1, 0, 0, 0);
				}
			}
		}

		CImg_Axis->draw_image(GridExtendX - 1, GridExtendY - 1, *CImg_Grid);
	}


	void SetYLabel()
	{
		int DatLen = y_data[0].size();

		float MinY = 0;
		float MaxY = 0;

		if (DatLen > 0)
		{
			MinY = y_data[0][0];
			MaxY = y_data[0][0];

			for (int cnt = 0; cnt < DatLen; cnt++)
			{
				MinY = min(MinY, y_data[0][cnt]);
				MaxY = max(MaxY, y_data[0][cnt]);
			}
		}

		SetLabelValue(MinY, MaxY, 1);
	}


	void InitImage()
	{
		// draw grid and data label
		ResetImage();



		int Grid_XPos[LABEL_NUMBER];
		int Grid_YPos[LABEL_NUMBER];

		CImg<int> CImg_Grid_XPos(LABEL_NUMBER, 1, 1, 1, 0);
		CImg<int> CImg_Grid_YPos(LABEL_NUMBER, 1, 1, 1, 0);

		char txtbuf[200];

		for (int cnt = 0; cnt < LABEL_NUMBER; cnt++)
		{
			Grid_XPos[cnt] = cnt*XLineGap;
			Grid_YPos[cnt] = cnt*YLineGap;

				
			if (X_LABEL_NUMBER_EN)
			{
				if (LABEL_TYPE == 0)
				{
					sprintf(txtbuf, "%.2f", (float)XLabelData[cnt]);
				}
				else
				{
					sprintf(txtbuf, "%d", (int)XLabelData[cnt]);
				}

				string str1 = txtbuf;
				int str1Len = str1.size();

				CImg_Axis->draw_text(GridExtendX + Grid_XPos[cnt] - CImgDisp_FontSize*(str1Len / 3.5f), GridImgHigh + GridExtendY + 2, txtbuf, CImg_Red, 0, 1, CImgDisp_FontSize);

			}			
			if (Y_LABEL_NUMBER_EN)
			{
				if (LABEL_TYPE == 0)
				{
					sprintf(txtbuf, "%.2f", (float)YLabelData[LABEL_NUMBER - 1 - cnt]);
				}
				else
				{
					sprintf(txtbuf, "%d", (int)YLabelData[LABEL_NUMBER - 1 - cnt]);
				}
				string str1 = txtbuf;
				int str1Len = str1.size();

				CImg_Axis->draw_text(GridExtendX - str1Len * CImgDisp_FontSize / 2.0f, GridExtendY - CImgDisp_FontSize / 2 + Grid_YPos[cnt], txtbuf, CImg_Red, 0, 1, CImgDisp_FontSize);
			}
		}
			
		CImg_Grid_XPos.assign(Grid_XPos, LABEL_NUMBER);
		CImg_Grid_YPos.assign(Grid_YPos, LABEL_NUMBER);
		CImg_Grid->draw_grid(CImg_Grid_XPos, CImg_Grid_YPos, CImg_Red);

	}

	void ResetImage()
	{
		// clear old image to avoid overlap display
		memset(CImg_Grid->data(), 0, GridImgWidth*GridImgHigh * 3);
		memset(CImg_Axis->data(), 0, AxesImgWidth*AxesImgHigh * 3);

	}


	void InsertStartValue(unsigned int Channel)
	{
		float iDat0;
		float iDat1;

		if (IS_LOG_AXIS_Y)
		{
			iDat0 = logf(max(StartPosValue0, 1.0f));
			iDat1 = logf(max(StartPosValue1, 1.0f));
		}
		else
		{
			iDat0 = StartPosValue0;
			iDat1 = StartPosValue1;
		}

		if (StartPosValueEn)
		{
			if (y_data[Channel].size() > 0)y_data[Channel][0] = iDat0;
			if (y_data[Channel].size() > 1)y_data[Channel][1] = iDat1;
		}

	}
};


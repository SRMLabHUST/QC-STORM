#include "stdafx.h"
#include "DepthMapCalc.h"


/*
only developed for beads imaging for 3D imaging plan inclination adjustment


*/

#define ZDepthTh		500

#define MaxDispDepth	200
#define MinDispDepth	-200


void DepthMapCalc_TypeDef::FitZPlane(float *h_iLocArray, int iFluoNum)
{
	if (iFluoNum < 1)return;

	iFluoNum = min(iFluoNum, PointNumTh);
	FluoNum = iFluoNum;

	memcpy(h_LocArray, h_iLocArray, FluoNum*OutParaNumGS2D*sizeof(float));


	DataFilter();

	FitZPlane_core();


	CreateDepthMap();

}

void DepthMapCalc_TypeDef::CreateDepthMap()
{
	int AddrOffset = 0;

	int r, g, b;

	for (int ycnt = 0; ycnt < DepthMap_ImageHigh; ycnt++)
	{
		for (int xcnt = 0; xcnt < DepthMap_ImageWidth; xcnt++)
		{
			float CurX = xcnt*PixelGap;
			float CurY = ycnt*PixelGap;

			float CurZDepth = kx*CurX + ky*CurY + b0;


			
			GetRGBColor(&r, &g, &b,  CurZDepth);


			AddrOffset = (ycnt*DepthMap_ImageWidth + xcnt)*3;

			DepthMap[AddrOffset + 0] = r;
			DepthMap[AddrOffset + 1] = g;
			DepthMap[AddrOffset + 2] = b;

		}
	}


	char *DepthMap;


}
void DepthMapCalc_TypeDef::GetRGBColor(int *r, int *g, int *b, float CurZDepth)
{
	int hval, sval, vval;

#define MaxVVal		100.0f

	sval = 100; // 100 saturation
	vval = MaxVVal;
	
	if (CurZDepth > MaxDispDepth)CurZDepth = MaxDispDepth;
	if (CurZDepth < MinDispDepth)CurZDepth = MinDispDepth;

	float MaxZRange = MaxDispDepth - MinDispDepth;

	// 0-240
	hval = (CurZDepth - MinDispDepth) * 240 / MaxZRange; // -130-0,0-130

	hval = 240 - hval; // HSV model, v 0 is red, 240 is blue


	if (hval < 0)hval = 0;
	if (hval > 240)hval = 240;

	HSVtoRGB(r, g, b, hval, sval, vval);

}

void DepthMapCalc_TypeDef::FitZPlane_core()
{
	#define DataArrayNum  3
	float *DataArray[DataArrayNum] = {x, y, z};

	
	BFGSOptimizer_TypeDef< float, DataArrayNum, 3, 8, 12> BFGSOptimizer(ZPlaneFit_PreFitting, ZPlaneFit_TargerF);
		
	BFGSOptimizer.BFGSOptimize(DataArray, ValidDataNum);

	BFGSOptimizer.PrintfFitPara();

	kx = BFGSOptimizer.FitPara[0];
	ky = BFGSOptimizer.FitPara[1];
	b0 = BFGSOptimizer.FitPara[2];

}

void DepthMapCalc_TypeDef::DataFilter()
{
	float(*pLocArry)[OutParaNumGS2D];
	pLocArry = (float(*)[OutParaNumGS2D])h_LocArray;

	ValidDataNum = 0;
	float CurX = 0;
	float CurY = 0;
	float CurZ = 0;
	
		
		
	for (int cnt = 0; cnt < FluoNum; cnt++)
	{
		CurX = pLocArry[cnt][Pos_XPos];
		CurY = pLocArry[cnt][Pos_YPos];
		CurZ = pLocArry[cnt][Pos_ZPos];

		if ((CurX > 0) && (abs(CurZ) < ZDepthTh))
		{
			x[ValidDataNum] = CurX;
			y[ValidDataNum] = CurY;
			z[ValidDataNum] = CurZ;

			ValidDataNum++;
		}
	}

//	printf("ValidDataNum:%d\n", ValidDataNum);

}


DepthMapCalc_TypeDef::DepthMapCalc_TypeDef(int ImageWidth, int ImageHigh)
{
	h_LocArray = new float [PointNumTh*OutParaNumGS2D];
	FluoNum = 0;

	x = new float[PointNumTh];
	y = new float[PointNumTh];
	z = new float[PointNumTh];

	ValidDataNum = 0;

	//
	DepthMap_ImageWidth = ImageWidth / PixelGap;
	DepthMap_ImageHigh = ImageHigh / PixelGap;

	DepthMap = new unsigned char[DepthMap_ImageWidth*DepthMap_ImageHigh * 3];

	memset(DepthMap, 0, DepthMap_ImageWidth*DepthMap_ImageHigh * 3);

}

DepthMapCalc_TypeDef::~DepthMapCalc_TypeDef()
{
	delete[] h_LocArray;

	delete[] x;
	delete[] y;
	delete[] z;

	//
	delete[] DepthMap;
}


void DepthMapCalc_TypeDef::HSVtoRGB(int *r, int *g, int *b, int h, int s, int v)
{
	// convert from HSV/HSB to RGB color
	// R,G,B from 0-255, H from 0-260, S, V from 0-100
	// ref http://colorizer.org/

	// The hue (H) of a color refers to which pure color it resembles
	// The saturation (S) of a color describes how white the color is
	// The value (V) of a color, also called its lightness, describes how dark the color is

	if (h < 0)h = 0;
	if (s < 0)s = 0;
	if (v < 0)v = 0;

	int i;

	float RGB_min, RGB_max;
	RGB_max = v*2.55f;
	RGB_min = RGB_max*(100 - s) / 100.0f;

	i = h / 60;
	int difs = h % 60; // factorial part of h

					   // RGB adjustment amount by hue 
	float RGB_Adj = (RGB_max - RGB_min)*difs / 60.0f;

	switch (i) {
	case 0:
		*r = RGB_max + 0.5f; // + 0.5f for round off
		*g = RGB_min + RGB_Adj + 0.5f;
		*b = RGB_min + 0.5f;
		break;
	case 1:
		*r = RGB_max - RGB_Adj + 0.5f;
		*g = RGB_max + 0.5f;
		*b = RGB_min + 0.5f;
		break;
	case 2:
		*r = RGB_min + 0.5f;
		*g = RGB_max + 0.5f;
		*b = RGB_min + RGB_Adj + 0.5f;
		break;
	case 3:
		*r = RGB_min + 0.5f;
		*g = RGB_max - RGB_Adj + 0.5f;
		*b = RGB_max + 0.5f;
		break;
	case 4:
		*r = RGB_min + RGB_Adj + 0.5f;
		*g = RGB_min + 0.5f;
		*b = RGB_max + 0.5f;
		break;
	default:		// case 5:
		*r = RGB_max + 0.5f;
		*g = RGB_min + 0.5f;
		*b = RGB_max - RGB_Adj + 0.5f;
		break;
	}

	*r = max(*r, 0);
	*g = max(*g, 0);
	*b = max(*b, 0);


	*r = min(*r, 255);
	*g = min(*g, 255);
	*b = min(*b, 255);

}

void ConvertRGBToCImg(unsigned char*poImgData, unsigned char*piImgData, int ImageWidth, int ImageHigh)
{

	int PixelNum = ImageWidth*ImageHigh;
	int rOffsetr = 0;
	int rOffsetg = ImageWidth*ImageHigh;
	int rOffsetb = ImageWidth*ImageHigh * 2;

	int cnt = 0;
	int cpos = 0;

	for (cnt = 0; cnt < PixelNum; cnt++)
	{
		cpos = cnt * 3;

		poImgData[rOffsetr + cnt] = piImgData[cpos + 0]; // r
		poImgData[rOffsetg + cnt] = piImgData[cpos + 1]; // g
		poImgData[rOffsetb + cnt] = piImgData[cpos + 2]; // b

	}
}


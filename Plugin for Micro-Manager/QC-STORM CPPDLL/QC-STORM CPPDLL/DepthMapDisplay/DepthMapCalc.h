#pragma once

#include <iostream>
#include <algorithm>
using namespace std;

/*
only developed for beads imaging for 3D imaging plan inclination adjustment


*/


#include "bfgs_MLE_para.h"

#include "CurveFitting.h"


#define PixelGap		2


class DepthMapCalc_TypeDef {

public:
	float *h_LocArray;
	int FluoNum;

	//
	float *x;
	float *y;
	float *z;

	int ValidDataNum;

	//
	unsigned char *DepthMap;
	int DepthMap_ImageWidth;
	int DepthMap_ImageHigh;


public:

	DepthMapCalc_TypeDef(int ImageWidth, int ImageHigh);
	~DepthMapCalc_TypeDef();

	void FitZPlane(float *h_iLocArray, int iFluoNum);


private:
	void DataFilter();
	void FitZPlane_core();

	void CreateDepthMap();
	
	void GetRGBColor(int *r, int *g, int *b, float CurZDepth);

	void HSVtoRGB(int *r, int *g, int *b, int h, int s, int v);

private:
	float kx;
	float ky;
	float b0;

};

void ConvertRGBToCImg(unsigned char*poImgData, unsigned char*piImgData, int ImageWidth, int ImageHigh);

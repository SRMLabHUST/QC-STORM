#pragma once

#include <iostream>
using namespace std;



class FloatImage {

public:
	unsigned int ImageWidth; 
	unsigned int ImageHigh;

	float *pImage;

public:
	FloatImage(unsigned int ImageWidth, unsigned int ImageHigh)
	{
		this->ImageWidth = ImageWidth;
		this->ImageHigh = ImageHigh;
		AllocMemory();
	}

	FloatImage(float *iImage, unsigned int ImageWidth, unsigned int ImageHigh)
	{
		this->ImageWidth = ImageWidth;
		this->ImageHigh = ImageHigh;

		AllocMemory();

		CopyImageIn(iImage);
	}

	FloatImage(FloatImage* iImage)
	{
		this->ImageWidth = iImage->ImageWidth;
		this->ImageHigh = iImage->ImageHigh;

		AllocMemory();

		CopyImageIn(iImage->pImage);
	}


	void CopyImageIn(float *iImage)
	{
		memcpy(pImage, iImage, ImageWidth*ImageHigh * sizeof(float));
	}

	float PixelValue(int x, int y)
	{
		if (x < 0)x = 0;
		if (y < 0)y = 0;
		if (x >= ImageWidth)x = ImageWidth-1;
		if (y >= ImageHigh)y = ImageHigh - 1;

		return pImage[y * ImageWidth + x];
	}

	~FloatImage()
	{
		delete[]pImage;
	}
private:
	void AllocMemory()
	{
		pImage = new float[ImageWidth*ImageHigh];

	}
};



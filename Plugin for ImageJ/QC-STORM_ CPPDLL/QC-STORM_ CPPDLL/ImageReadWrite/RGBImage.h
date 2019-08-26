#pragma once


#include <iostream>
using namespace std;


#define RGBImage_BytesPerPixel	3


class RGBImage {

public:
	unsigned int ImageWidth;
	unsigned int ImageHigh;

	
	unsigned char *pImage;

public:
	RGBImage(unsigned int ImageWidth, unsigned int ImageHigh)
	{
		this->ImageWidth = ImageWidth;
		this->ImageHigh = ImageHigh;

		AllocMemory();
	}

	RGBImage(unsigned char *iImage, unsigned int ImageWidth, unsigned int ImageHigh)
	{

		this->ImageWidth = ImageWidth;
		this->ImageHigh = ImageHigh;

		AllocMemory();

		CopyImageIn(iImage);
	}

	RGBImage(RGBImage* iImage)
	{
		this->ImageWidth = iImage->ImageWidth;
		this->ImageHigh = iImage->ImageHigh;

		AllocMemory();


		CopyImageIn(iImage->pImage);
	}


	void CopyImageIn(unsigned char *iImage)
	{
		memcpy(pImage, iImage, ImageWidth * ImageHigh * RGBImage_BytesPerPixel);
	}

	unsigned char PixelValue(int x, int y, int id)
	{
		if (x < 0)x = 0;
		if (y < 0)y = 0;
		if (x >= ImageWidth)x = ImageWidth - 1;
		if (y >= ImageHigh)y = ImageHigh - 1;

		return pImage[(y * ImageWidth + x)*RGBImage_BytesPerPixel + id];
	}

	~RGBImage()
	{
		delete[]pImage;
	}

private:

	void AllocMemory()
	{
		pImage = new unsigned char[ImageWidth*ImageHigh*RGBImage_BytesPerPixel];

	}

};

#include "stdafx.h"
#include "ReadWriteTiffImage.h"


FloatImage* ReadFloatImage(CString Imagename)
{

	char pFilePath[1024];

	WideCharToMultiByte(CP_ACP, 0, Imagename.GetBuffer(), -1, pFilePath, 1024, NULL, NULL);

	TinyTIFFReaderFile* tiffr = NULL;


	tiffr = TinyTIFFReader_open(pFilePath);

	int ImageWidth = 0;
	int ImageHigh = 0;

	if (tiffr)
	{
		ImageWidth = TinyTIFFReader_getWidth(tiffr);
		ImageHigh = TinyTIFFReader_getHeight(tiffr);
	}

	FloatImage* oImage = new FloatImage(ImageWidth, ImageHigh);


	if (tiffr)
	{
		TinyTIFFReader_getSampleData(tiffr, oImage->pImage, 0); // get image data
	}

	//	printf("%f %f %f\n", oImage->pImage[267* ImageWidth+ 187], oImage->pImage[268 * ImageWidth + 188], oImage->pImage[269 * ImageWidth + 189]);


	TinyTIFFReader_close(tiffr);

	return oImage;
}


void WriteFloatImage(CString Imagename, FloatImage *iImage)
{
	char pFilePath[1024];

	WideCharToMultiByte(CP_ACP, 0, Imagename.GetBuffer(), -1, pFilePath, 1024, NULL, NULL);

	TinyTIFFFile*tiff2 = TinyTIFFWriter_open(pFilePath, sizeof(float) * 8, iImage->ImageWidth, iImage->ImageHigh);

	TinyTIFFWriter_writeImage(tiff2, iImage->pImage);

	TinyTIFFWriter_close(tiff2);

}


RGBImage* ReadRGBImage(CString Imagename)
{

	char pFilePath[1024];

	WideCharToMultiByte(CP_ACP, 0, Imagename.GetBuffer(), -1, pFilePath, 1024, NULL, NULL);

	TinyTIFFReaderFile* tiffr = NULL;


	tiffr = TinyTIFFReader_open(pFilePath);

	int ImageWidth = 0;
	int ImageHigh = 0;

	if (tiffr)
	{
		ImageWidth = TinyTIFFReader_getWidth(tiffr);
		ImageHigh = TinyTIFFReader_getHeight(tiffr);

		int bits = TinyTIFFReader_getBitsPerSample(tiffr);
		int SamplesPerPixel = TinyTIFFReader_getSamplesPerPixel(tiffr);

		printf("image size:%d %d %d %d\n", ImageWidth, ImageHigh, bits, SamplesPerPixel);
	}

	RGBImage* oImage = new RGBImage(ImageWidth, ImageHigh);


	if (tiffr)
	{
		TinyTIFFReader_getSampleData(tiffr, oImage->pImage, 0); // get image data
	}



	TinyTIFFReader_close(tiffr);

	return oImage;
}

void WriteRGBImage(CString Imagename, RGBImage *iImage)
{
	char pFilePath[1024];

	WideCharToMultiByte(CP_ACP, 0, Imagename.GetBuffer(), -1, pFilePath, 1024, NULL, NULL);

	TinyTIFFFile*tiff2 = TinyTIFFWriter_open(pFilePath,  8, iImage->ImageWidth, iImage->ImageHigh);

	TinyTIFFWriter_writeImage(tiff2, iImage->pImage, 3);

	TinyTIFFWriter_close(tiff2);

}



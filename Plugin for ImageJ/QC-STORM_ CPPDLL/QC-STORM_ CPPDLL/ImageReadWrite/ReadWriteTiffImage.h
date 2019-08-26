#pragma once


#include <iostream>
using namespace std;
#include <afx.h>

#include "FloatImage.h"
#include "RGBImage.h"


#include "tinytiffreader.h"

#include "tinytiffwriter.h"



FloatImage* ReadFloatImage(CString Imagename);

void WriteFloatImage(CString Imagename, FloatImage *iImage);

RGBImage* ReadRGBImage(CString Imagename);

void WriteRGBImage(CString Imagename, RGBImage *iImage);


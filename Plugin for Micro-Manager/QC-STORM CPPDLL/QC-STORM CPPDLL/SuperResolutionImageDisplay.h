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

#include "OnlineLocalizationLD.h"

#include "CImg.h"
using namespace cimg_library;



#define SR_IMAGE_DISPLAY_WIDTH			880
#define SR_IMAGE_DISPLAY_HIGH			880

// super resolution image display
extern CImg<unsigned char> CImg_SRImage;

extern CImgDisplay CImgDisp_SRImage;


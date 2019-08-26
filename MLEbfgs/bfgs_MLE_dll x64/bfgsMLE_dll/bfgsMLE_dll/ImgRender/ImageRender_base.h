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





// for image render

#define ThreadsPerBlock					32 //Threads Per Block

#define ImageRender_Debug					0


#define ImgRend_MaxDatHistLen			400


/*
RenderingMode 0 : peak photon as weight for each molecule, rendered by localization precision calculated by CRLB
RenderingMode 1 : 1 as weight for each molecule, rendered by localization precision calculated by CRLB
RenderingMode 2 : 1 as weight for each molecule, rendered by fixed localization precision
*/

#define RenderingMode_FittedPhoton_CalculatedLocPrec		0
#define RenderingMode_1Photon_CalculatedLocPrec				1
#define RenderingMode_1Photon_FixedLocPrec					2
#define RenderingMode_FittedPhoton_1Pixel					3


// regular RGB image
#define RGBImage_EncodeMode_3B_BGR			0 // for bitmap image
#define RGBImage_EncodeMode_3B_RGB			1 // for tiff RGB image

// for ImageJ display,  imageJ RGB image is BGRA
#define RGBImage_EncodeMode_4B_BRGA			2
#define RGBImage_EncodeMode_4B_RGBA			3



// blue is low depth or reg is low depth
#define ImgRend_ColorMode_3D_BlueToRed		0
#define ImgRend_ColorMode_3D_RedToBlue		1


#define EncodeMode_Is_3B(x)	((x==RGBImage_EncodeMode_3B_BGR)||(x==RGBImage_EncodeMode_3B_RGB))

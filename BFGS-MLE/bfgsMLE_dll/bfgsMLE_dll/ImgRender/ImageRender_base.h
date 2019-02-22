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


// regular RGB image
#define RGBImage_EncodeMode_3Bytes			0
// for ImageJ display
#define RGBImage_EncodeMode_4Bytes			1

// blue is low depth or reg is low depth
#define ImgRend_ColorMode_3D_BlueToRed		0
#define ImgRend_ColorMode_3D_RedToBlue		1



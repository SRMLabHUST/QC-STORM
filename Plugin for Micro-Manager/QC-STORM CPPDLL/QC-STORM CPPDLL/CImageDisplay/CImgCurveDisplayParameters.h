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

#define GridImgWidth			800
#define GridImgHigh				500

#define GridExtendX				60
#define GridExtendY				30

#define AxesImgWidth			(GridImgWidth + GridExtendX + GridExtendY)
#define AxesImgHigh				(GridImgHigh + GridExtendY + GridExtendY)


#define LABEL_SPACE_NUMBER		5

#define LABEL_NUMBER			(LABEL_SPACE_NUMBER + 1)

#define XLineGap				((float)(GridImgWidth - 1) / LABEL_SPACE_NUMBER)
#define YLineGap				((float)(GridImgHigh - 1) / LABEL_SPACE_NUMBER)



#define CImgDisp_FontSize		16



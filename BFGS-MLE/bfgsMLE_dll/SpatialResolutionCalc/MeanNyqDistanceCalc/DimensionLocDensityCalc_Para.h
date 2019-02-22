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


// parameter number for each molecule and x,y,z,frame posisition

#include "bfgs_CommonPara.h"




// max fluo mum for a group with 50 frame 2048*2048 images
#define MAX_FLUO_NUM_PER_GROUP								(10240*3*50)

// max images in a group
#define MAX_FRAME_NUM_PER_GROUP								2000


// calculate only min neighboring distance of some molecules, calculation of all molecules is not necessary and time consuming
#define NEIGHBOR_DISTANCE_CALC_DATA_SELECT_NUMBER			20000


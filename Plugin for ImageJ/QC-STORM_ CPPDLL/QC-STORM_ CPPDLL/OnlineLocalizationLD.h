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

#ifndef __LM_ONLINE_LOC_H
#define __LM_ONLINE_LOC_H

#include "QC_STORM_CPPDLL.h"

// for multi-thread program useage
#include "tbb.h"
using namespace tbb;


#include <tchar.h>  

#include "LocResource.h"
#include "OnlineRendDisp.h"




void OpenConsole();


UINT th_OnlineLocalizationLD(LPVOID params);


UINT th_OnlineRendDispLD(LPVOID params);



#endif // __LM_ONLINE_LOC_H

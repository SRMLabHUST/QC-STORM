#pragma once
#include <iostream>
#include <vector>
using namespace std;

#include <afx.h>

#include "LocResource.h"


extern wstring BatchProc_FolderName;

extern wstring BatchProc_FileName_Prefix;
extern wstring BatchProc_FileName_Postfix;
extern wstring BatchProc_SavePath;
extern int BatchProc_MergeLocEn;



extern int TotalFileNum;
extern int CurFileCnt;


// thread of batch localization
UINT th_BatchLocalization(LPVOID params);


void FormatBatchLocSRImageName(CString & SRImageName);

void BatchProc_WriteSRImage2D(CString SRImageName);
void BatchProc_WriteSRImage3D(CString SRImageName);

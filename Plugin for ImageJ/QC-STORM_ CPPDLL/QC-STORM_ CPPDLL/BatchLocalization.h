#pragma once
#include <iostream>
#include <vector>
using namespace std;

#include <afx.h>

#include "LocResource.h"


extern wstring BatchProc_FolderName;
extern wstring BatchProc_Postfix;
extern wstring BatchProc_SavePath;



extern int TotalFileNum;
extern int CurFileCnt;


// thread of batch localization
UINT th_BatchLocalization(LPVOID params);



bool IsStringEndWithPostfix(wstring DirName, wstring PostFix);

void SearchFilesInDir(wstring DirName, wstring PostFix, vector<wstring> & FileNameList);


void BatchProc_WriteSRImage2D(CString SRFileName);
void BatchProc_WriteSRImage3D(CString SRFileName);

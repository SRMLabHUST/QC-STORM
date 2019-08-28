#pragma once

#include <iostream>
#include <vector>
using namespace std;

bool IsStringStartWithPrefix(wstring DirName, wstring PreFix);
bool IsStringEndWithPostfix(wstring DirName, wstring PostFix);
long GetNumberInFileName(wstring Str1);
wstring GetAbsoluteFileName(wstring Str1);

void SearchFilesInDir(wstring DirName, wstring PreFix, wstring PostFix, vector<wstring> & FileNameList);

void FilesSort(vector<wstring> & iFileNameList, vector<wstring> & oFileNameList);


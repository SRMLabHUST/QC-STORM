#include "stdafx.h"
#include "FileNameProcess.h"


bool IsStringStartWithPrefix(wstring DirName, wstring PreFix)
{

	int pos = DirName.find(PreFix);

	bool Valid = false;

	if (pos == 0)
	{
		Valid = true;
	}
	return Valid;
}

bool IsStringEndWithPostfix(wstring DirName, wstring PostFix)
{
	int str_len1 = DirName.size();
	int str_len2 = PostFix.size();

	int pos = DirName.rfind(PostFix);

	bool Valid = false;

	if ((pos >= 0) && (pos == str_len1 - str_len2))
	{
		Valid = true;
	}
	return Valid;
}

long GetNumberInFileName(wstring Str1)
{
	char NumberBuf[1024];

	long Number = 0;

	for (int cnt = 0; cnt < Str1.size(); cnt++)
	{
		char c = Str1[cnt];
		if ((c >= '0') && (c <= '9'))
		{
			c = c - '0';
			Number = Number * 10 + (int)c;
		}
	}

	return Number;
}

wstring GetAbsoluteFileName(wstring Str1)
{
	int pos = Str1.rfind(L"\\");

	return Str1.substr(pos);
}

void SearchFilesInDir(wstring DirName, wstring PreFix, wstring PostFix, vector<wstring> & FileNameList)
{
	wstring pattern(DirName);

	pattern.append(L"\\*");

	WIN32_FIND_DATA data;
	HANDLE hFind;

	vector<wstring> FolderNameList;

	if ((hFind = FindFirstFile(pattern.c_str(), &data)) != INVALID_HANDLE_VALUE) {
		do {
			wstring CurFileName = data.cFileName;

			bool valid_0 = !IsStringEndWithPostfix(CurFileName, L".");
			if (valid_0)
			{
				if (data.dwFileAttributes == FILE_ATTRIBUTE_DIRECTORY)
				{
					// is folder
					FolderNameList.push_back(DirName + L"\\" + CurFileName);
				}
				else
				{
					// is file
					bool Valid = IsStringEndWithPostfix(CurFileName, L".tif");

					if (PostFix.size() > 0)
					{
						Valid &= IsStringEndWithPostfix(CurFileName, PostFix);
					}
					if (PreFix.size() > 0)
					{
						Valid &= IsStringStartWithPrefix(CurFileName, PreFix);
					}


					if (Valid)
					{
						FileNameList.push_back(DirName + L"\\" + CurFileName); //
					}

					//					wprintf(L"search file:%s\n",  data.cFileName);
				}
			}

		} while (FindNextFile(hFind, &data) != 0);
		FindClose(hFind);
	}

	// search sub-folders
	for (int i = 0; i < FolderNameList.size(); i++)
	{
		SearchFilesInDir(FolderNameList[i], PreFix, PostFix, FileNameList);
	}
}


void FilesSort(vector<wstring> & iFileNameList, vector<wstring> & oFileNameList)
{
	oFileNameList.clear();


	while (iFileNameList.size() > 0)
	{
		long Number = GetNumberInFileName(iFileNameList[0]);
		int k = 0;

		for (int i = 0; i < iFileNameList.size(); i++)
		{
			long CurNum = GetNumberInFileName(GetAbsoluteFileName(iFileNameList[i]));
			if (CurNum < Number)
			{
				Number = CurNum;
				k = i;
			}
		}

		oFileNameList.push_back(iFileNameList[k]);
		iFileNameList.erase(iFileNameList.begin() + k);
	}
}


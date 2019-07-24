
#include <iostream>
using namespace std;

#include <afx.h>
#include "TxtConvert.h"




int main(int argc, char* argv[])
{

	int cnt = 0;
	int scnt = 0;

	CString CurFile;

	CString oFileName;

	TCHAR *tbuf = new TCHAR[512];


	int namelen = 0;
	int filetype = 0;

	cout << "converting" << endl;

	for (cnt = 1; cnt < argc; cnt++)
	{

		MultiByteToWideChar(CP_ACP, MB_COMPOSITE, argv[cnt], -1, tbuf, 512);

		CurFile = tbuf;
		oFileName = CurFile;

		oFileName.TrimRight(L".txt");
		oFileName = oFileName + L"_txt.csv";

		wprintf(L"CurFile:%s\n", CurFile);



		ConvertLocResultLD(CurFile.GetBuffer(), oFileName.GetBuffer());

	

	}   

	delete[] tbuf;
	cout << "convert end" << endl;

	system("pause");
	return 0;
}


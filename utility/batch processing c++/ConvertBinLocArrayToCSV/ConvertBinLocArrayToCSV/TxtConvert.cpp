#include <iostream>
using namespace std;

#include <afx.h>
#include "TxtConvert.h"


void ConvertLocResultLD(TCHAR* ifname, TCHAR* ofname)
{
	CString CurFile = ifname;

	CString oFileName = ofname;


	int cnt;


	CFile fraw(ifname, CFile::modeRead);

	CFile fText(oFileName, CFile::modeWrite | CFile::modeCreate);



	int FluoNum = fraw.GetLength() / sizeof(float) / OutParaNumGS2D;

	printf("%d\n", FluoNum);

	float * pRawData = new float[OutParaNumGS2D];



	char *tbuf = new char[2048];

	sprintf(tbuf, "%s\n", "peak intensity (photon),x (pixel),y (pixel),z (nm),PSFSigmaX (pixel),PSFSigmaY (pixel),Total intensity (photon),background (photon),SNR (peak to background e-),CRLBx (nm),CRLBy(nm),frame");

	fText.Write(tbuf, strlen(tbuf));

	int DispNum = 10;
	int DispGap = FluoNum / DispNum;

	for (cnt = 0; cnt<FluoNum; cnt++)
	{
		if (cnt % DispGap == 0)
		{
			printf("progress %d percentage\n", cnt / DispGap*DispNum);
		}

		fraw.Read(pRawData, OutParaNumGS2D * sizeof(float));

		sprintf(tbuf, "%10.4f,%10.4f,%10.4f,%10.4f,%10.4f,%10.4f,%10.4f,%10.4f,%10.4f,%10.4f,%10.4f,%10.4f\n", pRawData[0], pRawData[1], pRawData[2], pRawData[3], pRawData[4], pRawData[5], pRawData[6], pRawData[7], pRawData[8], pRawData[9], pRawData[10], pRawData[11]);

		fText.Write(tbuf, strlen(tbuf));

	}




	delete [] pRawData;
	delete [] tbuf;

	
	fraw.Close();
	fText.Close();



}


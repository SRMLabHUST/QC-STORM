#include "stdafx.h"

#include "BatchLocalization.h"

#include "ReadWriteTiffImage.h"

// call th_OnlineLocalizationLD etc.

wstring BatchProc_FolderName;
wstring BatchProc_Postfix;
wstring BatchProc_SavePath;


int TotalFileNum = 0;
int CurFileCnt = 0;

#include "tinytiffreader.h"
#include "tinytiffwriter.h"


// allocate resource
// feed images
// save localizations

UINT th_BatchLocalization(LPVOID params)
{
	IsBatchLocRunning = true;

	vector<wstring>  FileNameList;
	SearchFilesInDir(BatchProc_FolderName, BatchProc_Postfix, FileNameList);


	TotalFileNum = FileNameList.size();


	printf("File number:%d\n", TotalFileNum);

	if (TotalFileNum <= 0)return 0;


	char TiffImgName_Buff[1024];


	for (int cnt = 0; cnt < TotalFileNum; cnt++)
	{
		ImageName = FileNameList[cnt].c_str();



		CurFileCnt = cnt + 1;

		wprintf(L"Cur file: %d in %d, name: %s\n", CurFileCnt, TotalFileNum, ImageName);


		WideCharToMultiByte(CP_ACP, 0, ImageName.GetBuffer(), -1, TiffImgName_Buff, 1024, NULL, NULL);

		TinyTIFFReaderFile* tiffr = TinyTIFFReader_open(TiffImgName_Buff);


		if (!tiffr)
		{
			printf("image open error\n");
			TinyTIFFReader_close(tiffr);

			continue;
		}


		// update image size and frame number
		printf("Reading image info ...\n");

		int TotalFrame = TinyTIFFReader_countFrames(tiffr);

		int ImageWidth = TinyTIFFReader_getWidth(tiffr);
		int ImageHigh = TinyTIFFReader_getHeight(tiffr);
		printf("Image size:%d x %d x %d frame\n", ImageWidth, ImageHigh, TotalFrame);

		if (TotalFrame <= 1)
		{
			TinyTIFFReader_close(tiffr);
			continue;
		}

		LocPara_Global.ImageWidth = ImageWidth;
		LocPara_Global.ImageHigh = ImageHigh;
		LocPara_Global.TotalFrameNum = TotalFrame;

		LocPara_Global.UpdateSRImageSize();


		ImgDataQueue.clear();

		// start batch localization, allocate resource and start parallel threads
		SetLocDataFileName_BatchProc(); // set proper localization data save name
		StartLocalizationThread();



		// feed images, just as imagej to cpp dll
		int BatchedImgNum = (2048 * 2048 / ImageWidth / ImageHigh);
		BatchedImgNum = (BatchedImgNum / 2) * 2;

		BatchedImgNum = max(BatchedImgNum, 1);

		int GroupNum = (TotalFrame + BatchedImgNum - 1) / BatchedImgNum;;


		unsigned short *RawImgBuf = new unsigned short[BatchedImgNum*ImageWidth*ImageHigh];


		qImgData CurImgInf;

		int CurImgNum = 0;

		for (int gcnt = 0; gcnt < GroupNum; gcnt++)
		{
			if (gcnt == GroupNum - 1)CurImgNum = TotalFrame - gcnt*BatchedImgNum;
			else CurImgNum = BatchedImgNum;

			// wait localization
			unsigned int MaxBufferImageNum = 200 * 2048 * 2048 / (ImageWidth*ImageHigh);

			while (ImgDataQueue.unsafe_size() > MaxBufferImageNum)
			{
				Sleep(2);
			}

			// read images and feed to ImgDataQueue
			for (int fcnt = 0; fcnt < CurImgNum; fcnt++)
			{
				// read images
				TinyTIFFReader_getSampleData(tiffr, &RawImgBuf[fcnt*ImageWidth*ImageHigh], 0); // get image data
				TinyTIFFReader_readNext(tiffr);

			}

			CreateFeedImgMemory(CurImgInf, RawImgBuf, LocPara_Global.ImageWidth, LocPara_Global.ImageHigh, CurImgNum);

			ImgDataQueue.push(CurImgInf);
		}


		delete [] RawImgBuf;

		// finish localization
		TinyTIFFReader_close(tiffr);


		// wait localization thread to finish
		while (IsLocRunning == true)
		{
			Sleep(2);
		}

		// save super-resolution image
		CString SRFileName = ImageName;
		SRFileName.TrimRight(L".tif");


		CString MultiFitStr;
		if (LocPara_Global.MultiEmitterFitEn) MultiFitStr = "_M";
		else MultiFitStr = "_S";

		CString ConsecFitStr;
		if (LocPara_Global.ConsecFitEn) ConsecFitStr.Format(L"_Consec%.0fnm", LocPara_Global.ConsecFit_DistanceTh_nm);
		else ConsecFitStr = "";

		SRFileName.Format(L"%s_result%dD%d%s%s_rend%.2fnm.tif", SRFileName, LocPara_Global.LocType + 2, LocPara_Global.ROISize, MultiFitStr, ConsecFitStr, LocPara_Global.PixelSize / LocPara_Global.PixelZoom);
		
		// modify save path
		SRFileName.TrimLeft(BatchProc_FolderName.c_str());
		SRFileName.TrimLeft(L"\\");
		SRFileName.Replace(L"\\", L"__");

		CString SavePath = BatchProc_SavePath.c_str();
		SavePath = SavePath + L"\\";
		SRFileName = SavePath + SRFileName;

		printf("Writing super-resolution image...\n");


		if (LocPara_Global.LocType == LocType_GS2D)
		{
			BatchProc_WriteSRImage2D(SRFileName);
		}
		else
		{
			BatchProc_WriteSRImage3D(SRFileName);
		}

		// release resources
		FinishLocalizationThread();

	}


	IsBatchLocRunning = false;


	printf("Batch localization finish\n");

	return 0;
}


void BatchProc_WriteSRImage2D(CString SRFileName)
{
	float* h_RendFloatImage2D = new float[LocPara_Global.SRImageWidth*LocPara_Global.SRImageHigh];


	// create 2d image 
	cudaMemcpyAsync(h_RendFloatImage2D, RendData.d_SRIntensityImg, LocPara_Global.SRImageWidth* LocPara_Global.SRImageHigh * sizeof(float), cudaMemcpyDeviceToHost, loc_stream1);
	WaitGPUStream(loc_stream1);


#define BufLen		1024
	char SRImgFileName[BufLen];
	//		sprintf_s(SRImgFileName, "ReRendered sr image_%d-%dnm.tif", GroupCnt, (int)RenderingPixelSize);

	WideCharToMultiByte(CP_OEMCP, NULL, SRFileName.GetBuffer(), -1, SRImgFileName, BufLen, NULL, FALSE);

	TinyTIFFFile* tif = TinyTIFFWriter_open(SRImgFileName, 32, LocPara_Global.SRImageWidth, LocPara_Global.SRImageHigh); // 32 bit float image
	TinyTIFFWriter_writeImage(tif, h_RendFloatImage2D);

	TinyTIFFWriter_close(tif);


	delete[]h_RendFloatImage2D;

}

void BatchProc_WriteSRImage3D(CString SRFileName)
{
	RGBImage SRImg_RGB((unsigned char *)RendData.h_SaveRendImg, LocPara_Global.SRImageWidth, LocPara_Global.SRImageHigh);

	WriteRGBImage(SRFileName, &SRImg_RGB);
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

void SearchFilesInDir(wstring DirName, wstring PostFix, vector<wstring> & FileNameList)
{
	std::wstring pattern(DirName);

	pattern.append(L"\\*");

	WIN32_FIND_DATA data;
	HANDLE hFind;

	vector<wstring> FolderNameList;

	if ((hFind = FindFirstFile(pattern.c_str(), &data)) != INVALID_HANDLE_VALUE) {
		do {
			bool valid0 = !IsStringEndWithPostfix(data.cFileName, L".");
			if (valid0)
			{
				if (data.dwFileAttributes == FILE_ATTRIBUTE_DIRECTORY)
				{
					// is folder
					FolderNameList.push_back(DirName + L"\\" + data.cFileName);
				}
				else
				{
					// is file
					bool valid1 = IsStringEndWithPostfix(data.cFileName, PostFix);
					bool valid2 = IsStringEndWithPostfix(data.cFileName, L".tif");
					if (valid1 && valid2)
					{
						FileNameList.push_back(DirName + L"\\" + data.cFileName);
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
		SearchFilesInDir(FolderNameList[i], PostFix, FileNameList);
	}
}

#include "stdafx.h"

#include "BatchLocalization.h"

#include "ReadWriteTiffImage.h"

// call th_OnlineLocalizationLD etc.

wstring BatchProc_FolderName;
wstring BatchProc_FileName_Prefix;
wstring BatchProc_FileName_Postfix;
wstring BatchProc_SavePath;
int BatchProc_MergeLocEn = 0;

int TotalFileNum = 0;
int CurFileCnt = 0;

#include "tinytiffreader.h"
#include "tinytiffwriter.h"

#include "FileNameProcess.h"

// allocate resource
// feed images
// save localizations

UINT th_BatchLocalization(LPVOID params)
{
	IsBatchLocRunning = true;

	vector<wstring>  FileNameList;
	vector<wstring>  FileNameList_Sorted;

	SearchFilesInDir(BatchProc_FolderName, BatchProc_FileName_Prefix, BatchProc_FileName_Postfix, FileNameList);

	FilesSort(FileNameList, FileNameList_Sorted); // sort file name


	TotalFileNum = FileNameList_Sorted.size();


	printf("File number:%d\n", TotalFileNum);

	if (TotalFileNum <= 0)return 0;


	char TiffImgName_Buff[1024];

	int StartBatchEn = 0;
	int FinishBatchEn = 0;

	CString SRImageName;

	LocPara_Global.TotalFrameNum = 0;

	for (int fcnt = 0; fcnt < TotalFileNum; fcnt++)
	{
		ImageName = FileNameList_Sorted[fcnt].c_str();


		CurFileCnt = fcnt + 1;

		wprintf(L"Cur file: %d in %d, name: %s\n", CurFileCnt, TotalFileNum, ImageName);

		ImgDataQueue.clear();

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

		// don't loc if there is only 1 image
		if (TotalFrame <= 1)
		{
			TinyTIFFReader_close(tiffr);
			continue;
		}

		// image size of batch loc doesn't correct
		if (BatchProc_MergeLocEn)
		{
			if (fcnt > 0)
			{
				if ((LocPara_Global.ImageWidth != ImageWidth) || (LocPara_Global.ImageHigh != ImageHigh))
				{
					printf("\nbatch loc image size error\n");
					break;
				}
			}
		}

		// update sr image size
		LocPara_Global.ImageWidth = ImageWidth;
		LocPara_Global.ImageHigh = ImageHigh;
		LocPara_Global.UpdateSRImageSize();

		if (BatchProc_MergeLocEn)
		{
			LocPara_Global.TotalFrameNum += TotalFrame;
		}
		else
		{
			LocPara_Global.TotalFrameNum = TotalFrame;
		}

		// start batch localization, allocate resource and start parallel threads

		if (BatchProc_MergeLocEn)
		{
			if (fcnt == 0)StartBatchEn = 1;
			else StartBatchEn = 0;
		}
		else
		{
			StartBatchEn = 1;
		}

		if (StartBatchEn)
		{
			SetLocDataFileName_BatchProc(); // set proper localization data save name
			StartLocalizationThread();
		}


		// feed images, just as imagej to cpp dll
		int BatchedImgNum = (2048 * 2048 / ImageWidth / ImageHigh);
		BatchedImgNum = (BatchedImgNum / 2) * 2;

		BatchedImgNum = max(BatchedImgNum, 1);

		int GroupNum = (TotalFrame + BatchedImgNum - 1) / BatchedImgNum;;


		unsigned short *RawImgBuf = new unsigned short[BatchedImgNum*ImageWidth*ImageHigh];



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
			for (int i = 0; i < CurImgNum; i++)
			{
				// read images
				TinyTIFFReader_getSampleData(tiffr, &RawImgBuf[i*ImageWidth*ImageHigh], 0); // get image data
				TinyTIFFReader_readNext(tiffr);

			}

			qImgData CurImgInf;
			CreateFeedImgMemory(CurImgInf, RawImgBuf, LocPara_Global.ImageWidth, LocPara_Global.ImageHigh, CurImgNum);

			ImgDataQueue.push(CurImgInf);

		}


		while (!ImgDataQueue.empty())
		{
			Sleep(2);
		}

		delete [] RawImgBuf;

		// finish localization
		TinyTIFFReader_close(tiffr);



		// save super-resolution image
		if (BatchProc_MergeLocEn)
		{
			if (fcnt == 0)SRImageName = ImageName;
		}
		else
		{
			SRImageName = ImageName;
		}

		if (BatchProc_MergeLocEn)
		{
			if (fcnt == TotalFileNum - 1)FinishBatchEn = 1;
			else FinishBatchEn = 0;
		}
		else
		{
			FinishBatchEn = 1;
		}

		if (FinishBatchEn)
		{
			StopLocThread();

			// wait localization thread to finish
			while (IsLocRunning == true)
			{
				Sleep(2);
			}


			printf("Writing super-resolution image...\n");

			FormatBatchLocSRImageName(SRImageName);

			if (LocPara_Global.LocType == LocType_GS2D)
			{
				BatchProc_WriteSRImage2D(SRImageName);
			}
			else
			{
				BatchProc_WriteSRImage3D(SRImageName);
			}

			// release resources
			FinishLocalizationThread();

		}
	}


	IsBatchLocRunning = false;


	printf("Batch localization finish\n");

	return 0;
}

void FormatBatchLocSRImageName(CString & SRImageName)
{
	// format sr image name
	SRImageName.TrimRight(L".tif");


	CString MultiFitStr;
	if (LocPara_Global.MultiEmitterFitEn) MultiFitStr = "_M";
	else MultiFitStr = "_S";

	CString ConsecFitStr;
	if (LocPara_Global.ConsecFitEn) ConsecFitStr.Format(L"_Consec%.0fnm", LocPara_Global.ConsecFit_DistanceTh_nm);
	else ConsecFitStr = "";

	SRImageName.Format(L"%s_result%dD%d%s%s_rend%.2fnm.tif", SRImageName, LocPara_Global.LocType + 2, LocPara_Global.ROISize, MultiFitStr, ConsecFitStr, LocPara_Global.PixelSize / LocPara_Global.PixelZoom);

	// modify save path
	SRImageName.TrimLeft(BatchProc_FolderName.c_str());
	SRImageName.TrimLeft(L"\\");
	SRImageName.Replace(L"\\", L"__");

	CString SavePath = BatchProc_SavePath.c_str();
	SavePath = SavePath + L"\\";
	SRImageName = SavePath + SRImageName;


}

void BatchProc_WriteSRImage2D(CString SRImageName)
{
	float* h_RendFloatImage2D = new float[LocPara_Global.SRImageWidth*LocPara_Global.SRImageHigh];


	// create 2d image 
	cudaMemcpyAsync(h_RendFloatImage2D, RendData.d_SRIntensityImg, LocPara_Global.SRImageWidth* LocPara_Global.SRImageHigh * sizeof(float), cudaMemcpyDeviceToHost, loc_stream1);
	WaitGPUStream(loc_stream1);


#define BufLen		1024
	char SRImgFileName[BufLen];
	//		sprintf_s(SRImgFileName, "ReRendered sr image_%d-%dnm.tif", GroupCnt, (int)RenderingPixelSize);

	WideCharToMultiByte(CP_OEMCP, NULL, SRImageName.GetBuffer(), -1, SRImgFileName, BufLen, NULL, FALSE);

	TinyTIFFFile* tif = TinyTIFFWriter_open(SRImgFileName, 32, LocPara_Global.SRImageWidth, LocPara_Global.SRImageHigh); // 32 bit float image
	TinyTIFFWriter_writeImage(tif, h_RendFloatImage2D);

	TinyTIFFWriter_close(tif);


	delete[]h_RendFloatImage2D;

}

void BatchProc_WriteSRImage3D(CString SRImageName)
{
	RGBImage SRImg_RGB((unsigned char *)RendData.h_SaveRendImg, LocPara_Global.SRImageWidth, LocPara_Global.SRImageHigh);

	WriteRGBImage(SRImageName, &SRImg_RGB);
}


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


#include "stdafx.h"
#include "QC_STORM_CPPDLL.h"
#include "OnlineLocalizationLD.h"

#include "OnlineLocalizationLD.h"
#include "StatInfDisplay.h"

#include "PostProcess.h"

#include "BatchLocalization.h"

/*
* Class:     QC_STORM_
* Method:    lm_SetImagePara
* Signature: (III[C)V
*/
JNIEXPORT void JNICALL Java_QC_1STORM_1_lm_1SetImagePara
(JNIEnv *env, jclass obj, jint iImageWidth, jint iImageHigh, jint SRImageWidth, jint SRImageHigh, jint iFrameNum, jstring iJImageName)
{


	LocPara_Global.ImageWidth = iImageWidth;
	LocPara_Global.ImageHigh = iImageHigh;
	LocPara_Global.SRImageWidth = SRImageWidth;
	LocPara_Global.SRImageHigh = SRImageHigh;

	LocPara_Global.TotalFrameNum = iFrameNum;



	// access string name
	wchar_t * str_java = (wchar_t *)env->GetStringChars(iJImageName, NULL);

	ImageName = str_java;
	env->ReleaseStringChars(iJImageName, (jchar*)str_java);

}

/*
* Class:     QC_STORM_
* Method:    lm_SetLocPara
* Signature: (FFFIIIFFF)V
*/
JNIEXPORT void JNICALL Java_QC_1STORM_1_lm_1SetLocPara
(JNIEnv *env, jclass obj, jfloat KAdc, jfloat Offset, jfloat QE, jint ROISizeI, jint LocTypeI, jint MultiEmitterFitEn, jint WLEEn, jint ConsecutiveFitEnI, jfloat ConsecFilterRadiusF, jfloat RawPixelSizeF, jfloat RenderPixelZoomF, jfloat SNR_th)
{
	// for both 2d and 3d
	LocPara_Global.KAdc = KAdc;
	LocPara_Global.Offset = Offset;
	LocPara_Global.QE = QE;

	LocPara_Global.ROISize = ROISizeI;
	LocPara_Global.LocType = LocTypeI;

	LocPara_Global.MultiEmitterFitEn = MultiEmitterFitEn;

	LocPara_Global.WLEEn = WLEEn;


	LocPara_Global.ConsecFit_DistanceTh_nm = ConsecFilterRadiusF;
	LocPara_Global.ConsecFitEn = ConsecutiveFitEnI;

	LocPara_Global.PixelSize = RawPixelSizeF;
	LocPara_Global.PixelZoom = RenderPixelZoomF;

	// image rendering
	LocPara_Global.SNR_th = SNR_th;

}
/*
* Class:     QC_STORM_
* Method:    lm_SetLocPara3D
* Signature: (FFFFFFFFI)V
*/
JNIEXPORT void JNICALL Java_QC_1STORM_1_lm_1SetLocPara3D
(JNIEnv *env, jclass obj, jfloat MinZDepthF, jfloat MaxZDepthF, jfloat ZDepthCorrFactor, jfloat p4_XGY, jfloat p3_XGY, jfloat p2_XGY, jfloat p1_XGY, jfloat p0_XGY, jfloat p4_XLY, jfloat p3_XLY, jfloat p2_XLY, jfloat p1_XLY, jfloat p0_XLY)
{
	// astigmatism and double helix
	LocPara_Global.MinZDepth = MinZDepthF;
	LocPara_Global.MaxZDepth = MaxZDepthF;

	LocPara_Global.ZDepthCorrFactor = ZDepthCorrFactor;


	LocPara_Global.p4_XGY = p4_XGY;
	LocPara_Global.p3_XGY = p3_XGY;
	LocPara_Global.p2_XGY = p2_XGY;
	LocPara_Global.p1_XGY = p1_XGY;
	LocPara_Global.p0_XGY = p0_XGY;

	LocPara_Global.p4_XLY = p4_XLY;
	LocPara_Global.p3_XLY = p3_XLY;
	LocPara_Global.p2_XLY = p2_XLY;
	LocPara_Global.p1_XLY = p1_XLY;
	LocPara_Global.p0_XLY = p0_XLY;


	LocPara_Global.ColorMode_3D = ImgRend_ColorMode_3D_BlueToRed;

}

/*
* Class:     QC_STORM_
* Method:    lm_SetLocParaDH3D
* Signature: (I)V
*/
JNIEXPORT void JNICALL Java_QC_1STORM_1_lm_1SetLocParaDH3D
(JNIEnv *, jclass, jint RotationType, jfloat MeanDistance, jfloat DistanceTh)
{
	LocPara_Global.DH_RotateType = RotationType;
	LocPara_Global.DH_MeanDistance = MeanDistance;
	LocPara_Global.DH_DistanceTh = DistanceTh;

}


/*
* Class:     hust_whno_SMLM_QC_STORM_Plug
* Method:    lm_SetStatInfSelection
* Signature: (I)V
*/
JNIEXPORT void JNICALL Java_QC_1STORM_1_lm_1SetStatInfSelection
(JNIEnv *env, jclass obj, jint DispSel, jint SpatialResolutionEn)
{
	StatDispSel = DispSel;

	LocPara_Global.SpatialResolutionCalcEn = SpatialResolutionEn;

}

/*
* Class:     QC_STORM_
* Method:    lm_StartLocThread
* Signature: ()V
*/
JNIEXPORT void JNICALL Java_QC_1STORM_1_lm_1StartLocThread
(JNIEnv *env, jclass obj)
{
	OpenConsole(); // open console window to display printf

	IsBatchLocRunning = false;
	IsLocRunning = true;

	SetLocDataFileName(); // set proper localization data save name

	StartLocalizationThread();

}

/*
* Class:     QC_STORM_
* Method:    lm_StopLocThread
* Signature: ()V
*/
JNIEXPORT void JNICALL Java_QC_1STORM_1_lm_1StopLocThread
(JNIEnv *env, jclass obj)
{

	FinishLocalizationThread();

}

/*
* Class:     QC_STORM_
* Method:    lm_GetMaxDispVal
* Signature: ()I
*/
JNIEXPORT jint JNICALL Java_QC_1STORM_1_lm_1GetMaxDispVal
(JNIEnv *env, jclass obj)
{
	int DispMaxVal = 1;

	if (IsLocResourceAllocated == 1)
	{
		DispMaxVal = RendData.GetDispMaxVal();
	}
	
	return DispMaxVal;
}

/*
* Class:     QC_STORM_
* Method:    lm_FeedImageData
* Signature: ([S)V
*/
JNIEXPORT void JNICALL Java_QC_1STORM_1_lm_1FeedImageData
(JNIEnv *env, jclass obj, jshortArray jImgArry, jint FrameNumI)
{
	// jni environment
	int len = (*env).GetArrayLength(jImgArry);

	jboolean iscopy = false;
	jshort * elems = (jshort *)env->GetPrimitiveArrayCritical(jImgArry, &iscopy);


	// cur image info
	int BatchedImgSize = FrameNumI*LocPara_Global.ImageWidth*LocPara_Global.ImageHigh;

	unsigned int MaxBufferImageNum = 200 * 2048 * 2048 / BatchedImgSize;

	while (ImgDataQueue.unsafe_size() > MaxBufferImageNum)
	{
		Sleep(2);
	}


	qImgData CurImgInf;

	CreateFeedImgMemory(CurImgInf, (unsigned short*)elems, LocPara_Global.ImageWidth, LocPara_Global.ImageHigh, FrameNumI);


	//
	ImgDataQueue.push(CurImgInf);


	//	env->ReleaseShortArrayElements(jImgArry, elems, 0);
	env->ReleasePrimitiveArrayCritical(jImgArry, elems, JNI_ABORT);

}

/*
* Class:     QC_STORM_
* Method:    lm_GetSMLMImage
* Signature: ()[F
*/
JNIEXPORT jfloatArray JNICALL Java_QC_1STORM_1_lm_1GetSMLMImage
(JNIEnv *env, jclass obj)
{
	int ImgSize;
	jfloatArray result = NULL;


	// for 2d imaging
	if (IsLocResourceAllocated == 1)
	{
		// disp image
		ImgSize = LocPara_Global.SRImageWidth* LocPara_Global.SRImageHigh;
		result = env->NewFloatArray(ImgSize);
		env->SetFloatArrayRegion(result, 0, ImgSize, h_RendFloatImage2D);

	}
	else
	{
		// puodo
		ImgSize = 4;
		result = env->NewFloatArray(4);

	}

	return result;
}
/*
* Class:     QC_STORM_
* Method:    lm_GetSMLMImage3D
* Signature: ()[I
*/
JNIEXPORT jintArray JNICALL Java_QC_1STORM_1_lm_1GetSMLMImage3D
(JNIEnv *env, jclass obj)
{
	int ImgSize;
	jintArray result = NULL;


	// for 3d imaging
	if (IsLocResourceAllocated == 1)
	{
		// disp image
		ImgSize = LocPara_Global.SRImageWidth* LocPara_Global.SRImageHigh;
		result = env->NewIntArray(ImgSize);
		env->SetIntArrayRegion(result, 0, ImgSize, (jint*)(RendData.h_SaveRendImg));

	}
	else
	{
		// puodo
		ImgSize = 4;
		result = env->NewIntArray(4);

	}


	return result;
}

/*
* Class:     hust_whno_SMLM_QC_STORM_Plug
* Method:    lm_SetSpatialResolutionInf
* Signature: (I)V
*/
JNIEXPORT void JNICALL Java_QC_1STORM_1_lm_1SetSpatialResolutionInf
(JNIEnv *env, jclass obj, jfloat iStructureSize)
{
	LocPara_Global.StrucuteSize_2D = iStructureSize;

}


/*
* Class:     QC_STORM_
* Method:    lm_GetStatInfImageSize
* Signature: (I)[I
*/
JNIEXPORT jintArray JNICALL Java_QC_1STORM_1_lm_1GetStatInfImageSize
(JNIEnv *env, jclass obj)
{

	jintArray result = env->NewIntArray(2);

	int ImgInf[2];

	ImgInf[0] = AxesImgWidth;
	ImgInf[1] = AxesImgHigh;

	env->SetIntArrayRegion(result, 0, 2, (jint*)(ImgInf));

	return result;
}



/*
* Class:     QC_STORM_
* Method:    lm_GetStatInfImage
* Signature: (I)[I
*/
JNIEXPORT jintArray JNICALL Java_QC_1STORM_1_lm_1GetStatInfImage
(JNIEnv *env, jclass obj, jint n)
{
	int ImgSize;
	jintArray result = NULL;

	if (IsLocResourceAllocated == 1)
	{
		// stat inf disp image
		ImgSize = AxesImgWidth* AxesImgHigh;
		result = env->NewIntArray(ImgSize);


		UpdateStatInfDisplay();

		switch (n)
		{
		case 0:			
			ConvertCImgToImageJ(ConvertInfImageBuf, StatInfDisplay.InfDisp_Hist_TotalPhoton->GetAxisImageData(), AxesImgWidth, AxesImgHigh);
			break;

		case 1:
			ConvertCImgToImageJ(ConvertInfImageBuf, StatInfDisplay.InfDisp_Hist_LocPrec->GetAxisImageData(), AxesImgWidth, AxesImgHigh);
			break;

		case 2:
			ConvertCImgToImageJ(ConvertInfImageBuf, StatInfDisplay.InfDisp_Hist_SNR->GetAxisImageData(), AxesImgWidth, AxesImgHigh);
			break;

		case 3:
			ConvertCImgToImageJ(ConvertInfImageBuf, StatInfDisplay.InfDisp_Hist_PSFWidth->GetAxisImageData(), AxesImgWidth, AxesImgHigh);
			break;

		// elapse by time
		case 4:
			ConvertCImgToImageJ(ConvertInfImageBuf, StatInfDisplay.InfDisp_Curve_TotalPhoton->GetAxisImageData(), AxesImgWidth, AxesImgHigh);
			break;

		case 5:
			ConvertCImgToImageJ(ConvertInfImageBuf, StatInfDisplay.InfDisp_Curve_LocPrec->GetAxisImageData(), AxesImgWidth, AxesImgHigh);
			break;
			
		case 6:
			ConvertCImgToImageJ(ConvertInfImageBuf, StatInfDisplay.InfDisp_Curve_Ontime->GetAxisImageData(), AxesImgWidth, AxesImgHigh);
			break;

		case 7:
			ConvertCImgToImageJ(ConvertInfImageBuf, StatInfDisplay.InfDisp_Curve_SNR->GetAxisImageData(), AxesImgWidth, AxesImgHigh);
			break;

		case 8:
			ConvertCImgToImageJ(ConvertInfImageBuf, StatInfDisplay.InfDisp_Curve_PSFWidth->GetAxisImageData(), AxesImgWidth, AxesImgHigh);
			break;

		case 9:
			ConvertCImgToImageJ(ConvertInfImageBuf, StatInfDisplay.InfDisp_Curve_LocDensity2D->GetAxisImageData(), AxesImgWidth, AxesImgHigh);
			break;

		case 10:
			ConvertCImgToImageJ(ConvertInfImageBuf, StatInfDisplay.InfDisp_Curve_Background->GetAxisImageData(), AxesImgWidth, AxesImgHigh);
			break;

		// spatial resolution variation
		case 11:
			ConvertCImgToImageJ(ConvertInfImageBuf, StatInfDisplay.InfDisp_Curve_SpatialResolution->GetAxisImageData(), AxesImgWidth, AxesImgHigh);
			break;
		case 12:
			ConvertCImgToImageJ(ConvertInfImageBuf, StatInfDisplay.InfDisp_Curve_NyquistResolution->GetAxisImageData(), AxesImgWidth, AxesImgHigh);
			break;
		case 13:
			ConvertCImgToImageJ(ConvertInfImageBuf, StatInfDisplay.InfDisp_Curve_DimensionFD->GetAxisImageData(), AxesImgWidth, AxesImgHigh);
			break;
		case 14:
			ConvertCImgToImageJ(ConvertInfImageBuf, StatInfDisplay.InfDisp_Curve_LocDensityFD->GetAxisImageData(), AxesImgWidth, AxesImgHigh);
			break;

		default:

			break;
		}
		
		
		env->SetIntArrayRegion(result, 0, ImgSize, (jint*)(ConvertInfImageBuf));

	}
	else
	{
		// puodo
		ImgSize = 4;
		result = env->NewIntArray(4);

	}

	return result;
}

/*
* Class:     QC_STORM_
* Method:    lm_IsLocFinish
* Signature: ()I
*/
JNIEXPORT jint JNICALL Java_QC_1STORM_1_lm_1IsLocFinish
(JNIEnv *env, jclass obj)
{
	int IsFinish = (IsLocRunning == false) && (IsBatchLocRunning == false);


	return IsFinish;
}


// rerend image

/*
* Class:     QC_STORM_
* Method:    lm_SetRerendImagePara
* Signature: ([C)V
*/
JNIEXPORT void JNICALL Java_QC_1STORM_1_lm_1SetRerendImagePara
(JNIEnv *env, jclass obj, jstring DataPath, jint IsDriftCorrectionI, jint DriftCorrGroupFrameNumI)
{

	// access string name
	wchar_t * str_java = (wchar_t *)env->GetStringChars(DataPath, NULL);

	RerendDataPath = str_java;
	env->ReleaseStringChars(DataPath, (jchar*)str_java);

	IsDriftCorrection = IsDriftCorrectionI;
	DriftCorrGroupFrameNum = DriftCorrGroupFrameNumI;

}
/*
* Class:     QC_STORM_
* Method:    lm_StartRerend
* Signature: ()V
*/
JNIEXPORT void JNICALL Java_QC_1STORM_1_lm_1StartRerend
(JNIEnv *, jclass)
{

	RerendProgress = 0;

	IsLocRunning = true;
	OnlineRendAlive = true;

	IsBatchLocRunning = false;


	AfxBeginThread(th_RerendImage, NULL);

	AfxBeginThread(th_OnlineRendDispLD, NULL);


}
/*
* Class:     QC_STORM_
* Method:    lm_GetRerendImageInf
* Signature: ()[I
*/
JNIEXPORT jintArray JNICALL Java_QC_1STORM_1_lm_1GetRerendImageInf
(JNIEnv *env, jclass obj)
{

#define RerendImageInfLen	8

	jintArray result = NULL;

	result = env->NewIntArray(RerendImageInfLen);

	int ImgInf[RerendImageInfLen];

	ImgInf[0] = RerendProgress;
	ImgInf[1] = LocPara_Global.ImageWidth;
	ImgInf[2] = LocPara_Global.ImageHigh;
	ImgInf[3] = LocPara_Global.SRImageWidth;
	ImgInf[4] = LocPara_Global.SRImageHigh;


	env->SetIntArrayRegion(result, 0, RerendImageInfLen, (jint*)(ImgInf));

	return result;

}

/*
* Class:     QC_STORM_
* Method:    lm_ReleaseRerendResource
* Signature: ()V
*/
JNIEXPORT void JNICALL Java_QC_1STORM_1_lm_1ReleaseRerendResource
(JNIEnv *, jclass)
{

	DeinitAllLocResource(1);

}

/*
* Class:     QC_STORM_
* Method:    lm_StartBatchImageLoc
* Signature: (III[C)V
*/
JNIEXPORT void JNICALL Java_QC_1STORM_1_lm_1StartBatchImageLoc
(JNIEnv *env, jclass obj, jstring ImageFolderPath, jstring FileExtension, jstring ResultsPath)
{
	
	OpenConsole();

	// access string name
	wchar_t * str_java = (wchar_t *)env->GetStringChars(ImageFolderPath, NULL);

	BatchProc_FolderName = str_java;
	env->ReleaseStringChars(ImageFolderPath, (jchar*)str_java);


	// access string name
	str_java = (wchar_t *)env->GetStringChars(FileExtension, NULL);

	BatchProc_Postfix = str_java;
	env->ReleaseStringChars(FileExtension, (jchar*)str_java);

	// access string name
	str_java = (wchar_t *)env->GetStringChars(ResultsPath, NULL);

	BatchProc_SavePath = str_java;
	env->ReleaseStringChars(ResultsPath, (jchar*)str_java);

	//
	IsBatchLocRunning = true;
	IsLocRunning = true;

	AfxBeginThread(th_BatchLocalization, NULL);


//	wprintf(L"%s %s\n", BatchProc_FolderName.c_str(), BatchProc_Postfix.c_str());
}


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

// OnlineSMLM CPPDLL.cpp : 定义 DLL 应用程序的导出函数。
//

#include "stdafx.h"
#include "QC_STORM_CPPDLL.h"


#include "StatInfDisplay.h"

#include "OnlineFeedback_rc.h"

#include "LocResource_ZDriftMeasure.h"

#define DLL_CALL_DEBUG		0



/*
* Class:     hust_whno_SMLM_QC_STORM_Plug
* Method:    lm_SetProcessorID
* Signature: ([C)V
*/
JNIEXPORT void JNICALL Java_hust_whno_SMLM_QC_1STORM_1Plug_lm_1SetProcessorID
(JNIEnv *, jclass, jint id)
{
	ProcessorID = id;
}

/*
* Class:     hust_whno_SMLM_QC_STORM_Plug
* Method:    lm_GetProcessorID
* Signature: ([C)V
*/
JNIEXPORT jint JNICALL Java_hust_whno_SMLM_QC_1STORM_1Plug_lm_1GetProcessorID
(JNIEnv *, jclass)
{
	return ProcessorID;
}

/*
* Class:     hust_whno_SMLM_QC_STORM_Plug
* Method:    lm_SetSavePath
* Signature: ([C)V
*/
JNIEXPORT void JNICALL Java_hust_whno_SMLM_QC_1STORM_1Plug_lm_1SetSavePath
(JNIEnv *env, jclass obj, jcharArray JSavePath)
{
#if DLL_CALL_DEBUG
	printf("lm_SetSavePath called\n");
#endif // DLL_CALL_DEBUG

	// jni environment
	int len = (*env).GetArrayLength(JSavePath);
	jboolean IsCopy = true;
	jchar * elems = (*env).GetCharArrayElements(JSavePath, &IsCopy);

	wchar_t *NameBuf = new wchar_t[200];
	memcpy(NameBuf, elems, 2 * len);
	NameBuf[len + 0] = 0;
	NameBuf[len + 1] = 0;


	ResultSavePath = NameBuf;


	delete[] NameBuf;

	env->ReleaseCharArrayElements(JSavePath, elems, 0);

}



/*
* Class:     hust_whno_SMLM_QC_STORM_Plug
* Method:    lm_SetLocPara
* Signature: (FFFIIIFFF)V
*/
JNIEXPORT void JNICALL Java_hust_whno_SMLM_QC_1STORM_1Plug_lm_1SetLocPara
(JNIEnv *env, jclass obj, jfloat KAdc, jfloat Offset, jfloat QE, jint ROISizeI, jint LocTypeI, jint MultiEmitterFitEn, jint WLEEn, jint ConsecutiveFitEnI, jfloat ConsecFilterRadiusF, jfloat RawPixelSizeF, jfloat RenderPixelZoomF, jfloat SNR_th)
{
	OpenConsole();

#if DLL_CALL_DEBUG
	printf("lm_SetLocPara called\n");
#endif // DLL_CALL_DEBUG

	//
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
* Class:     hust_whno_SMLM_QC_STORM_Plug
* Method:    lm_SetLocPara3D
* Signature: (FFFFFFFFI)V
*/
JNIEXPORT void JNICALL Java_hust_whno_SMLM_QC_1STORM_1Plug_lm_1SetLocPara3D
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
* Class:     hust_whno_SMLM_QC_STORM_Plug
* Method:    lm_SetImagePara
* Signature: (IIII)V
*/
JNIEXPORT void JNICALL Java_hust_whno_SMLM_QC_1STORM_1Plug_lm_1SetImagePara
(JNIEnv *env, jclass obj, jint iImageWidth, jint iImageHigh, jint SRImageWidth, jint SRImageHigh)
{
	LocPara_Global.ImageWidth = iImageWidth;
	LocPara_Global.ImageHigh = iImageHigh;
	LocPara_Global.SRImageWidth = SRImageWidth;
	LocPara_Global.SRImageHigh = SRImageHigh;

}


/*
* Class:     hust_whno_SMLM_QC_STORM_Plug
* Method:    lm_SetAcquisitionPara
* Signature: (I[C)V
*/
JNIEXPORT void JNICALL Java_hust_whno_SMLM_QC_1STORM_1Plug_lm_1SetAcquisitionPara
(JNIEnv *env, jclass obj, jcharArray iCreateTimeIdxS)
{
#if DLL_CALL_DEBUG
	printf("lm_SetAcquisitionPara called\n");
#endif // DLL_CALL_DEBUG

	// jni environment
	int len = (*env).GetArrayLength(iCreateTimeIdxS);
	jboolean IsCopy = true;
	jchar * elems = (*env).GetCharArrayElements(iCreateTimeIdxS, &IsCopy);

	wchar_t *NameBuf = new wchar_t[200];
	memcpy(NameBuf, elems, 2 * len);
	NameBuf[len + 0] = 0;
	NameBuf[len + 1] = 0;


	CreateTimeIdxStr = NameBuf;

	delete[] NameBuf;

	env->ReleaseCharArrayElements(iCreateTimeIdxS, elems, 0);

}


/*
* Class:     hust_whno_SMLM_QC_STORM_Plug
* Method:    lm_SetStatInfSelection
* Signature: (I)V
*/
JNIEXPORT void JNICALL Java_hust_whno_SMLM_QC_1STORM_1Plug_lm_1SetStatInfSelection
(JNIEnv *env, jclass obj, jint DispSel, jint SpatialResolutionEn)
{
	StatInfDispSel = DispSel;
	
	LocPara_Global.SpatialResolutionCalcEn = SpatialResolutionEn;

	RenderingState.MakeAProcess();
}

/*
* Class:     hust_whno_SMLM_QC_STORM_Plug
* Method:    lm_FeedImageData
* Signature: ([S)V
*/
JNIEXPORT void JNICALL Java_hust_whno_SMLM_QC_1STORM_1Plug_lm_1FeedImageData
(JNIEnv *env, jclass obj, jshortArray jImgArry, jint FrameNumI)
{
#if DLL_CALL_DEBUG
	printf("lm_FeedImageData called\n");
#endif // DLL_CALL_DEBUG

	// jni environment
	int len = (*env).GetArrayLength(jImgArry);

	jboolean iscopy = false;
	jshort * elems = (jshort *)env->GetPrimitiveArrayCritical(jImgArry, &iscopy);



	qImgData CurImgInf;

	// cur image info
	FrameNumI = 1;
	int BatchedImgSize = FrameNumI*LocPara_Global.ImageWidth*LocPara_Global.ImageHigh;


	CurImgInf.BatchFrameNum = FrameNumI;


	cudaError_t err = cudaErrorMemoryAllocation;

	// try allocate CPU memory
	err = cudaMallocHost((void **)&CurImgInf.pImgData, BatchedImgSize * sizeof(short));


	if (err == cudaSuccess)
	{
		//	printf("cuda suc:%s\n", str);
		CurImgInf.ImageSource = ImageSource_CPU_Pinned;
	}
	else
	{
		while (1)
		{
			// if GPU memory is allocate error, then try CPU memory
			CurImgInf.pImgData = new unsigned short[BatchedImgSize];
			CurImgInf.ImageSource = ImageSource_CPU_Normal;

			if (CurImgInf.pImgData != NULL)break;
		}
	}

	memcpy(CurImgInf.pImgData, elems, BatchedImgSize * sizeof(short));

	//
	ImgDataQueue.push(CurImgInf);


	env->ReleasePrimitiveArrayCritical(jImgArry, elems, JNI_ABORT);

	// 
	if (IsLocRunning == false)
	{
		IsLocRunning = true;

		Java_hust_whno_SMLM_QC_1STORM_1Plug_lm_1StartLocThread(env, obj);
	}

}

/*
* Class:     hust_whno_SMLM_QC_STORM_Plug
* Method:    lm_StartLocThread
* Signature: ()V
*/
JNIEXPORT void JNICALL Java_hust_whno_SMLM_QC_1STORM_1Plug_lm_1StartLocThread
(JNIEnv *env, jclass obj)
{
#if DLL_CALL_DEBUG
	printf("lm_StartLocThread called %d %d\n", OnlineLocAlive, IsLocRunning);
#endif // DLL_CALL_DEBUG

	if (OnlineLocAlive == false)
	{

		SetAllOnlineThreadAlive();

		InitAllLocResource();


		// start online rendering and feedback threads
		AfxBeginThread(th_OnlineRendDispLD, NULL);
		AfxBeginThread(th_OnlineFeedback, NULL);

		// later to start localization thread
		AfxBeginThread(th_OnlineLocalizationLD, NULL);

	}
}

/*
* Class:     hust_whno_SMLM_QC_STORM_Plug
* Method:    lm_StopLocThread
* Signature: ()V
*/
JNIEXPORT void JNICALL Java_hust_whno_SMLM_QC_1STORM_1Plug_lm_1StopLocThread
(JNIEnv *env, jclass obj)
{
#if DLL_CALL_DEBUG
	printf("lm_StopLocThread called\n");
#endif // DLL_CALL_DEBUG


	// stop only online loc, others will be stoped by online loc
	OnlineLocAlive = false;

}

/*
* Class:     hust_whno_SMLM_QC_STORM_Plug
* Method:    lm_ReleaseLocResource
* Signature: ()V
*/
JNIEXPORT void JNICALL Java_hust_whno_SMLM_QC_1STORM_1Plug_lm_1ReleaseLocResource
(JNIEnv *env, jclass obj)
{
#if DLL_CALL_DEBUG
	printf("lm_ReleaseLocResource called\n");
#endif // DLL_CALL_DEBUG

	ClearAllOnlineThreadAlive();

	DeinitAllLocResource();

}

/*
* Class:     hust_whno_SMLM_QC_STORM_Plug
* Method:    lm_IsLocFinish
* Signature: ()I
*/
JNIEXPORT jint JNICALL Java_hust_whno_SMLM_QC_1STORM_1Plug_lm_1IsLocFinish
(JNIEnv *env, jclass obj)
{
#if DLL_CALL_DEBUG
	printf("lm_IsLocFinish called %d %d\n", OnlineLocAlive, IsLocRunning);
#endif // DLL_CALL_DEBUG

	int IsFinish = IsLocRunning == false;

	return IsFinish;
}


/*
* Class:     hust_whno_SMLM_QC_STORM_Plug
* Method:    lm_GetMaxDispVal
* Signature: ()I
*/
JNIEXPORT jint JNICALL Java_hust_whno_SMLM_QC_1STORM_1Plug_lm_1GetMaxDispVal
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
* Class:     hust_whno_SMLM_QC_STORM_Plug
* Method:    lm_GetWaitImageNum
* Signature: ()I
*/
JNIEXPORT jint JNICALL Java_hust_whno_SMLM_QC_1STORM_1Plug_lm_1GetWaitImageNum
(JNIEnv *env, jclass obj)
{
	int WaitImgNum;

	WaitImgNum = ImgDataQueue.unsafe_size();

	return WaitImgNum;
}

/*
* Class:     hust_whno_SMLM_QC_STORM_Plug
* Method:    lm_GetSMLMImage
* Signature: ()[F
*/
JNIEXPORT jfloatArray JNICALL Java_hust_whno_SMLM_QC_1STORM_1Plug_lm_1GetSMLMImage
(JNIEnv *env, jclass obj)
{

	int ImgSize;
	jfloatArray result=NULL;


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
* Class:     hust_whno_SMLM_QC_STORM_Plug
* Method:    lm_GetSMLMImage3D
* Signature: ()[I
*/
JNIEXPORT jintArray JNICALL Java_hust_whno_SMLM_QC_1STORM_1Plug_lm_1GetSMLMImage3D
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
JNIEXPORT void JNICALL Java_hust_whno_SMLM_QC_1STORM_1Plug_lm_1SetSpatialResolutionInf
(JNIEnv *env, jclass obj, jint iFramePerGroup, jint iIsHollowTube, jfloat iStructureSize, jfloat RSCResolutionTh)
{
	LocPara_Global.ImagesPerGroup = iFramePerGroup;
	LocPara_Global.StrucuteSize_2D = iStructureSize;
	LocPara_Global.IsHollowTube = iIsHollowTube;
	LocPara_Global.RSCResolutionTh = RSCResolutionTh;

}

/*
* Class:     hust_whno_SMLM_QC_STORM_Plug
* Method:    float lm_GetCurSpatialResolution
* Signature: ()I
*/
JNIEXPORT jfloat JNICALL Java_hust_whno_SMLM_QC_1STORM_1Plug_lm_1GetCurSpatialResolution
(JNIEnv *env, jclass obj)
{
	float CurResolution = 1000;

	if (OnlineLocAlive)
	{
		CurResolution = SpatialResolutionCalc.GetCurSpatialResolution();
	}

	return CurResolution;
}

/*
* Class:     hust_whno_SMLM_QC_STORM_Plug
* Method:    float lm_GetMeanLocPrec
* Signature: ()I
*/
JNIEXPORT jfloat JNICALL Java_hust_whno_SMLM_QC_1STORM_1Plug_lm_1GetMeanLocPrec
(JNIEnv *env, jclass obj)
{
	float Mean_LocPrecisionXY = 0;

	if (OnlineLocAlive)
	{
		Mean_LocPrecisionXY = FluoStatisticData_TypeDef::GetTimeVaryMean(FluoStatData.TimeVary_LocPrecisionXY);
	}

	return Mean_LocPrecisionXY;
}


// feedback control port set

int SetUARTPort(serial::Serial & UARTPort, jint UARTID, jint BaudRateI, jint IsEnable)
{

	char buf[64];
	sprintf(buf, "COM%d", UARTID);

	string buf_str = buf;

	int IsValidPort = IsSerialPortValid(buf_str);

	printf("set uart:%d %d %d\n", UARTID, IsValidPort, BaudRateI);



	if (IsValidPort == 0)
	{
		printf("is not valid port\n");

		return 0;
	}

	if ((buf_str != UARTPort.getPort())||(BaudRateI != UARTPort.getBaudrate()))
	{
		if (UARTPort.isOpen())
		{
			UARTPort.close();
		}
	}

	if (IsEnable && (!UARTPort.isOpen()))
	{
		UARTPort.setPort(buf);
		UARTPort.setBaudrate(BaudRateI);
		UARTPort.setTimeout(serial::Timeout::simpleTimeout(1000));
		UARTPort.open();
		
	}


	if (IsEnable == 0)
	{
		if (UARTPort.isOpen())
		{
			UARTPort.close();
		}
	}

	return 1;
}


/*
* Class:     hust_whno_SMLM_QC_STORM_Plug
* Method:    lm_SetFeedbackDevice
* Signature: (III)V
*/
JNIEXPORT void JNICALL Java_hust_whno_SMLM_QC_1STORM_1Plug_lm_1SetFeedbackDevice
(JNIEnv *env, jclass obj, jint ControlParaId, jint UARTID, jint DataRateI, jint IsEnableB)
{
#if DLL_CALL_DEBUG
	printf("lm_SetFeedbackDevice called\n");
#endif // DLL_CALL_DEBUG

	int suc = 0;

	switch (ControlParaId)
	{
	case 0:
		// activation laser control
		// z drift control
		// share the same UART

		suc = SetUARTPort(FeedbackCmdUART, UARTID, DataRateI, IsEnableB);

		if (suc)
		{
			OnlineFeedbackAlive = true;
			StartSendCmdUARTRunning();
		}

		break;

	case 1:

		break;

	}
}



/*
* Class:     hust_whno_SMLM_QC_STORM_Plug
* Method:    lm_SetFeedbackEnable
* Signature: (III)V
*/
JNIEXPORT void JNICALL Java_hust_whno_SMLM_QC_1STORM_1Plug_lm_1SetFeedbackEnable
(JNIEnv *env, jclass obj, jint ControlParaId, jint IsEnableB)
{
//	printf("feedback en:%d %d\n", ID, IsEnableB);

	switch (ControlParaId)
	{
	case 0:
		// activation laser control
		LocDensity_CtlData.SetFeedbackEn(IsEnableB);
		break;

	case 1:

		break;
	}

}





/*
* Class:     hust_whno_SMLM_QC_STORM_Plug
* Method:    lm_ResetFeedback
* Signature: ()V
*/
JNIEXPORT void JNICALL Java_hust_whno_SMLM_QC_1STORM_1Plug_lm_1ResetFeedback
(JNIEnv *, jclass)
{
	ResetFeedbackCtl();
}


/*
* Class:     hust_whno_SMLM_QC_STORM_Plug
* Method:    lm_SetFeedbackManualTarget
* Signature: (IIF)V
*/
JNIEXPORT void JNICALL Java_hust_whno_SMLM_QC_1STORM_1Plug_lm_1SetFeedbackManualTarget
(JNIEnv *, jclass, jint ControlParaId, jfloat Value, jint Enable)
{
//	printf("set ManualTarget:%d %f %d\n", id, Value, Enable);

	switch (ControlParaId)
	{
	case 0:
		// density

		if (Enable)
		{
			LocDensity_CtlData.SetManualTargetValue(Value);
		}
		else
		{
			LocDensity_CtlData.SetCtlMode(0);
			LocDensity_CtlData.ResetAutoTarget();
		}

		break;
	case 1:


		break;


	default:
		break;
	}
}

/*
* Class:     hust_whno_SMLM_QC_STORM_Plug
* Method:    lm_SetFeedbackPIDParameters
* Signature: (IIF)V
*/
JNIEXPORT void JNICALL Java_hust_whno_SMLM_QC_1STORM_1Plug_lm_1SetFeedbackPIDParameters
(JNIEnv *, jclass, jint ControlParaId, jfloat ProportionF, jfloat IntegralF, jfloat DerivativeF)
{
//	printf("set  pid:%d %f %f \n", id, ProportionF, IntegralF);


	switch (ControlParaId)
	{
	case 0:
		// density

		LocDensity_PID.SetPIDPara(ProportionF, IntegralF, 0);

		break;
	case 1:


		break;

	default:
		break;
	}
}



/*
* Class:     hust_whno_SMLM_QC_STORM_Plug
* Method:    lm_GetFirstUARTId
* Signature: ()I
*/
JNIEXPORT jint JNICALL Java_hust_whno_SMLM_QC_1STORM_1Plug_lm_1GetFirstUARTId
(JNIEnv *env, jclass obj)
{

	int id = GetFirstPortId();
	return id;
}

/*
* Class:     hust_whno_SMLM_QC_STORM_Plug
* Method:    lm_ZDepthSMMove
* Signature: (I)V
*/
JNIEXPORT void JNICALL Java_hust_whno_SMLM_QC_1STORM_1Plug_lm_1ZDepthSMMove
(JNIEnv *env, jclass obj, jint MoveSteps)
{
#if DLL_CALL_DEBUG
	printf("lm_ZDepthSMMove called\n");
#endif // DLL_CALL_DEBUG

	// send to MCU
	char *UARTBuf = new char[64];
	string cmdString;

	// sent to the MCU to drive the step motor to adjust laser power by rotating the ND filter	
	sprintf(UARTBuf, "@z%d#", MoveSteps);
	cmdString = UARTBuf;

	UARTCmdQueue.push(cmdString);
	delete[]UARTBuf;
	
}


/*
* Class:     hust_whno_SMLM_QC_STORM_Plug
* Method:    lm_FeedbackCtlTest
* Signature: (II)V
*/
JNIEXPORT void JNICALL Java_hust_whno_SMLM_QC_1STORM_1Plug_lm_1FeedbackCtlTest
(JNIEnv *env, jclass obj, jint ControlParaId, jint SetState)
{
	switch (ControlParaId)
	{
	case 0:
		// density
		if (SetState)
		{
			LocDensityTest_Set();
		}
		else
		{
			LocDensityTest_Reset();
		}

		break;
	case 1:


		break;

	default:
		break;
	}


}

// multi ROI acquisition, transtlation stage control

/*
* Class:     hust_whno_SMLM_QC_STORM_Plug
* Method:    lm_SetTranslationStage
* Signature: (III)V
*/
JNIEXPORT void JNICALL Java_hust_whno_SMLM_QC_1STORM_1Plug_lm_1SetTranslationStage
(JNIEnv *env, jclass obj, jint UARTID, jint DataRateI, jint IsEnableB)
{

	SetUARTPort(TranslationStageUART, UARTID, DataRateI, IsEnableB);

}


/*
* Class:     hust_whno_SMLM_QC_STORM_Plug
* Method:    lm_TranslationStageMove
* Signature: (III)V
*/
JNIEXPORT void JNICALL Java_hust_whno_SMLM_QC_1STORM_1Plug_lm_1TranslationStageMove
(JNIEnv *env, jclass obj, jint XSteps, jint YSteps, jint ZSteps)
{
	// send to MCU
	char *UARTBuf = new char[64];
	string cmdString;

	int TimeDelay;
	
	// x move
	if (abs(XSteps) > 0)
	{
		if (XSteps > 0)sprintf(UARTBuf, "#+X %d#", abs(XSteps));
		else sprintf(UARTBuf, "#-X %d#", abs(XSteps));

		cmdString = UARTBuf;

		// send cmd
		if (TranslationStageUART.isOpen())
		{
			TranslationStageUART.write(cmdString);
		}

		printf("stage: %s\n", cmdString);

		TimeDelay = 5 * abs(XSteps); //
		
		Sleep(TimeDelay);
	}

	// y move
	if (abs(YSteps) > 0)
	{
		if (YSteps > 0)sprintf(UARTBuf, "#+Y %d#", abs(YSteps));
		else sprintf(UARTBuf, "#-Y %d#", abs(YSteps));

		cmdString = UARTBuf;

		// send cmd
		if (TranslationStageUART.isOpen())
		{
			TranslationStageUART.write(cmdString);
		}

		printf("stage: %s\n", cmdString);

		TimeDelay = 5 * abs(YSteps); // 
		
		Sleep(TimeDelay);

	}
}


/*
* Class:     hust_whno_SMLM_QC_STORM_Plug
* Method:    lm_SetMultiFOVAcqParameters
* Signature: (III)V
*/
JNIEXPORT void JNICALL Java_hust_whno_SMLM_QC_1STORM_1Plug_lm_1SetMultiFOVAcqParameters
(JNIEnv *env, jclass obj, jfloat iFOVOverlapPercent)
{
	FOVOverlapPercent = iFOVOverlapPercent;

}

/*
* Class:     hust_whno_SMLM_QC_STORM_Plug
* Method:    lm_LocBatchedImg
* Signature: ([SI)[F
*/
JNIEXPORT jfloatArray JNICALL Java_hust_whno_SMLM_QC_1STORM_1Plug_lm_1LocBatchedImg
(JNIEnv *env, jclass obj, jshortArray jImgArry, jint BatchedImgNum)
{
#if DLL_CALL_DEBUG
	printf("lm_LocBatchedImg called\n");
#endif // DLL_CALL_DEBUG

	// get raw images	// jni environment

	int len = (*env).GetArrayLength(jImgArry);

	jboolean iscopy = false;
	jshort * elems = (jshort *)env->GetPrimitiveArrayCritical(jImgArry, &iscopy);
	//	jshort * elems = (*env).GetShortArrayElements(jImgArry, NULL);

	// only work at the first time
	InitAllLocResource();


	int BatchedImgSize = BatchedImgNum*LocPara_Global.ImageWidth*LocPara_Global.ImageHigh;

//	printf("rec dat:%d %d %d\n", len, BatchedImgSize, BatchedImgNum);



	memcpy(ZDriftCtl.GetRawImageMem(), elems, BatchedImgSize*sizeof(short));

	//	env->ReleaseShortArrayElements(jImgArry, elems, 0);
	env->ReleasePrimitiveArrayCritical(jImgArry, elems, JNI_ABORT);



	int FluoNum = ZDriftCtl.MLELocalization(NULL, BatchedImgNum);



	// return results 
	int ReturnInfSize = 8;
	float *oInf = new float[ReturnInfSize];


	float Mean_PSFWidth = ZDriftCtl.FluoStatData_.MeanPSFWidth;
	float Mean_SNR = ZDriftCtl.FluoStatData_.MeanPeakSNR;
	float Mean_PSFWidth_Ctl = ZDriftCtl.FluoStatData_.MeanPSFWidth_Ctl;
	float Mean_LocDensity2D = ZDriftCtl.FluoStatData_.MeanLocDensity2D;


	oInf[0] = Mean_PSFWidth;
	oInf[1] = -Mean_SNR;

	oInf[2] = Mean_PSFWidth - Mean_SNR;
	oInf[3] = Mean_PSFWidth*(-Mean_SNR);
	oInf[4] = FluoNum;

	// set results contents
	jfloatArray result = env->NewFloatArray(ReturnInfSize);

	env->SetFloatArrayRegion(result, 0, ReturnInfSize, oInf);

	return result;
}


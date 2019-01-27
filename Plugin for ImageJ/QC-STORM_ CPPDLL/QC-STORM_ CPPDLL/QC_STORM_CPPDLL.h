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

/* DO NOT EDIT THIS FILE - it is machine generated */
#include <jni.h>
/* Header for class QC_STORM_ */

#ifndef _Included_QC_STORM_
#define _Included_QC_STORM_
#ifdef __cplusplus
extern "C" {
#endif
/*
* Class:     QC_STORM_
* Method:    lm_SetImagePara
* Signature: (IIIII[C)V
*/
JNIEXPORT void JNICALL Java_QC_1STORM_1_lm_1SetImagePara
	(JNIEnv *, jclass, jint, jint, jint, jint, jint, jcharArray);

/*
* Class:     QC_STORM_
* Method:    lm_SetLocPara
* Signature: (FFFIIIFFF)V
*/
JNIEXPORT void JNICALL Java_QC_1STORM_1_lm_1SetLocPara
	(JNIEnv *, jclass, jfloat, jfloat, jfloat, jint, jint, jint, jint, jfloat, jfloat, jfloat, jfloat);

/*
* Class:     QC_STORM_
* Method:    lm_SetLocPara3D
* Signature: (FFFFFFFFI)V
*/
JNIEXPORT void JNICALL Java_QC_1STORM_1_lm_1SetLocPara3D
	(JNIEnv *, jclass, jfloat, jfloat, jfloat, jfloat, jfloat, jfloat, jfloat, jfloat);

/*
* Class:     QC_STORM_
* Method:    lm_SetStatInfSelection
* Signature: (I)V
*/
JNIEXPORT void JNICALL Java_QC_1STORM_1_lm_1SetStatInfSelection
(JNIEnv *, jclass, jint, jint, jint);

/*
* Class:     QC_STORM_
* Method:    lm_StartLocThread
* Signature: ()V
*/
JNIEXPORT void JNICALL Java_QC_1STORM_1_lm_1StartLocThread
	(JNIEnv *, jclass);

/*
* Class:     QC_STORM_
* Method:    lm_StopLocThread
* Signature: ()V
*/
JNIEXPORT void JNICALL Java_QC_1STORM_1_lm_1StopLocThread
	(JNIEnv *, jclass);

/*
* Class:     QC_STORM_
* Method:    lm_GetMaxDispVal
* Signature: ()I
*/
JNIEXPORT jint JNICALL Java_QC_1STORM_1_lm_1GetMaxDispVal
	(JNIEnv *, jclass);

/*
* Class:     QC_STORM_
* Method:    lm_FeedImageData
* Signature: ([SI)V
*/
JNIEXPORT void JNICALL Java_QC_1STORM_1_lm_1FeedImageData
	(JNIEnv *, jclass, jshortArray, jint);

/*
* Class:     QC_STORM_
* Method:    lm_GetSMLMImage
* Signature: ()[F
*/
JNIEXPORT jfloatArray JNICALL Java_QC_1STORM_1_lm_1GetSMLMImage
	(JNIEnv *, jclass);

/*
* Class:     QC_STORM_
* Method:    lm_GetSMLMImage3D
* Signature: ()[I
*/
JNIEXPORT jintArray JNICALL Java_QC_1STORM_1_lm_1GetSMLMImage3D
	(JNIEnv *, jclass);

/*
* Class:     QC_STORM_
* Method:    lm_SetSpatialResolutionInf
* Signature: (IF)V
*/
JNIEXPORT void JNICALL Java_QC_1STORM_1_lm_1SetSpatialResolutionInf
	(JNIEnv *, jclass, jint, jint, jfloat, jfloat);


/*
* Class:     QC_STORM_
* Method:    lm_GetStatInfImageSize
* Signature: (I)[I
*/
JNIEXPORT jintArray JNICALL Java_QC_1STORM_1_lm_1GetStatInfImageSize
(JNIEnv *env, jclass obj);



/*
* Class:     QC_STORM_
* Method:    lm_GetStatInfImage
* Signature: (I)[I
*/
JNIEXPORT jintArray JNICALL Java_QC_1STORM_1_lm_1GetStatInfImage
	(JNIEnv *, jclass, jint);

/*
* Class:     QC_STORM_
* Method:    lm_IsLocFinish
* Signature: ()I
*/
JNIEXPORT jint JNICALL Java_QC_1STORM_1_lm_1IsLocFinish
	(JNIEnv *, jclass);

/*
* Class:     QC_STORM_
* Method:    lm_SetRerendImagePara
* Signature: ([CII)V
*/
JNIEXPORT void JNICALL Java_QC_1STORM_1_lm_1SetRerendImagePara
	(JNIEnv *, jclass, jcharArray, jint, jint);

/*
* Class:     QC_STORM_
* Method:    lm_StartRerend
* Signature: ()V
*/
JNIEXPORT void JNICALL Java_QC_1STORM_1_lm_1StartRerend
	(JNIEnv *, jclass);

/*
* Class:     QC_STORM_
* Method:    lm_GetRerendImageInf
* Signature: ()[I
*/
JNIEXPORT jintArray JNICALL Java_QC_1STORM_1_lm_1GetRerendImageInf
	(JNIEnv *, jclass);

/*
* Class:     QC_STORM_
* Method:    lm_ReleaseRerendResource
* Signature: ()V
*/
JNIEXPORT void JNICALL Java_QC_1STORM_1_lm_1ReleaseRerendResource
	(JNIEnv *, jclass);

#ifdef __cplusplus
}
#endif
#endif
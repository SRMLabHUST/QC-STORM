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


package hust.whno.SMLM;



import org.micromanager.data.ProcessorConfigurator;
import org.micromanager.data.ProcessorPlugin;
import org.micromanager.data.ProcessorFactory;

import org.micromanager.PropertyMap;
import org.micromanager.Studio;

import org.scijava.plugin.Plugin;
import org.scijava.plugin.SciJavaPlugin;


@Plugin(type = ProcessorPlugin.class)
public class QC_STORM_Plug implements ProcessorPlugin, SciJavaPlugin{
    private Studio studio_;
    
    public QC_STORM_Configurator myCinfigConfigurator;
    
    public static int ProcessorID = 0;
    
    // wrapper of CPP & CUDA library
    public native static void lm_SetProcessorID(int id);
    public native static int lm_GetProcessorID();
    public native static void lm_SetSavePath(char PathCs[]);
    // for both 2d and 3d
    public native static void lm_SetLocPara(float KAdc, float Offset, float QE, int ROISizeI, int LocTypeI, int MultiEmitterFitEn, int ConsecutiveFitEnI, float ConsecFilterRadiusF, float RawPixelSizeF, float RenderPixelZoomF, float SNR_th);
    public native static void lm_SetLocPara3D(float MinZDepthF, float MaxZDepthF, float ZDepthCorrFactor, float p4, float p3, float p2, float p1, float p0);
    
    public native static void lm_SetImagePara(int ImageWidth, int ImageHigh, int SRImageWidth, int SRImageHigh);
    public native static void lm_SetAcquisitionPara(char CreateTimeIdx[]);
    
    public native static void lm_SetStatInfSelection(int DispSel, int OnTimeEn, int SpatialResolutionEn);
    
    // send image to be processed by GPU
    public native static void lm_FeedImageData(short ImgDataS[], int FrameNumI);

    public native static void lm_StartLocThread();
    public native static void lm_StopLocThread();
    public native static void lm_ReleaseLocResource();
    public native static int lm_IsLocFinish();
            
    public native static int lm_GetMaxDispVal();     
    public native static int lm_GetWaitImageNum();     

    public native static float [] lm_GetSMLMImage();
    public native static int [] lm_GetSMLMImage3D();
   
    public native static void lm_SetSpatialResolutionInf(int FramePerGroup, int IsHollowTube, float StructureSize, float RSCResolutionTh);
    public native static float lm_GetCurSpatialResolution();
    public native static float lm_GetMeanLocPrec();
    
    // feedback control
    public native static void lm_SetFeedbackDevice(int ControlParaId, int UARTID, int DataRateI, int IsEnableB);
    public native static void lm_SetFeedbackEnable(int ControlParaId, int IsEnableB);
    
    public native static void lm_SetFeedbackManualTarget(int ControlParaId, float ValueF, int EnableI);
    public native static void lm_SetFeedbackPIDParameters(int ControlParaId, float ProportionF, float IntegralF, float DerivativeF);
    public native static void lm_ResetFeedback();
    public native static int lm_GetFirstUARTId();
    
    public native static void lm_ZDepthSMMove(int MoveSteps);
    
    public native static void lm_FeedbackCtlTest(int ControlParaId, int SetState);
    
    // multi FOV acquisition, transtlation stage control
    public native static void lm_SetTranslationStage(int UARTID, int DataRateI, int IsEnableB);
    public native static void lm_TranslationStageMove(int XSteps, int YSteps, int ZSteps);
    public native static void lm_SetMultiFOVAcqParameters(float FOVOverlapPercent);
    
    // for z drift control, independent localization of several images
    public native static float [] lm_LocBatchedImg(short ImgDataS[], int BatchedImgNum);
    
    
    public static int SendProcessorID()
    {
        int CurId = ProcessorID;
        ProcessorID++;
        lm_SetProcessorID(CurId);
        
        return CurId;
    }
    public static void CloseAndWaitCurrentLocThread()
    {
        lm_StopLocThread();
        
        while(lm_IsLocFinish()==0)
        {
            // wait GPU localization finsih
        }
        
    }
    
    @Override
    public ProcessorConfigurator createConfigurator(PropertyMap pm) {
          
        System.loadLibrary ("QC-STORM CPPDLL");

        myCinfigConfigurator = new QC_STORM_Configurator(studio_, pm);
        
        return myCinfigConfigurator;
    }

    @Override
    public ProcessorFactory createFactory(PropertyMap pm) {
        return new QC_STORM_Factory (studio_, myCinfigConfigurator);
    }

    @Override
    public void setContext(Studio studio) {
        studio_ = studio;
                
    }
    

    @Override
    public String getName() {
        return "QC-STORM";
    }

    @Override
    public String getHelpText() {
        return "online processing SMLM";
    }

    @Override
    public String getVersion() {
        return "1.0";
    }

    @Override
    public String getCopyright() {
        return "hust.whno";
    }
    
}

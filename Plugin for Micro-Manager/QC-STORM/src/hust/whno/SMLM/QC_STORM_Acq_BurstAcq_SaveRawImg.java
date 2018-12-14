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

import java.util.ArrayList;
import org.micromanager.ScriptController;
import org.micromanager.Studio;
import mmcorej.CMMCore;
import mmcorej.TaggedImage;
import org.micromanager.data.Datastore;
import org.micromanager.data.Image;
import org.micromanager.data.Pipeline;
import org.micromanager.data.ProcessorFactory;
import org.micromanager.display.DisplaySettings;
import org.micromanager.display.DisplayWindow;



/**
 *
 * another method to display and process is use datastore and pipeline
 * however, cpu memory useage continously increasing and result in error in multi ROI acquisition
 */
public class QC_STORM_Acq_BurstAcq_SaveRawImg {
    private Studio studio_;
    CMMCore mmc;
    ScriptController gui;
   
    public QC_STORM_Configurator MyConfigurator; // reference to configurator
    
    
    public BurstLiveAcqThread BurstAcqThread=null;
    
    int BurstAcqFrameNum;
    float SpatialResolutionTh;
        
    Datastore MyDatastore;
    QC_STORM_Factory MyFactory;
    java.util.List<ProcessorFactory> MyFactoryList;
    DisplayWindow Mydisp;
    Pipeline MyPipeline;
    

      
    QC_STORM_Acq_BurstAcq_SaveRawImg(Studio studio, QC_STORM_Configurator iConfigurator) throws Exception
    {
        studio_=studio;
        mmc=studio_.getCMMCore();
        gui=studio_.scripter();
        
        MyConfigurator=iConfigurator;

        BurstAcqFrameNum=MyConfigurator.GetBurstAcqFrameNum();
        SpatialResolutionTh=MyConfigurator.GetSpatialResolutionTh();


        // create data store and processing pipeline
        MyDatastore = studio.data().createRAMDatastore();
        MyFactoryList = new ArrayList<ProcessorFactory>();
        
        MyFactory = new QC_STORM_Factory(studio, MyConfigurator);
        
        MyFactoryList.add(MyFactory);
        MyPipeline=studio_.data().createPipeline(MyFactoryList, MyDatastore, true);

        Mydisp=studio_.displays().createDisplay(MyDatastore);
        
        DisplaySettings MyDispSetting=Mydisp.getDisplaySettings();
        Mydisp.setDisplaySettings(MyDispSetting.copy().channelColorMode(DisplaySettings.ColorMode.GRAYSCALE).build());
           
        studio_.displays().manage(MyDatastore);

//		gui.message("into QC_STORM_BurstLiveProc");

    }

       
    public void StartBurstAcq()
    {			
        // start loc thread
        BurstAcqThread = new BurstLiveAcqThread();
        BurstAcqThread.start();
        
    }
    public void WaitBurstAcqFinish()
    {
        if(BurstAcqThread!=null)
        {
            try {
                BurstAcqThread.join();
                
            } catch (Exception ex) {
                gui.message("except BurstAcqThread join");
            }
        }
    }
    
    public class BurstLiveAcqThread extends Thread
    {

        @Override
        public void run()
        {            
            
            float CurSpatialResolution;
            
            MyConfigurator.UpdateBurstAcqFrameNum(0);
            
            // set localization parameters 
            MyConfigurator.SendLocalizationPara();

			//
			try {

				int curFrame = 0;
                
                int CorrImgNumI = Math.max(MyConfigurator.GetZDriftCorrFrameNum(), 100);
                int GroupNumI = (BurstAcqFrameNum + CorrImgNumI - 1) /CorrImgNumI;
                
                boolean IsBreak=false;
                
                for(int gcnt=0; (gcnt< GroupNumI) && MyConfigurator.CurBurstLiveActive; gcnt++)
                {
                    int nrFrames = Math.min(CorrImgNumI, BurstAcqFrameNum - gcnt*CorrImgNumI);
                   
                    // Start burst acq
                    
                    double exposureMs = mmc.getExposure();
                    mmc.startSequenceAcquisition(nrFrames, 0, true);

                    
                    while (MyConfigurator.CurBurstLiveActive &&((mmc.getRemainingImageCount() > 0)|| mmc.isSequenceRunning(mmc.getCameraDevice()))) 
                    {
                        if(IsBreak)break;
                        
                        // reference scriper bust acquisition
                        if (mmc.getRemainingImageCount() > 0) 
                        {
                            TaggedImage taggedImg = mmc.popNextTaggedImage();

                            Image image1 = studio_.data().convertTaggedImage(taggedImg, studio_.data().getCoordsBuilder().time(curFrame).build(), null);
                            
                            MyDatastore.putImage(image1);
                            MyPipeline.insertImage(image1);
                            
                            
                            // update acquired frame
                            curFrame++;
                            
                            
                            // don't update GUI to frequently
                            if(curFrame % 10 == 0)
                            {
                                CurSpatialResolution = QC_STORM_Plug.lm_GetCurSpatialResolution();
                                float CurMeanLocPrec = QC_STORM_Plug.lm_GetMeanLocPrec();
                                
                                MyConfigurator.SetBufferedImgNum(QC_STORM_Plug.lm_GetWaitImageNum(), mmc.getRemainingImageCount(), (int)CurSpatialResolution); // show buffered image num in cuda fifo
                                MyConfigurator.UpdateBurstAcqFrameNum(curFrame);  
                                
                                // normal acquisition
                                if(SpatialResolutionTh > 0)
                                {
                                    // when use 3 oversampling, can achieve same locprecision in x and 
                                    if((CurSpatialResolution <= SpatialResolutionTh)|| (CurSpatialResolution < 1.1f*2.35f*CurMeanLocPrec)) // 
                                    {
                                        IsBreak = true;
                                    }
                                    
                                    // after many frames, the spatial resolution is still high, thus could be blank images
                                    if((curFrame >= 500) && (CurSpatialResolution > 500.0f))
                                    {
                                        IsBreak = true;
                                    }
                                }
                            }
                        }
                        else {
                            // Wait for another image to arrive.
                            mmc.sleep(Math.min(.4 * exposureMs, 10));
                        }
                    }
                    mmc.stopSequenceAcquisition();
                    
                    
                    // z drift correct
                    if((gcnt < GroupNumI-1)  && MyConfigurator.CurBurstLiveActive && MyConfigurator.IsZDriftCtlEn())
                    {
                        // correct z drift by move to the place where the PSF width is minimum
                        QC_STORM_ZDriftCorrection ZdriftCorrThread = new QC_STORM_ZDriftCorrection(studio_, MyConfigurator, 0, 1, 2, 1); 
                        ZdriftCorrThread.start();
                        ZdriftCorrThread.join();
                    }
                }
                
                MyFactory.MyProcessor.cleanup(null);
                
			}
            catch (Exception ex) 
            {
                gui.message("get BurstLiveAcqThread exception");
			}
            // finish
			MyConfigurator.ResetBurstLiveBtn();

		}
    }
}

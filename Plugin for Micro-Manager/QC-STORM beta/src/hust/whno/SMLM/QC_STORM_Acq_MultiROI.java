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

import mmcorej.CMMCore;
import org.micromanager.ScriptController;
import org.micromanager.Studio;



public class QC_STORM_Acq_MultiROI {
    
    private Studio studio_;
    CMMCore mmc;
    ScriptController gui;
   
    public QC_STORM_Configurator MyConfigurator; // reference to configurator

    
    // multi ROI acquisition
    private MultiROIAcqThread MultiAcqThread;

    public int ImageWidthI,ImageHighI;    
    
    
    QC_STORM_Acq_MultiROI(Studio studio, QC_STORM_Configurator iConfigurator)
    {
        studio_=studio;
        mmc=studio_.getCMMCore();
        gui=studio_.scripter();
        
        MyConfigurator=iConfigurator;
        
        ImageWidthI = (int) mmc.getImageWidth();
        ImageHighI = (int) mmc.getImageHeight();              
        
    
//		gui.message("into QC_STORM_BurstLiveProc");

    }

    public void StartMultiROIAcq()
    {
                
        // start loc thread
        MultiAcqThread = new MultiROIAcqThread();
        MultiAcqThread.start();
   
    }    

    public class MultiROIAcqThread extends Thread
    {
        
        @Override
        public void run()
        {
            try {
                int ROINum_X=MyConfigurator.MultiROINum_X;
                int ROINum_Y=MyConfigurator.MultiROINum_Y;
                int ROISteps_X=MyConfigurator.ROIMoveSteps_X;
                int ROISteps_Y=MyConfigurator.ROIMoveSteps_Y;   
                
                gui.message("into MultiROIAcqThread");

                for(int ycnt=0; ycnt < ROINum_Y; ycnt++)
                {
                    for(int xcnt=0; xcnt < ROINum_X; xcnt++)
                    {
                       
                        // if this is selected, don't perform acquisition initialization
                        if(!MyConfigurator.IsNoAcqIni())
                        {
                            // initial before acquisition
                            NewAcqInitial acqInitial = new NewAcqInitial(ycnt*ROINum_X+xcnt, ROINum_X*ROINum_Y);
                            acqInitial.start();
                            acqInitial.join();
                        }
                        // start acquisition
                        if(!MyConfigurator.CurMultiROIAcqActive)break;
                        
                        MyConfigurator.CurBurstLiveActive = true;
                        QC_STORM_Acq_BurstAcq BurstLiveProc = new QC_STORM_Acq_BurstAcq(studio_, MyConfigurator);
                        
                        BurstLiveProc.StartBurstAcq();
                        BurstLiveProc.WaitBurstAcqFinish();
                       
                        
                        // move ROI first and then save raw images to save time, such as wait density down
                        if(xcnt < ROINum_X-1)
                        {
                            // x move a ROI
                            QC_STORM_Plug.lm_TranslationStageMove(ROISteps_X, 0, 0);   
                        }
                        // a line is finished
                        if(xcnt == ROINum_X-1)
                        {
                            if(ycnt < ROINum_Y-1)
                            {
                                for(int i = 0; i < ROINum_X-1; i++)
                                {
                                    // x return to the front
                                    QC_STORM_Plug.lm_TranslationStageMove(-ROISteps_X, 0, 0);
                                }

                                // y move a ROI
                                QC_STORM_Plug.lm_TranslationStageMove(0, ROISteps_Y, 0);       
                            }
                        }
                    }
                }  
            } catch (Exception ex) {
                gui.message("get exception");
            }
            
            MyConfigurator.ResetMultiROIAcqBtn();                  
        }
    }
    
    public class NewAcqInitial extends Thread
    {
        int CurId;
        int TotalROINum;
        
        public NewAcqInitial(int CurId, int TotalROINum)
        {
            this.CurId = CurId;
            this.TotalROINum = TotalROINum;
        }

        @Override
        public void run()
        {
            // initial acquisition befor acquisition of each ROI

            QC_STORM_ZDriftCorrection ZDriftCorrThread;
            WaitLocDensityDown WaitLocDensityThread;
            
            try {
                    QC_STORM_Plug.lm_ResetFeedback();

                    // make sure z drift is corrected and the density is ok  
                    // correct z drift befor acquistion
                    
                    ZDriftCorrThread = new QC_STORM_ZDriftCorrection(studio_, MyConfigurator, 4, 4, 3, 3);
                    ZDriftCorrThread.start();
                    ZDriftCorrThread.join();

                    
                    ZDriftCorrThread = new QC_STORM_ZDriftCorrection(studio_, MyConfigurator, MyConfigurator.ZCorrMode(), 1, 3, 2);
                    ZDriftCorrThread.start();
                    ZDriftCorrThread.join();

                    float MaxTolerableDensity = MyConfigurator.GetMaxTolerableLocDensity();
                    
                    if(MyConfigurator.IsRandomDensity())
                    {
                        MaxTolerableDensity = 0.35f - (float)CurId / (float)TotalROINum * 0.25f;
 //                       MaxTolerableDensity = 0.15f + (float)Math.random() * 0.2f;
                    }
                            
                    // wait localization density down to software's ability
                    WaitLocDensityThread = new WaitLocDensityDown(true, MaxTolerableDensity);
                    WaitLocDensityThread.start();
                    WaitLocDensityThread.join();

                } catch (InterruptedException ex) {              
            }
        }
    }
    
    public class WaitLocDensityDown extends Thread
    {
        boolean WaitDensityEn;
        float DensityTh;
        
        WaitLocDensityDown(boolean iWaitDensityEn, float iDensityTh)
        {
            WaitDensityEn=iWaitDensityEn;
            DensityTh=iDensityTh;
        }
        
        @Override
        public void run()
        {
            while(WaitDensityEn)
            {
                if(!MyConfigurator.CurMultiROIAcqActive)break;
                
                // also correct PSF when wait density down
                try {
                    
                    QC_STORM_ZDriftCorrection ZdriftCorrThread = new QC_STORM_ZDriftCorrection(studio_, MyConfigurator, MyConfigurator.ZCorrMode(), 1, 2, 1); 
                    
                    // get localization density
                    int BatchedImgNum = ZdriftCorrThread.GetBatchedImgNum();
                    short [] BatchedImgDat = ZdriftCorrThread.GetBatchedImgDat(BatchedImgNum);
                    
                    float [] LocResults = QC_STORM_Plug.lm_LocBatchedImg(BatchedImgDat, BatchedImgNum);
                    
                    float CurLocDensity = LocResults[3];

                    
                    if(CurLocDensity <= DensityTh)
                    {
                        break;
                    }
                    else
                    {
                        // z drift correction in the wait process
                        ZdriftCorrThread.start();
                        ZdriftCorrThread.join();
                    }
                    
                } catch (InterruptedException ex) {
                    
                }
            }
        }
    }
}

  /*  
    Datastore MyDatastore;
    DisplayWindow Mydisp;
     
        // create data store and processing pipeline
        MyDatastore = studio.data().createRAMDatastore();

        Mydisp=studio_.displays().createDisplay(MyDatastore);
        studio_.displays().manage(MyDatastore);
  
    */
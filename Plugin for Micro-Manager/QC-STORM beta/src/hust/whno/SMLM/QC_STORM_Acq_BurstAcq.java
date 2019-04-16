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

import ij.ImagePlus;
import ij.gui.ImageWindow;
import ij.process.ShortProcessor;
import mmcorej.CMMCore;
import mmcorej.TaggedImage;
import org.micromanager.ScriptController;
import org.micromanager.Studio;
import org.micromanager.data.Datastore;
import org.micromanager.data.Image;
import org.micromanager.display.DisplayWindow;



public class QC_STORM_Acq_BurstAcq {

    private Studio studio_;
    CMMCore mmc;
    ScriptController gui;
   
    public QC_STORM_Configurator MyConfigurator; // reference to configurator
    
    QC_STORM_Processor MyDataProcessor;
    
    public BurstLiveAcqThread BurstAcqThread=null;
    
    int BurstAcqFrameNum;
    float SpatialResolutionTh;
        
    public int ImageWidthI,ImageHighI;    

    ShortProcessor RawImageProcessor;
    public ImagePlus RawImagePlus;
    public ImageWindow RawImageDispWindow;
    
    boolean IsResolutionAchieved = false;
    
    // save raw image
    boolean IsSaveRawImage = false;
    Datastore MyDatastore;
    DisplayWindow Mydisp;
    String ImageSavePath;
    
    
    
    QC_STORM_Acq_BurstAcq(Studio studio, QC_STORM_Configurator iConfigurator, String iNamePostFix) throws Exception
    {
        studio_ = studio;
        mmc = studio_.getCMMCore();
        gui = studio_.scripter();
        
        MyConfigurator = iConfigurator;
        
        
        
        // send localization para and image data for processing
        MyDataProcessor = new QC_STORM_Processor(studio, iConfigurator, iNamePostFix);

        BurstAcqFrameNum = MyConfigurator.GetBurstAcqFrameNum();
        SpatialResolutionTh = MyConfigurator.GetSpatialResolutionTh();


        ImageWidthI = (int)mmc.getImageWidth();
        ImageHighI = (int)mmc.getImageHeight();


        RawImageProcessor = new ShortProcessor(ImageWidthI, ImageHighI);
        RawImagePlus = new ImagePlus("raw image", RawImageProcessor);
        RawImageDispWindow = new ImageWindow(RawImagePlus); 
        
        
        // save raw images
        IsSaveRawImage = MyConfigurator.IsRawImageSave();
        
        if(IsSaveRawImage)
        {
            String ImageNamePreFix = String.format("RawImage_%dD_", MyConfigurator.LocPara.LocType + 2);
            ImageSavePath = MyConfigurator.GetResultsSavePath() + ImageNamePreFix + MyDataProcessor.CreateTimeIdxStr;
            
            MyDatastore = studio_.data().createMultipageTIFFDatastore(ImageSavePath, true, true);
            
            
            // don't managed by micro-manager, or the memory will always increase and not enough
//            studio_.displays().manage(MyDatastore);
//           Mydisp = studio_.displays().createDisplay(MyDatastore);
            
        }

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
                    
                    int DispGap = Math.max((int)(1000/5/exposureMs), 1);// display 5 images per second
                    
                    
                    while (MyConfigurator.CurBurstLiveActive &&((mmc.getRemainingImageCount() > 0)|| mmc.isSequenceRunning(mmc.getCameraDevice()))) 
                    {
                        if(IsBreak)break;
                        
                        // reference scriper bust acquisition
                        if (mmc.getRemainingImageCount() > 0) 
                        {
                            TaggedImage taggedImg = mmc.popNextTaggedImage();

                            Image image1 = studio_.data().convertTaggedImage(taggedImg, studio_.data().getCoordsBuilder().time(curFrame).build(), null);
                            curFrame++;  
                            
                            // manually mimic pipeline, process images by GPU
                            MyDataProcessor.ProcessAImage(image1);
                
                            
                            if(IsSaveRawImage)
                            {
                                MyDatastore.putImage(image1);
                            }

                            
                            // manually display raw images, don't display every image to improve efficiency
                            if((curFrame-1)%DispGap==0)
                            {
                                RawImageProcessor.setPixels((short[])image1.getRawPixels());

                                RawImagePlus.setProcessor(RawImageProcessor);
                                RawImagePlus.updateImage();
                                RawImagePlus.resetDisplayRange(); 

                                if(RawImageDispWindow.isClosed())
                                {
                                    RawImageDispWindow = new ImageWindow(RawImagePlus); 
                                }
                                RawImagePlus.draw();         
                            }

                            
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
                                        IsResolutionAchieved = true;
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
                            mmc.sleep(Math.max(.4 * exposureMs, 1));
                        }
                    }
                    mmc.stopSequenceAcquisition();
                    
                    if(IsResolutionAchieved)break;
                    
                    // z drift correct
                    if((gcnt < GroupNumI-1)  && MyConfigurator.CurBurstLiveActive && MyConfigurator.IsZDriftCtlEn())
                    {
                        // correct z drift by move to the place where the PSF width is minimum
                        QC_STORM_ZDriftCorrection ZdriftCorrThread = new QC_STORM_ZDriftCorrection(studio_, MyConfigurator, MyConfigurator.ZCorrMode(), 1, 2, 1); 
                        ZdriftCorrThread.start();
                        ZdriftCorrThread.join();
                    }
                }
                
                // save super resolution image
                MyDataProcessor.CleanProcessing();

                if(IsSaveRawImage)
                {
                    MyDatastore.freeze();
                    MyDatastore.setSavePath(ImageSavePath);
                }
            }
            
            catch (Exception ex) 
            {
                gui.message("get BurstLiveAcqThread exception");
			}
            
            // finish
			MyConfigurator.ResetBurstLiveBtn();

            MyConfigurator.ResetFeedbackCtl();
            
            RawImagePlus.close();
		}
    }
}


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

import static java.lang.Thread.sleep;
import mmcorej.CMMCore;
import mmcorej.TaggedImage;
import org.micromanager.ScriptController;
import org.micromanager.Studio;
import org.micromanager.data.Image;

/**
 *
 * z drift correction by move the focus knob to minimize the PSF width, may also consider molecule number and SNR
 */

public class QC_STORM_ZDriftCorrection  extends Thread
{
        private Studio MyStudio;
        private CMMCore MMCore;
        private ScriptController MMScript;
    
        QC_STORM_Configurator myConfigurator;
                
        int ImageWidth;
        int ImageHigh;
        
        int ZCorrMode;
        int StepZoom;
        // move step motor to correct z drift
        int SMMoveNum_half;
        int SMMoveNum;
        int IterationNum;
        
                
        QC_STORM_ZDriftCorrection(Studio iStudio, QC_STORM_Configurator iConfigurator, int iZCorrMode, int iStepZoom, int iSMMoveNum_half, int iIterationNum)
        {
        	MyStudio=iStudio;
        	MMCore=MyStudio.getCMMCore();
        	MMScript=MyStudio.scripter();
            
            myConfigurator = iConfigurator;
            
            ImageWidth = (int) MMCore.getImageWidth();
            ImageHigh = (int) MMCore.getImageHeight();
            
            ZCorrMode = iZCorrMode;
            StepZoom = iStepZoom;
        	//
        	SMMoveNum_half = iSMMoveNum_half;
        	SMMoveNum = SMMoveNum_half*2+1;
        	IterationNum = iIterationNum;
        }

        @Override
        public void run()
        {
            // set localization parameters 
            myConfigurator.SendLocalizationPara();
            myConfigurator.SetFeedbackDevicePort();
            
            
            int BatchedImgNum = GetBatchedImgNum(); 
            int SMMoveSteps = myConfigurator.GetZMoveSteps() ;
            
            SMMoveSteps = SMMoveSteps * StepZoom;


            float [] PSFVaryInf = new float[SMMoveNum];
            
            for(int icnt=0; icnt<IterationNum; icnt++)
            {
                // pre move
                ZDriftSMMoveAndWait(-SMMoveNum_half*SMMoveSteps);

                for(int mcnt=0; mcnt<SMMoveNum; mcnt++)
                {
                    short [] BatchedImgDat = GetBatchedImgDat(BatchedImgNum);
                    float [] LocResults = QC_STORM_Plug.lm_LocBatchedImg(BatchedImgDat, BatchedImgNum);
                    
                    // PSF width as a reference to adjust focus knob
    /*                
oInf[0] = Mean_PSFWidth;
oInf[1] = -Mean_SNR;

oInf[2] = Mean_PSFWidth - Mean_SNR;
oInf[3] = Mean_PSFWidth*(-Mean_SNR);
oInf[4] = FluoNum;*/
                    
                    PSFVaryInf[mcnt] = LocResults[ZCorrMode]; // mean PSF width
 

                    // move to a new depth
                    ZDriftSMMoveAndWait(SMMoveSteps);
                }
                
                int MinPos = FindArryMinPos(PSFVaryInf);
                
//                MMScript.message("move target:"+Float.toString(PSFVaryInf[MinPos])+" "+Integer.toString(MinPos));
                
                // move to the place with minimum PSF width
                ZDriftSMMoveAndWait(-SMMoveSteps*(SMMoveNum-MinPos));

                
                if(MinPos == SMMoveNum_half)
                {
                    // already in the center focus
                    break;
                }
            }
        }
        
        int FindArryMinPos(float iArry[])
        {
            int MinPos=0;
            for (int cnt=0; cnt<iArry.length;cnt++)
            {
                if(iArry[cnt]<iArry[MinPos])
                {
                    MinPos=cnt;
                }
            }
            return MinPos;
        }    
        
        void ZDriftSMMoveAndWait(int MoveSteps)
        {
            QC_STORM_Plug.lm_ZDepthSMMove(MoveSteps);

            try {
                sleep(12);
            } catch (InterruptedException ex) {

            }
        }
        
        public int GetBatchedImgNum()
        {
            int BatchedImgNum = 2048 * 2048 / ImageWidth / ImageHigh;
            
            BatchedImgNum = Math.max(BatchedImgNum, 4);
            BatchedImgNum = Math.min(BatchedImgNum, 16);
            
            return BatchedImgNum;
        }
        
        public short [] GetBatchedImgDat(int BatchedImgNum)
        {
            int ImgPixelNum=ImageWidth*ImageHigh;
            short [] BatchedImgDat = new short[BatchedImgNum*ImgPixelNum];
            short [] RawImgDatS;

            try {
                int CurFrame=0;

                // Start burst acq
                double exposureMs = MMCore.getExposure();

                MMCore.startSequenceAcquisition(BatchedImgNum, 0, true);

                while ((MMCore.getRemainingImageCount() > 0) || MMCore.isSequenceRunning(MMCore.getCameraDevice())) 
                {
                    if (MMCore.getRemainingImageCount() > 0) 
                    {
                        TaggedImage taggedImg = MMCore.popNextTaggedImage();

                        Image image1 = MyStudio.data().convertTaggedImage(taggedImg, MyStudio.data().getCoordsBuilder().time(CurFrame).build(), null);

                        RawImgDatS = (short [])image1.getRawPixels();

                        // debug
                        // simulation, read from hard disk
//                        RawImgDatS = QC_STORM_Plug.GetSimuImage512x512();

                        // group several images into one
                        System.arraycopy(RawImgDatS, 0 ,BatchedImgDat, ImgPixelNum*CurFrame, ImgPixelNum);

                        CurFrame++;
                    }
                    else 
                    {
                        // Wait for another image to arrive.
                        MMCore.sleep(Math.max(0.4 * exposureMs, 1));
                    }
                }       
                MMCore.stopSequenceAcquisition();
                
            } catch (Exception ex) {

            }
            return BatchedImgDat;
        }    
}

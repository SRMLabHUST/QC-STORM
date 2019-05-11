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
import mmcorej.TaggedImage;
import org.micromanager.ScriptController;
import org.micromanager.Studio;
import org.micromanager.data.Image;

/**
 *
 * @author luchang li
 */

// capture small batch of images for drift correction ...

public class QC_STORM_SmallBatchImageAcq {

    private Studio MyStudio;
    private CMMCore MMCore;
    private ScriptController gui;

    QC_STORM_Configurator myConfigurator;


    int ImageWidth;
    int ImageHigh;
        
    QC_STORM_SmallBatchImageAcq(Studio iStudio, QC_STORM_Configurator iConfigurator)           
    {
        MyStudio=iStudio;
        MMCore=MyStudio.getCMMCore();
        gui=MyStudio.scripter();


        myConfigurator = iConfigurator;

        ImageWidth = (int) MMCore.getImageWidth();
        ImageHigh = (int) MMCore.getImageHeight();
         
    }
    

    public int GetBatchedImgNum()
    {
        int BatchedImgNum = 2000 * 2048 / ImageWidth / ImageHigh;

        BatchedImgNum = Math.max(BatchedImgNum, 2);
        BatchedImgNum = Math.min(BatchedImgNum, 10);

        return BatchedImgNum;
    }

    public short [] GetBatchedImgDat(int BatchedImgNum)
    {
        int ImgPixelNum = ImageWidth*ImageHigh;
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
    
    public short [] GetBatchedImgDat_Default()
    {
        int BatchedImgNum = GetBatchedImgNum();
        
        return GetBatchedImgDat(BatchedImgNum);
    }
    
    public float [] GetLocResultsOfBatchedImageDat_Default(int ImageNumOffset)
    {
        if(ImageNumOffset<0)ImageNumOffset=0;
        
        int BatchedImgNum = GetBatchedImgNum() + ImageNumOffset;
        
        short [] BatchedImgDat = GetBatchedImgDat(BatchedImgNum);
        
        float [] LocResults = QC_STORM_Plug.lm_LocBatchedImg(BatchedImgDat, BatchedImgNum);
        
        return LocResults;
    }
    
}

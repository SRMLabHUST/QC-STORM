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


import org.micromanager.data.Image;
import org.micromanager.data.Processor;
import org.micromanager.data.ProcessorContext;
import org.micromanager.Studio;
import ij.ImagePlus;
import ij.gui.ImageWindow;
import ij.io.FileSaver;
import ij.io.Opener;
import ij.process.ColorProcessor;
import ij.process.FloatProcessor;
import java.awt.image.ColorModel;
import mmcorej.CMMCore;
import org.micromanager.ScriptController;



public class QC_STORM_Processor  extends Processor{
    
    private Studio studio_;
    CMMCore mmc;
    ScriptController gui;


    public QC_STORM_Configurator MyConfigurator; // reference to configurator

    public QC_STORM_Parameters.LocalizationPara CurLocPara;
                
    private int MyProcessorId = 0; // avoid system clean after a class is destroyed when users' processor is running
    
    // create time of this processor
    public String CreateTimeIdxStr;
    public String NamePostFix;
    
    public int ImageWidthI, ImageHighI;
    public int SRImageWidthI, SRImageHighI;
    
    public String ResultSavePathStr;
        


    // display of super resolution image
    // for 2d image, it's a float image, 3d is a color(RGB) image
    public ImagePlus CurSRImagePlus;
    public ImageWindow CurSRImageWindow;
    public FloatProcessor CurSRImageProcessor2D;
    public ColorProcessor CurSRImageProcessor3D;


    volatile boolean IsAcquisitionB = true;    
    volatile boolean IsSaved = false;
    
    boolean IsMultiFOVAcqB = false;
    

    QC_STORM_Processor(Studio studio, QC_STORM_Configurator iConfigurator, String iNamePostFix)
    {
        studio_ = studio;
        mmc = studio_.getCMMCore();
        gui = studio_.scripter();
        
        MyProcessorId = QC_STORM_Plug.SendProcessorID();
        
        MyConfigurator = iConfigurator;
        
        CurLocPara = MyConfigurator.GetLocalizationPara();
        
        NamePostFix = iNamePostFix;
        
        if(NamePostFix.length()>0)
        {
            IsMultiFOVAcqB = true;
        }else
        {
            IsMultiFOVAcqB = false;
        }
        
        IsAcquisitionB = true;
        IsSaved = false;
        
        // use time to mark each acquistion
        CreateTimeIdxStr = QC_STORM_Parameters.GetCreateTimeIdx() + NamePostFix;
        
        CurSRImagePlus = new ImagePlus();
        CurSRImagePlus.setTitle("rec SR image " + CreateTimeIdxStr);
        
        
        ImageWidthI = (int) mmc.getImageWidth();
        ImageHighI = (int) mmc.getImageHeight();
                   
  
        // set super resolution image size
        SRImageWidthI = (int) (ImageWidthI*CurLocPara.RenderingPixelZoom);
        SRImageHighI = (int) (ImageHighI*CurLocPara.RenderingPixelZoom);
        SRImageWidthI = (SRImageWidthI+3)/4*4;
        SRImageHighI = (SRImageHighI+3)/4*4;
        
     
        // start loc thread
        
        MyConfigurator.SendLocalizationPara();
        
        ResultSavePathStr = MyConfigurator.GetResultsSavePath();
        QC_STORM_Plug.lm_SetSavePath(ResultSavePathStr.toCharArray());
        QC_STORM_Plug.lm_SetAcquisitionPara(CreateTimeIdxStr.toCharArray());
        
        
        QC_STORM_Plug.lm_StartLocThread();

        
        gui.message("create new processor");
        
    }
    
    @Override
    public void cleanup(ProcessorContext pc)
    {
        gui.message("system clean processor");
        
        CleanProcessing();
    }
    public void CleanProcessing()
    {
        gui.message("manual clean processor");

        if(IsSaved==false)
        {
            IsSaved=true;
            
            // avoid system clean after a class is destroyed when users' processor is running
            if(MyProcessorId == QC_STORM_Plug.lm_GetProcessorID())
            {
                CloseAndSaveSRImage();

            }      
        }       
    }

    public void CloseAndSaveSRImage()
    {

        QC_STORM_Plug.CloseAndWaitCurrentLocThread();
        
        IsAcquisitionB = false;
        
        /*         
        // don't save image for multi fov acq for time savting
        if(!IsMultiFOVAcqB)
        {
        }
        */
        SaveAndDispSRImage();
        
        // close super-resolution image for save memory buffer of ImageJ platform
        if(IsMultiFOVAcqB)
        {
            CurSRImagePlus.close();
        }
    }
    
    public void SaveAndDispSRImage()
    {
        GetSRImage();
                
        CurSRImageWindow = new ImageWindow(CurSRImagePlus); 
        CurSRImagePlus.draw();     
        
        
        FileSaver ResultTifSaver = new FileSaver(CurSRImagePlus);
        
        String MultiFitStr;
		if(CurLocPara.MultiEmitterFitEn>0) MultiFitStr = "_M";
		else MultiFitStr = "_S";

		String ConsecFitStr;
		if(CurLocPara.ConsecutiveFitEn>0) ConsecFitStr = String.format("_Consec%.0fnm", CurLocPara.ConsecFilterRadius);
		else ConsecFitStr = "";
        
        String SaveImgName = String.format("loc_result%dD%d_%s%s%s_rend%.2fnm.tif", CurLocPara.LocType+2,CurLocPara.RegionSize, CreateTimeIdxStr,MultiFitStr,ConsecFitStr, CurLocPara.RenderingPixelSize);
        
        ResultTifSaver.saveAsTiff(ResultSavePathStr + SaveImgName);
            

    }
    
    public void GetSRImage()
    {
        if(CurLocPara.LocType==0)
        {
            CurSRImageProcessor2D = new FloatProcessor(SRImageWidthI, SRImageHighI);

            // get 2d display rended image
            float RecImgF[] = QC_STORM_Plug.lm_GetSMLMImage();  

            if(RecImgF.length>10)
            {
                ColorModel hotCM = QC_STORM_Parameters.GetHotColorModel();
                
                CurSRImageProcessor2D.setPixels(RecImgF);
                CurSRImageProcessor2D.setColorModel(hotCM);
                CurSRImagePlus.setProcessor(CurSRImageProcessor2D);

                CurSRImagePlus.setDisplayRange(0, QC_STORM_Plug.lm_GetMaxDispVal());
                CurSRImagePlus.updateImage();

            }
            else
            {
                CurSRImagePlus.setProcessor(CurSRImageProcessor2D);
                CurSRImagePlus.updateImage();
            }
        }
        else
        {
            CurSRImageProcessor3D = new ColorProcessor(SRImageWidthI, SRImageHighI);

            // get 3d display rended image
            int RecImgI[]=QC_STORM_Plug.lm_GetSMLMImage3D();  
            if(RecImgI.length>10)
            {
                CurSRImageProcessor3D.setPixels(RecImgI);

            }
            CurSRImagePlus.setProcessor(CurSRImageProcessor3D);
            CurSRImagePlus.updateImage();

        }
    }
    
    public void ProcessAImage(Image image)
    {
        
//        gui.message("processImage");
        short pImgData[] = (short[])image.getRawPixels();

        
        ///////////////////////////////////////////////////////////////////////////////
        // simulation, read from hard disk
//        pImgData = QC_STORM_Parameters.GetSimuImage512x512();
        
        // send image to c++ and processed by GPU
        QC_STORM_Plug.lm_FeedImageData(pImgData, 1);
   
    }
    
    @Override
    public void processImage(Image image, ProcessorContext pc) {

        ProcessAImage(image);
 
        pc.outputImage(image);      
    }

}

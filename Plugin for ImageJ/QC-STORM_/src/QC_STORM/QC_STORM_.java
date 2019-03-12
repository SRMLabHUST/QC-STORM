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


//package QC_STORM;


import ij.ImagePlus;
import ij.ImageStack;
import ij.gui.ImageWindow;
import ij.io.FileSaver;
import ij.plugin.filter.PlugInFilter;
import static ij.plugin.filter.PlugInFilter.DOES_16;
import static ij.plugin.filter.PlugInFilter.NO_IMAGE_REQUIRED;
import ij.process.*;

import java.awt.image.ColorModel;
import static java.lang.Thread.sleep;
import javax.swing.JOptionPane;


public class QC_STORM_ implements PlugInFilter{
    // parameters from imageJ 
    public ImagePlus CurImagePlus;
    public ImageProcessor CurImageProcessor;
    public ImageStack CurImageStack;
    
    // image info
    public int ImageWidthI, ImageHighI, FrameNum;
    public int SRImageWidthI,SRImageHighI;
    
    public String RawImgName;
    public String SaveFileName;
    
    
    // display of super resolution image
    // for 2d image, it's a float image, 3d is a color(RGB) image
    public ImagePlus CurSRImagePlus;
    public volatile ImageWindow CurSRImageWindow;
    public FloatProcessor CurSRImageProcessor2D;
    public ColorProcessor CurSRImageProcessor3D;
    public ColorModel hotCM;

    public ImageLocalizationThread LocThread=null;
    public  RerendThread CurRerendThread=null;

    
    // receive rendered image thread
    public String SRImgName;
    public RecSRImgThread recImgThread;
         
    public QC_STORM_Configurator MyConfigurator;

    
    short [] BatchedImgDat;
    int BatchedImgNum;
    int PixelNumI;
    
    
 
    public StatInfoDisplay CurStatInfoDisplay=null;
    public volatile int StatDispSelI;
    
    
    // wrapper of CPP & CUDA library
    public native static void lm_SetImagePara(int ImageWidth, int ImageHigh, int SRImageWidth, int SRImageHigh, int FrameNum, char ImageName[]);

    public native static void lm_SetLocPara(float KAdc, float Offset, float QE, int ROISize, int LocType, int MultiEmitterFitEn, int ConsecutiveFitEn, float ConsecFilterRadius, float RawPixelSize, float RenderPixelZoom, float SNR_th);
    public native static void lm_SetLocPara3D(float MinZDepth, float MaxZDepth, float ZDepthCorrFactor, float p4, float p3, float p2, float p1, float p0);
    
    public native static void lm_SetStatInfSelection(int DispSel, int OnTimeEn, int SpatialResolutionEn);

    public native static void lm_StartLocThread();
    public native static void lm_StopLocThread();
    public native static int lm_GetMaxDispVal();

    public native static void lm_FeedImageData(short ImgDataS[], int FrameNumI);
    

    public native static float [] lm_GetSMLMImage();
    public native static int [] lm_GetSMLMImage3D();
    
    public native static void lm_SetSpatialResolutionInf(int FramePerGroup, int IsHollowTube, float StructureSize, float RSCResolutionTh);
    
    public native static int [] lm_GetStatInfImageSize();
    public native static int [] lm_GetStatInfImage(int n);
    
    
    public native static int lm_IsLocFinish();

    // rerend image
    public native static void lm_SetRerendImagePara(char ImageName[], int IsDriftCorrectionI, int DriftCorrGroupFrameNum);
    public native static void lm_StartRerend();

    public native static int [] lm_GetRerendImageInf();

    public native static void lm_ReleaseRerendResource();


    
    
    @Override
    public int setup(String string, ImagePlus ip) {
        // initial
        System.loadLibrary ("QC-STORM_ CPPDLL");
        CurImagePlus=ip;

        
        if(CurImagePlus!=null)
        {
            CurImageStack = CurImagePlus.getImageStack();

            ImageWidthI = CurImagePlus.getWidth();// width of raw PALM image
            ImageHighI = CurImagePlus.getHeight();// height of raw PALM image
            FrameNum = CurImagePlus.getImageStackSize();// total frames number of raw PALM image
            RawImgName = CurImagePlus.getTitle();
            

            // batch processing
            BatchedImgNum=(2048*2048 / ImageWidthI / ImageHighI );
            BatchedImgNum = BatchedImgNum/2*2;
            
            if(BatchedImgNum<1)BatchedImgNum=1;
            if(BatchedImgNum>FrameNum)BatchedImgNum=FrameNum;

            //            BatchedImgNum=1;
            
            BatchedImgDat = new short[BatchedImgNum*ImageWidthI*ImageHighI];
            PixelNumI=ImageWidthI*ImageHighI;
            
        }

                
        CurSRImagePlus=new ImagePlus();
                           
        MyConfigurator=new QC_STORM_Configurator(this);
        MyConfigurator.setVisible(true);

		hotCM = QC_STORM_Parameters.GetHotColorModel();

        return NO_IMAGE_REQUIRED + DOES_16 + NO_CHANGES + SUPPORTS_MASKING; //open statement
    }

    @Override
    public void run(ImageProcessor ip) {
        CurImageProcessor=ip;
        
        while(MyConfigurator.PluginStatus!=2)
        {
            // wait the window closed
        }        
    }
    
    public void StartLocalization()
    {
        if((LocThread!=null)&&(LocThread.isAlive()))
        {
            return;
        }
        if(CurImagePlus==null)
        {
            return;
        }
        
        MyConfigurator.GetLocalizationPara();

        // set super resolution image size
        SRImageWidthI = (int) (ImageWidthI*MyConfigurator.LocPara.RenderingPixelZoom);
        SRImageHighI = (int) (ImageHighI*MyConfigurator.LocPara.RenderingPixelZoom);
        SRImageWidthI = (SRImageWidthI+3)/4*4;
        SRImageHighI = (SRImageHighI+3)/4*4;
        
        
        String MultiFitStr;
		if(MyConfigurator.LocPara.MultiEmitterFitEn>0) MultiFitStr = "_M";
		else MultiFitStr = "_S";

		String ConsecFitStr;
		if(MyConfigurator.LocPara.ConsecutiveFitEn>0) ConsecFitStr = String.format("_Consec%.0fnm", MyConfigurator.LocPara.ConsecFilterRadius);
		else ConsecFitStr = "";
        
        
        SRImgName = RawImgName;
        SRImgName = SRImgName.replace(".tif", "");

        SRImgName=String.format("%s_result%dD%d%s%s_rend%.2fnm.tif", SRImgName, MyConfigurator.LocPara.LocType+2, MyConfigurator.LocPara.RegionSize,MultiFitStr,ConsecFitStr, MyConfigurator.LocPara.RenderingPixelSize);


        CurSRImagePlus.setTitle(SRImgName);

        
        LocThread = new ImageLocalizationThread();
        LocThread.start();
        
       
    }
    private void InitSRDisplay()
    {
//            String str1=String.format("sr %d %d", SRImageWidthI, SRImageHighI);
//            JOptionPane.showMessageDialog(null, "InitSRDisplay "+Integer.toString(MyConfigurator.LocPara.LocType), "InitSRDisplay finish!", JOptionPane.PLAIN_MESSAGE);
           
        // for display rendered image, pre create them since create takes a lot of time
        if(MyConfigurator.LocPara.LocType == 0)
        {
            CurSRImageProcessor2D = new FloatProcessor(SRImageWidthI, SRImageHighI);
            CurSRImageProcessor3D = null;
        }
        else
        {
            CurSRImageProcessor2D = null;
            CurSRImageProcessor3D = new ColorProcessor(SRImageWidthI, SRImageHighI);
        }

    }
    public class ImageLocalizationThread extends Thread
    {
        @Override
        public void run()
        {
            
            SaveFileName=MyConfigurator.GetResultsSavePath() + RawImgName;


            // start to localize
            long startTime = System.currentTimeMillis(); //

            MyConfigurator.SendLocalizationPara();

            lm_SetImagePara(ImageWidthI, ImageHighI, SRImageWidthI, SRImageHighI, FrameNum, SaveFileName.toCharArray());

            lm_StartLocThread();    

            // stastic inf display
            StatDispSelI = MyConfigurator.LocPara.StatDispSel;

            InitSRDisplay();
            // receive rendered image        
            recImgThread = new RecSRImgThread();
            recImgThread.start();


            int curf=0;
            short RawImgDatS[];

            int ProcGroupNum = (FrameNum + BatchedImgNum-1)/BatchedImgNum;
            
            int gcntI;
            int CurFrameNum;
            int CurFramePos=1;
            int icntI;
            
            for(gcntI=0;gcntI<ProcGroupNum;gcntI++)
            {
                if(gcntI == ProcGroupNum-1)
                {
                    CurFrameNum = FrameNum - gcntI*BatchedImgNum;
                }else
                {
                    CurFrameNum = BatchedImgNum;
                }
                
                for(icntI=0;icntI<CurFrameNum;icntI++)
                {
                    CurFramePos = gcntI*BatchedImgNum + icntI + 1;
                    
                    RawImgDatS = (short [])CurImageStack.getPixels(CurFramePos);
                    
                    System.arraycopy(RawImgDatS, 0 ,BatchedImgDat, PixelNumI*icntI, PixelNumI);
                    
                }
                lm_FeedImageData(BatchedImgDat, CurFrameNum);
                
                CurImagePlus.setPosition(CurFramePos);
            }
            
            CurImagePlus.setPosition(FrameNum);


            while(lm_IsLocFinish()==0)
            {
                // wait GPU localization finsih
            }

            long endTime = System.currentTimeMillis();
            long runtimeL = endTime-startTime;
            

            try {
                recImgThread.join();
            } catch (InterruptedException ex) {
            }
            

            FileSaver ResultTifSaver = new FileSaver(CurSRImagePlus);
            String SaveImgName = MyConfigurator.GetResultsSavePath()+SRImgName;
            ResultTifSaver.saveAsTiff(SaveImgName);
            
            
            String str1=String.format("loc finish,runtime=%d ms", runtimeL);
            JOptionPane.showMessageDialog(null, str1, "loc finish!", JOptionPane.PLAIN_MESSAGE);
            
       
            lm_StopLocThread(); // finish localization and release resources
          
        }
    }
    public void StartRerend()
    {
        if((CurRerendThread!=null))
        {
            if(CurRerendThread.isAlive())
            {
                return;
            }
        }       
        MyConfigurator.GetLocalizationPara();
        MyConfigurator.SendLocalizationPara();        
        
        
        CurSRImagePlus.setTitle("Rerend image");

//        JOptionPane.showMessageDialog(null, "begin to rend", "loc finish!", JOptionPane.PLAIN_MESSAGE);
        
        CurRerendThread = new RerendThread();
        CurRerendThread.start();
           
    }

    public class RerendThread extends Thread
    {
        @Override
        public void run()
        {
            int [] ImgInf;

            lm_SetRerendImagePara(MyConfigurator.RerendLocDataPath.toCharArray(), MyConfigurator.DriftCorrEnableI, MyConfigurator.DriftCorrGroupFrameNum);  
            lm_StartRerend();
            

            while(true)
            {
                ImgInf = lm_GetRerendImageInf();
                
                if(ImgInf[0]>0)
                {
                    break;
                }  
            }

            // get image size from cpp codes
            ImageWidthI = ImgInf[1];
            ImageHighI = ImgInf[2];
            SRImageWidthI = ImgInf[3];
            SRImageHighI = ImgInf[4];


            InitSRDisplay();
            
             
            // receive rendered image        
            recImgThread = new RecSRImgThread();
            recImgThread.start();
            

            try {
                recImgThread.join();
            } catch (InterruptedException ex) {

            }
            

 //           while(recImgThread.isAlive());

            FileSaver ResultTifSaver=new FileSaver(CurSRImagePlus);
            String SaveImgName_PostFix=String.format("_rerend%.2fnm.tif", MyConfigurator.LocPara.RenderingPixelSize);

            String SaveImgName=MyConfigurator.RerendLocDataPath.substring(0, MyConfigurator.RerendLocDataPath.length()-4) + SaveImgName_PostFix;
            ResultTifSaver.saveAsTiff(SaveImgName);
                      
            
            lm_ReleaseRerendResource();
            
            JOptionPane.showMessageDialog(null, "rend finish", "loc finish!", JOptionPane.PLAIN_MESSAGE);
        }
    }

    public class RecSRImgThread extends Thread
    {
        @Override
        public void run()
        {
            boolean IsBreakB=false;
            boolean IsDispaly=false;

            while(true)
            {
                if(lm_IsLocFinish()==0)
                {
                    IsBreakB=false;
 //                   continue;// for speed testing
                }else
                {
                    IsBreakB=true;
                }              

                if(MyConfigurator.LocPara.LocType == 0)
                {
                    // get 2d display rended image
                    float RecImgF[] = lm_GetSMLMImage();  

                    if(RecImgF.length>10)
                    {
         
                        CurSRImageProcessor2D.setPixels(RecImgF);
                        CurSRImageProcessor2D.setColorModel(hotCM);
                        CurSRImagePlus.setProcessor(CurSRImageProcessor2D);
                        CurSRImagePlus.updateImage();
                        CurSRImagePlus.setDisplayRange(0, lm_GetMaxDispVal());


                        if(IsDispaly==false)
                        {
                            IsDispaly = true;
                            CurSRImageWindow = new ImageWindow(CurSRImagePlus); 

                        }
                        CurSRImagePlus.draw();   
                        CurSRImagePlus.updateImage();
                        CurSRImagePlus.setDisplayRange(0, lm_GetMaxDispVal());
                    }
                }
                else
                {
                    // get 3d display rended image
                    int RecImgI[] = lm_GetSMLMImage3D();  
                    if(RecImgI.length>10)
                    {
                        CurSRImageProcessor3D.setPixels(RecImgI);
                        CurSRImagePlus.setProcessor(CurSRImageProcessor3D);


                        if(IsDispaly==false)
                        {
                            IsDispaly=true;
                            CurSRImageWindow=new ImageWindow(CurSRImagePlus); 
                        }
                        CurSRImagePlus.draw();  
                        CurSRImagePlus.updateImage();
                    }
                }

                
                if(CurStatInfoDisplay==null)
                {
                    CurStatInfoDisplay = new StatInfoDisplay(StatDispSelI);
                }
                // get current para set
                MyConfigurator.LocPara = MyConfigurator.GetLocalizationPara();
                StatDispSelI = MyConfigurator.LocPara.StatDispSel;
            
                CurStatInfoDisplay.UpdateDisplay(StatDispSelI);
              
            
                if(IsBreakB)break;
                try {
                    sleep(2000);
                } catch (InterruptedException ex) {
                    
                }
           }
        }
    }
    
    public class StatInfoDisplay
    {
        // define in CPP codes
        int AxesImgWidth = 890;
        int AxesImgHigh = 560;

        public final int DispInfNum = 15;

        public ColorProcessor [] DispColorProcessor = new ColorProcessor[DispInfNum];
        public ImagePlus      [] DispImagePlus = new ImagePlus[DispInfNum];
		public ImageWindow    [] DispImageWindow = new ImageWindow[DispInfNum];

		String [] DispInfName={
                    "Total photon distribution",
                    "Localization precision distribution",
                    "SNR (peak photon to background) distribution",
                    "PSF width (Pixel) distribution",

                    "Total photon variation",
                    "Localization precision variation",
                    "Ontime variation", 
                    "SNR variation", 
                    "PSF width variation", 
                    "Localization density 2D variation",
                    "Background variation",
                    
                    "Spatial resolution variation (nm)",
                    "Nyquist resolution variation (nm)",
                    "Dimension FD variation",
                    "Localization density FD variation"
                    
		};
		
        public void UpdateDisplay(int DispSel)
        {
            int curselI=0x0001;
            int cntI;
              
            for(cntI=0;cntI<DispInfNum;cntI++)
            {
                if((DispSel&curselI)!=0)
                {
                    int RecImgI[] = lm_GetStatInfImage(cntI);
                    
       
                    DispColorProcessor[cntI].setPixels(RecImgI);
                    DispImagePlus[cntI].setProcessor(DispColorProcessor[cntI]);
                    DispImagePlus[cntI].updateImage();
                    DispImagePlus[cntI].draw();

                    if(DispImageWindow[cntI].isClosed())
                    {
                        DispImageWindow[cntI] = new ImageWindow(DispImagePlus[cntI] );
                    }
                    if(DispImageWindow[cntI].isVisible()==false)
                    {
                        DispImageWindow[cntI].setVisible(true);
                    }                        
                }
                else
                {
                    DispImageWindow[cntI].setVisible(false);
                }
                curselI<<=1;
            }
        }

        public StatInfoDisplay(int DispSel)
        {
            int ImageSize[] = lm_GetStatInfImageSize();
            
            AxesImgWidth = ImageSize[0];
            AxesImgHigh = ImageSize[1];

        
            int curselI=0x0001;
            int cntI;
            for(cntI=0;cntI<DispInfNum;cntI++)
            {
                DispColorProcessor[cntI] = new ColorProcessor(AxesImgWidth, AxesImgHigh);
                DispImagePlus[cntI] = new ImagePlus(DispInfName[cntI], DispColorProcessor[cntI]);
                DispImageWindow[cntI] = new ImageWindow(DispImagePlus[cntI] );

                if((DispSel&curselI)!=0)
                {
                    DispImageWindow[cntI].setVisible(true);
                }
                else
                {
	                DispImageWindow[cntI].setVisible(false);
                }
                curselI<<=1;
                
            }
        }
    }


}


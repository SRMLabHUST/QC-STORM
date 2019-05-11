/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package hust.whno.SMLM;

import mmcorej.CMMCore;
import org.micromanager.ScriptController;
import org.micromanager.Studio;
import org.micromanager.data.Datastore;
import org.micromanager.data.Image;
import org.micromanager.display.DisplayWindow;

/**
 *
 * @author lucha
 */
public class Calibration3DAcq {

    private Studio studio_;
    CMMCore mmc;
    ScriptController gui;
    
    Datastore MyDatastore;
    DisplayWindow Mydisp;
    
    String ImageSavePath;
    
    QC_STORM_Configurator MyConfigurator; // reference to configurator
    
    int PlaneNum_Half;
    int PlaneNum;
    
    
    Calibration3DAcq(Studio studio, QC_STORM_Configurator iConfigurator, int iPlaneNum_Half)
    {
        studio_ = studio;
        mmc = studio_.getCMMCore();
        gui = studio_.scripter();
        
        MyConfigurator = iConfigurator;
        
        PlaneNum_Half = iPlaneNum_Half;
        PlaneNum = 2*PlaneNum_Half+1;
        
        
        String ImageNamePreFix = String.format("CalibrationImage_Steps%d_", MyConfigurator.GetZMoveSteps());
        ImageSavePath = MyConfigurator.GetResultsSavePath() + ImageNamePreFix + QC_STORM_Parameters.GetCreateTimeIdx();
            
    }
    
    public void StartAcq()
    {
        MyConfigurator.SetCalibAcqEnable(false);
        
        Calib3DImagesAcq ImageAcq = new Calib3DImagesAcq();
        ImageAcq.start();
        
    }
    public class Calib3DImagesAcq extends Thread
    {
        @Override
        public void run()
        {
            try {
            
                MyDatastore = studio_.data().createMultipageTIFFDatastore(ImageSavePath, true, true);

                studio_.displays().manage(MyDatastore);

                Mydisp = studio_.displays().createDisplay(MyDatastore);

                
                // z steps move
                MyConfigurator.ZStepsMoveNPlane(-PlaneNum_Half);
                
                for(int cnt=0; cnt<PlaneNum;cnt++)
                {
                    java.util.List<Image> Images = studio_.live().snap(false);
                    
                    Image image1 = Images.get(0);
                    
                    Image image2 =  image1.copyAtCoords(studio_.data().getCoordsBuilder().time(cnt).build());
                    
                            
                    MyDatastore.putImage(image2);

                    MyConfigurator.ZStepsMoveNPlane(1);
                }
                
                MyDatastore.close();
                //MyDatastore.freeze();
                //MyDatastore.setSavePath(ImageSavePath);
                
            } catch (Exception ex) {

            }
            
            MyConfigurator.SetCalibAcqEnable(true);
        }
    }
}

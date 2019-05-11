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

import org.micromanager.ScriptController;
import org.micromanager.Studio;




/**
 *
 * z drift correction by move the focus knob to minimize the PSF width, optional the molecule number and SNR
 */

public class QC_STORM_ZDriftCorrection  extends Thread
{
    
    QC_STORM_SmallBatchImageAcq SmallBatchImageAcq;
    
    QC_STORM_Configurator myConfigurator;

    private ScriptController gui;

    int ZCorrMode;
    int StepZoom;

    // move step motor to correct z drift
    int SMMoveNum_half;
    int SMMoveNum;
    int IterationNum;
    
    float LocDensity_FocalPlane;
    
                
    QC_STORM_ZDriftCorrection(Studio iStudio, QC_STORM_Configurator iConfigurator, int iZCorrMode, int iStepZoom, int iSMMoveNum_half, int iIterationNum)
    {

        SmallBatchImageAcq = new QC_STORM_SmallBatchImageAcq(iStudio, iConfigurator);
        
        myConfigurator = iConfigurator;
        
        
        gui = iStudio.scripter();

        ZCorrMode = iZCorrMode;
        StepZoom = iStepZoom;
        //
        SMMoveNum_half = iSMMoveNum_half;
        SMMoveNum = SMMoveNum_half*2+1;
        IterationNum = iIterationNum;
        
        LocDensity_FocalPlane = 0;
    }

    @Override
    public void run()
    {
                
        // set localization parameters 
        myConfigurator.SendLocalizationPara();
        
        
        myConfigurator.SetFeedbackDevicePort();

        int SMMoveSteps = myConfigurator.GetZMoveSteps() ;

        SMMoveSteps = SMMoveSteps * StepZoom;


        float [] ZDriftCorrInfVary = new float[SMMoveNum];
        
        float [] LocDensityVary = new float[SMMoveNum];


        for(int itcnt=0; itcnt<IterationNum; itcnt++)
        {
            // pre move
            ZDriftSMMoveAndWait(-SMMoveNum_half*SMMoveSteps);

            for(int mcnt = 0; mcnt < SMMoveNum; mcnt++)
            {
                
                float [] LocResults = SmallBatchImageAcq.GetLocResultsOfBatchedImageDat_Default(0);

//                gui.message("QC_STORM_ZDriftCorrection:"+String.format("%f", LocResults[0]));

                // maximize
                ZDriftCorrInfVary[mcnt] = LocResults[ZCorrMode]; //  
                
                LocDensityVary[mcnt] = LocResults[QC_STORM_Plug.LocInfID_LocDensity]; //  
                
                
                if(mcnt<SMMoveNum-1)
                {
                    // move to a new depth
                    ZDriftSMMoveAndWait(SMMoveSteps);   
                }

                if(!myConfigurator.IsZDriftCorrActive())
                {
                    break;
                }
            }

            int MaxPos = FindArryMaxPos(ZDriftCorrInfVary);

//                MMScript.message("move target:"+Float.toString(PSFVaryInf[MinPos])+" "+Integer.toString(MinPos));

            // move to the place with max ZDriftCorrInfVary
            ZDriftSMMoveAndWait(-SMMoveSteps*(SMMoveNum - 1 - MaxPos));

            LocDensity_FocalPlane = LocDensityVary[MaxPos];
        }
    }
    
    int FindArryMaxPos(float iArry[])
    {
        int FindPos = iArry.length / 2;

        for (int cnt = 0; cnt < iArry.length; cnt++)
        {
            if(iArry[cnt] > iArry[FindPos])
            {
                FindPos = cnt;
            }
        }
        return FindPos;
    }    

    void ZDriftSMMoveAndWait(int MoveSteps)
    {
        QC_STORM_Plug.lm_ZDepthSMMove(MoveSteps);

    }
    
    float GetCurrentDensity()
    {
        return LocDensity_FocalPlane;
    }
}

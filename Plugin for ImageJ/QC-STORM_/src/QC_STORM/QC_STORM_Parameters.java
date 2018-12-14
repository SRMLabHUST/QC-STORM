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


// package QC_STORM;

import java.awt.image.ColorModel;
import java.awt.image.IndexColorModel;
import java.io.File;





public class QC_STORM_Parameters {
    
    // parameters for localization
    public static class LocalizationPara
    {   
		float KadcF    = 0.45f;
   		float OffsetF  = 100;
   		float QEF 	   = 0.72f;
   		int RegionSizeI= 7;

        float RawImgPixelSize = 100;
        float RenderingPixelSize =10;
        float RenderingPixelZoom = 10;
        
        
        int LocTypeI=0;
        int ConsecutiveFitEnI=0;
        float ConsecFilterRadius = 80;
        
    	float MinZDepthF = 0;
        float MaxZDepthF = 1000;

    	float p4F = 0;
    	float p3F = -4.217f;
    	float p2F = 13.08f;
    	float p1F = -278.79f;
    	float p0F = 0;
        
        // double helix 3d
        float MeanDistance=10.1f;
        float DistanceTh  =0.9f;
        int RotateType    =0;

        // rendering
        float SNR_thF = 5; // used in rendering
        
        // spatial resolution
        int ImagePerGroup = 50;
        float StructureSize2D = 40;
        
        
        // statastic information display enable
        int StatDispSelI=0;
        
        int OnTimeEn=1;
        int SpatialResolutionEn=1;


    }
        
    public static String SelectDisk()
    {
        final int DiskNum=8; 
        String [] DiskName={
        "D",
        "E",
        "F",
        "G",
        "H",
        "I",
        "J",
        "C"
        };
        
        int i;
        for(i=0;i<DiskNum;i++)
        {
            String CurDiskName=DiskName[i]+":\\";
            File f = new File(CurDiskName);
            if(f.exists()){
                return CurDiskName;
            }  
        }
        return null;
    }    
    public static ColorModel GetHotColorModel()
    {
                // construct hot colormap
		byte[] r = new byte[256], g = new byte[256], b = new byte[256];
		for (int q = 0; q < 85; q++) {
			r[q] = (byte) (3 * q);
			g[q] = 0;
			b[q] = 0;
		}
		for (int q = 85; q < 170; q++) {
			r[q] = (byte) (255);
			g[q] = (byte) (3 * (q - 85));
			b[q] = 0;
		}
		for (int q = 170; q < 256; q++) {
			r[q] = (byte) (255);
			g[q] = (byte) 255;
			b[q] = (byte) (3 * (q - 170));
		}
		 return new IndexColorModel(3, 256, r, g, b);
    }
        
}

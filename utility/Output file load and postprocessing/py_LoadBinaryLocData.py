import numpy as np

FilePath = 'D:\MMStack_Pos0_7.ome_result2D7_M.txt'

# parameters are 
# peak intensity (photon), x (pixel), y (pixel), z (nm), PSFSigmaX (pixel), PSFSigmaY (pixel), Total intensity (photon), background (photon), SNR (peak to background e-), CRLBx (nm), CRLBy (nm), frame

LocArray = np.fromfile(FilePath, dtype = np.float32).reshape(-1, 12)
        

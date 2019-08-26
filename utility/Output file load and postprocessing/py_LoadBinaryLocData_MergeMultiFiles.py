import numpy as np
import os

FolderPath = r'G:\test images\test'

# these files may not all exist, the codes auto terminates
FilePath = [
    'MMStack_Pos-1.ome_result2D7_M.txt',
    'MMStack_Pos-1_1.ome_result2D7_M.txt',
    'MMStack_Pos-1_2.ome_result2D7_M.txt',
    'MMStack_Pos-1_3.ome_result2D7_M.txt',
    'MMStack_Pos-1_4.ome_result2D7_M.txt',
    'MMStack_Pos-1_5.ome_result2D7_M.txt',
    'MMStack_Pos-1_6.ome_result2D7_M.txt',
    'MMStack_Pos-1_7.ome_result2D7_M.txt',
    'MMStack_Pos-1_8.ome_result2D7_M.txt',
    'MMStack_Pos-1_9.ome_result2D7_M.txt',
    'MMStack_Pos-1_10.ome_result2D7_M.txt',
    'MMStack_Pos-1_11.ome_result2D7_M.txt',
    'MMStack_Pos-1_12.ome_result2D7_M.txt',
]

FileNum = len(FilePath)


# parameters are 
# peak intensity (photon), x (pixel), y (pixel), z (nm), PSFSigmaX (pixel), PSFSigmaY (pixel), Total intensity (photon), background (photon), SNR (peak to background e-), CRLBx (nm), CRLBy (nm), frame


for i in range(FileNum):
    
    FileWholePath = FolderPath + '\\' + FilePath[i]
    
    if(i==0):
        LocArray = np.fromfile(FileWholePath, dtype = np.float32).reshape(-1, 12)
        EndFrame = LocArray[-1,-1]
    else:
        if(os.path.exists(FileWholePath)):
            LocArray_t = np.fromfile(FileWholePath, dtype = np.float32).reshape(-1, 12)
            LocArray_t[:,-1] += EndFrame
            
            LocArray= np.concatenate((LocArray, LocArray_t), axis = 0)
            EndFrame = LocArray[-1,-1]
        else:
            break

print('total frame:', LocArray[-1,-1])
      
import os
import glob
import numpy as np
import shutil
import tifffile as tf
''' using tifffile.py by Christoph (Version: 2014.02.05)
    (http://www.lfd.uci.edu/~gohlke/code/tifffile.py.html)
'''

basePath='E:\\2014 Data\\2014-12-02_timelapse_python'
os.chdir(basePath)

N=9
Nx_orig=3
Ny_orig=3

Nx=Ny_orig
Ny=Nx_orig

index=range(1, N+1)
locationMatrix=np.reshape(index, (Ny_orig,Nx_orig)).T
locationMatrixNew=np.reshape(index, (Nx_orig,Ny_orig))

print locationMatrix
print locationMatrixNew


for filename in glob.glob("alexa*.tif"):
    print filename

    fileNumber = int(filename[5:8])

    if fileNumber < N+1:
    
        fileNumStr=str(fileNumber).zfill(3)
      
        alexa=tf.imread('alexa'+fileNumStr+'.tif')
##        cy=tf.imread('cy'+fileNumStr+'.tif')
##        dapi=tf.imread('dapi'+fileNumStr+'.tif')
##        gfp=tf.imread('gfp'+fileNumStr+'.tif')
##        nir=tf.imread('nir'+fileNumStr+'.tif')
##        tmr=tf.imread('tmr'+fileNumStr+'.tif')
##        trans=tf.imread('trans'+fileNumStr+'.tif')

        all=np.zeros((7,1024,1024), dtype=np.uint16)
        all[0,:,:]=alexa
##        all[1,:,:]=cy
##        all[2,:,:]=dapi
##        all[3,:,:]=gfp
##        all[4,:,:]=nir
##        all[5,:,:]=tmr
##        all[6,:,:]=trans

        #find array position of image number, relies on unique numbers
        location=np.nonzero(locationMatrix == fileNumber)
        loc_x=location[0]
        loc_y=location[1]

        outFileNumber=locationMatrixNew[loc_x[0], loc_y[0]]
        
        outFileNumStr=str(outFileNumber).zfill(3)
        outFilename=outFileNumStr+'out.tif'
        tf.imsave(outFilename, all)

def copyFile(src, dest):
    try:
        shutil.copy(src, dest)
    # eg. src and dest are the same file
    except shutil.Error as e:
        print('Error: %s' % e)
    # eg. source or destination doesn't exist
    except IOError as e:
        print('Error: %s' % e.strerror)

for filename in glob.glob("*out.tif"):
    print filename

    fileNumber = int(filename[0:3])
    folderNumber=fileNumber


    folderPath=basePath + '\\' + str(folderNumber).zfill(3) + '\\'
      
    location=np.nonzero(locationMatrix == fileNumber)
    loc_x=location[0]
    loc_y=location[1]
    x=loc_x[0]
    y=loc_y[0]
    posX=[x-1, x, x+1, x-1, x+1, x-1, x, x+1]
    posY=[ y-1, y-1, y-1, y, y, y+1, y+1, y+1]
    
    # if surrounded on all sides        
    if min(posX)>-1 and min(posY)>-1 and max(posX)<Nx_orig and max(posY)<Ny_orig: #fileNumber==5:

        if not os.path.exists(folderPath):
            os.makedirs(folderPath)

            imgName=filename
            
            copyName=folderPath + filename
            copyFile(imgName, copyName)
            
        for q in range(0,8):
            posFileNumberStr=str(locationMatrixNew[posX[q], posY[q]]).zfill(3)
            posFilename=posFileNumberStr+'out.tif'
            copyFile(posFilename, folderPath + posFilename)
        
        
 
    
    

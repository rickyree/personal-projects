#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 20 22:56:57 2021

@author: hjlee
"""


import os 
import numpy as np
from natsort import natsorted
import nibabel as nib
import imageio
import shutil
import dicom2nifti as d2n
 


#pathways specific to my pc: 

#main path
os.chdir("/Users/hjlee/urop")
path_main = os.getcwd()

#pathway storing all 2d dicom files 
slices = path_main + "/slices"

#pathway to store 3d nifti files
images = path_main + "/images"


#####

#function to combine 2d dicom images into a 3d nifti image 

#inputs: 
#path_to_slices: path to the 2d dicom images
#path_to_images: path in which the 3d nifti image will be stored 

#in path_to_slices there should be dicom files, all of which are 
#concatenated into a single nifti file. 

def combine (path_to_slices, path_to_images):
    os.chdir(path_to_slices)
    

    #creating a temporary pathway to store only the dicom files: 
    if not os.path.exists('temporary'): 
        os.makedirs('temporary')
    path_to_temporary = path_to_slices + '/temporary'
    
    
    #copying only dicom files into temporary pathway
    for n in os.listdir(path_to_slices): 
        if (n.endswith('.dcm')): 
            
            shutil.copy(n, path_to_temporary)
            
    #stacking all the 2d dicom files in temporary into a 3d nifti file (to obtain .nii.gz, set compression as true)        
    d2n.convert_directory(path_to_temporary, path_to_images, compression = False)
    
    #remove temporary pathway
    shutil.rmtree(path_to_temporary)

#####    


combine(slices, images)





#End of function. 


################################################


#Below is the code that I have been practicing on. Please ignore: 
    
    
    
    
    


def combining (path_to_slices, path_to_images):
    os.chdir(path_to_slices)
    
    #A list to store the get_fdata arrays of the 2d images
    combined = []
    
    #A list of files to concatenate in path_to_slices
    files_to_analyze = []
    
    #obtaining files that are dicom 
    for count, n in enumerate(os.listdir(path_to_slices)): 
        if (n.endswith('.dcm')): 
            files_to_analyze.append(count)
    
    print(files_to_analyze)
    #for each dicom image: 
    for n in files_to_analyze: 
        
        image = nib.read(n)
        
        
        #obtain get_fdata arrays for each dicom image and store them in combined: 
        dimensions = image.get_fdata()
        combined.append(dimensions)
        print(dimensions.shape)
    
    #concatenate the dicom images into a 3d nifti image: 
    hi = nib.concat_images(combined, check_affines=False)
    
    #save file 
    nib.save(hi, os.path.join(path_to_images,'image' +'.nii'))
    
    
#perform function on 2 aforementioned pathways: 
#combining(slices, images)  






def slicing (path_to_images, path_to_slices_folder):
  os.chdir(path_to_images)
  
  
  
  i=0
  new_names=[]
  original_names=[]
  os.chdir(path_to_images)
  for filename in natsorted(os.listdir(path_to_images)):
    if filename.endswith('.nii'):
      original_names.append(filename)
      os.rename(filename, str(i+1)+'.nii')
      renamed='%d.nii'%i
      new_names.append(renamed)
      i = i + 1

 
  os.chdir(path_to_images)
  i=0
  for filename in natsorted(os.listdir(path_to_images)):
   img=nib.load(filename)
   print(filename)
   im1=img.get_fdata() #creating a numpy array
   print(im1.shape)
   print(im1.dtype)
   Slices=im1.shape[2]
   for n in range(0, Slices):
    image=np.rot90(im1[:, :,n]) #rotates image 90 degrees to the left, if this is not needed you can skip this step and run the code below, changing image to im[:,:,n]
    imageio.imwrite(os.path.join(path_to_slices_folder,str(n+1)+'.nii'), image)
   i=i+1

#slicing(images, slices)


      
            













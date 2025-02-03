#!/usr/bin/env python

# -*- coding: utf-8 -*-
"""
Created on Mon Feb  6 17:45:16 2012

@author: Sokratis Makrogiannis
"""

import itk,sys


def add_fixed_value(input_image):

    adder  = itk.AddImageFilter.IF3IF3IF3.New()
    adder.SetInput1( input_image )
    adder.SetConstant2( 1.0e-10 )
    adder.Update()
    return adder.GetOutput()


def fat_ratio_function(fat_dixon_image, water_dixon_image, mask_image, image_type, label_value=[]):

    # Compute F/W ratio image    .
    # Add 1.0e-10 to avoid NaN issues.
    fat_dixon_image = add_fixed_value(fat_dixon_image)
    water_dixon_image = add_fixed_value(water_dixon_image)
    # divider_by_constant = itk.DivideImageFilter.IF3IF3IF3.New()
    # divider_by_constant.SetInput1(water_dixon_image)
    # divider_by_constant.SetConstant(1.5992)
    # divider_by_constant.Update()
    
    # Median filtering to remove outliers.
    radius = [1, 1, 1]
    median_filter = itk.MedianImageFilter[image_type, image_type].New()
    median_filter.SetRadius( radius )
    median_filter.SetInput( water_dixon_image ) # divider_by_constant.GetOutput()
    median_filter.Update()

    # Compute F/W image.
    divider = itk.DivideImageFilter.IF3IF3IF3.New(fat_dixon_image, median_filter.GetOutput())
    divider.Update()
    
    # Write output to file.    
    itk.write(divider.GetOutput(), "FatOverWater.nii.gz")

    # Threshold for removing outliers.
    #thresholdFilter = itk.ThresholdImageFilter.IF3.New()
    #thresholdFilter.SetInput( divider.GetOutput() )
    #thresholdFilter.ThresholdAbove( 1.0e06 );
    #thresholdFilter.SetOutsideValue(0);
    #thresholdFilter.Update()
    
    # Compute F+W image.
    adder = itk.AddImageFilter.IF3IF3IF3.New(fat_dixon_image, water_dixon_image); # divider_by_constant.GetOutput()
    adder.Update()
    # Write output to file.
    itk.write(adder.GetOutput(), "FatPlusWater.nii.gz")
    
    # Compute F/F+W image.
    divider2 = itk.DivideImageFilter.IF3IF3IF3.New(fat_dixon_image, adder.GetOutput())
    divider2.Update()
    # Write output to file.
    itk.write(divider2.GetOutput(), "FatOverFatPlusWater.nii.gz")
        
    # Convert mask pixel values to [0,1].
    thresholdFilter2 = itk.BinaryThresholdImageFilter.IF3IF3.New()
    if label_value == []:
        thresholdFilter2.SetUpperThreshold( float( 255 ) )
        thresholdFilter2.SetLowerThreshold( float( 1 ) )
    else:
        thresholdFilter2.SetUpperThreshold( float( label_value ) )
        thresholdFilter2.SetLowerThreshold( float( label_value ) )
        
    thresholdFilter2.SetInsideValue(1);
    thresholdFilter2.SetOutsideValue(0);
    thresholdFilter2.SetInput( mask_image )
    thresholdFilter2.Update()
    
    # Multiply mask with F/W image.
    multiplier = itk.MultiplyImageFilter.IF3IF3IF3.New(divider.GetOutput(), thresholdFilter2.GetOutput())
    multiplier.Update()
    
    # Multiply mask with F/F+W image.
    multiplier2 = itk.MultiplyImageFilter.IF3IF3IF3.New(divider2.GetOutput(), thresholdFilter2.GetOutput())
    multiplier2.Update()
    
    # Write output to file.
    itk.write( multiplier.GetOutput(), "FatOverWaterMasked.nii.gz")
    itk.write( multiplier2.GetOutput(), "FatOverFatPlusWaterMasked.nii.gz")
    
    #
    #  def readimage( imagetype, filename ):
    #    rd1 = itk.ImageFileReader[imagetype].New()
    #    rd1.SetFileName(filename)
    #    rd1.Update()
    #    image = rd1.GetOutput()
    #    return image
    
    return 0


if __name__ == "__main__":

    if len(sys.argv) < 4:
        print('Usage:\n' + str(sys.argv[0]) + ' <Dixon fat image> <Dixon water image> <Mask image> (<Label value>)')
        sys.exit(0)

    if len(sys.argv) == 5:
        label_value = sys.argv[4]

    # Read fat, water and mask images.
    image_type = itk.Image[itk.F, 3]
    reader = itk.ImageFileReader[image_type].New()
    reader.SetFileName( sys.argv[1] )
    reader.Update()
    fat_dixon_image = reader.GetOutput()
    reader = 0

    reader = itk.ImageFileReader[image_type].New()
    reader.SetFileName( sys.argv[2] )
    reader.Update()
    water_dixon_image = reader.GetOutput()
    reader = 0

    reader = itk.ImageFileReader[image_type].New()
    reader.SetFileName( sys.argv[3] )
    reader.Update()
    mask_image = reader.GetOutput()
    reader = 0

    if len(sys.argv) == 4:
        fat_ratio_function(fat_dixon_image, water_dixon_image, mask_image, image_type)
    elif len(sys.argv) == 5:
        fat_ratio_function(fat_dixon_image, water_dixon_image, mask_image, image_type, label_value)
    else:
        print('Usage:\n' + str(sys.argv[0]) + ' <Dixon fat image> <Dixon water image> <Mask image> (<Label value>)')
        sys.exit(0)



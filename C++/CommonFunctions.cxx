/*===========================================================================

Program:   Segmentation/Classification of thymus tissues from MRI.
Module:    $RCSfile: CommonFunctions.cxx,v $
Language:  C++
Date:      $Date: 2010/04/15 10:42:32 $
Version:   $Revision: 0.1 $
Author:    S. K. Makrogiannis
3T MRI Facility National Institute on Aging/National Institutes of Health.

=============================================================================*/

#include "CommonFunctions.h"
#include <vcl_algorithm.h>

#include <itkImageRegionIterator.h>
#include <itkImageRegionIteratorWithIndex.h>
#include <itkCurvatureAnisotropicDiffusionImageFilter.h>
#include <itkMedianImageFilter.h>
#include <itkRecursiveGaussianImageFilter.h>
#include <itkAddImageFilter.h>


// Produce an ROI around the seed point that will be used for labeling and segmentation.
FloatImageType::RegionType GenerateROI ( FloatImageType::Pointer inputVolume,
					 FloatImageType::PointType seedPoint,
					 std::vector<float> halfROISize )
{

  FloatImageType::RegionType seedbasedROI;
  FloatImageType::IndexType seedIndex;
  FloatImageType::SizeType roiSize;
  InputImageType::IndexType m1, m2;

  inputVolume->TransformPhysicalPointToIndex(seedPoint, seedIndex);
  

  for (int i = 0; i < FloatImageType::ImageDimension; ++i)
    {
      m1[i] = seedIndex[i] - 
	int(floor( ( halfROISize[i] / inputVolume->GetSpacing()[i]) + 0.5) );
      m2[i] = seedIndex[i] + 
	int(floor( ( halfROISize[i] / inputVolume->GetSpacing()[i]) + 0.5) );
      roiSize[i] = m2[i] - m1[i] + 1;
    }

  seedbasedROI.SetIndex( m1 );
  seedbasedROI.SetSize( roiSize );

  return seedbasedROI;

}


// Determine Median ROI Label.
int DetermineMedianROILabel( MaskImageType::Pointer InputVolume,
			     MaskImageType::RegionType seedbasedROI)
{
  
  typedef itk::ImageRegionIterator<MaskImageType> MaskImageIteratorType;
  MaskImageIteratorType itMask(InputVolume, seedbasedROI);
  float labelvectorLength = ( seedbasedROI.GetSize()[0] *
			      seedbasedROI.GetSize()[1] *
			      seedbasedROI.GetSize()[2] );
  vnl_vector<int> labelVector((int)labelvectorLength, 0);

  // Determine the vector's median value and and assign the label.
  int count = 0;
  for(itMask.GoToBegin();!itMask.IsAtEnd();++itMask) {
    labelVector[count] = itMask.Get();
    count++;
  }
  //std::cout << labelVector << std::endl;
  vnl_vector<int>::iterator medianIter = labelVector.begin() + labelVector.size()/2;
  vcl_nth_element(labelVector.begin(), medianIter, labelVector.end());
  int medianROILabel = *medianIter;

  return medianROILabel;

}


// Smooth out the input image.
FloatImageType::Pointer SmoothInputVolume( FloatImageType::Pointer inputVolume,
					   int denoisingMethod )
{
  FloatImageType::Pointer outputVolume = NULL;

  switch( denoisingMethod ) {
  case GAUSSIAN:
    {
      // Linear gaussian blurring.

      // Type definition of recursive gaussian filtering.
      typedef itk::RecursiveGaussianImageFilter<FloatImageType,
	FloatImageType> 
	GaussianFilterType;

      // Instantiate filter.
      GaussianFilterType::Pointer smoothFilter = GaussianFilterType::New();

      // Set input and parameters.
      smoothFilter->SetInput( inputVolume );
      smoothFilter->SetSigma( (double)SMOOTHINGSIGMA *
			      (double)SMOOTHINGITERATIONS ); // previously: SMOOTHINGSIGMA
      // smoothFilter->SetNormalizeAcrossScale( true );
      smoothFilter->Update();

      // Display execution info.
      std::cout << "Recursive Gaussian filtering completed." << std::endl;
      outputVolume = smoothFilter->GetOutput();
    }
    break;
  case DIFFUSION: 
    {
    // Non-linear diffusion-based smoothing.
    int nIterations = SMOOTHINGITERATIONS;
    double timestep = SMOOTHINGSIGMA;
    double conductanceParameter = 2.0;

    typedef itk::CurvatureAnisotropicDiffusionImageFilter<FloatImageType, 
	FloatImageType> 
      CurvatureAnisotropicDiffusionImageFilterType;
    CurvatureAnisotropicDiffusionImageFilterType::Pointer smoothFilter = 
      CurvatureAnisotropicDiffusionImageFilterType::New();
    smoothFilter->SetInput( inputVolume );
    smoothFilter->SetNumberOfIterations( nIterations );
    smoothFilter->SetTimeStep( timestep );
    smoothFilter->SetConductanceParameter( conductanceParameter );
    smoothFilter->Update();
    std::cout << "Non-linear diffusion filtering completed." << std::endl;
    outputVolume = smoothFilter->GetOutput();
    }
    break;
  case MEDIAN:
    {
      // Statistical median filtering.
      typedef itk::MedianImageFilter<FloatImageType, 
	FloatImageType>
	MedianImageFilterType;
      MedianImageFilterType::Pointer smoothFilter = 
	MedianImageFilterType::New();
      smoothFilter->SetInput( inputVolume );
      FloatImageType::SizeType indexRadius;
      indexRadius[0] = (unsigned int)MEDIANFILTERRADIUS;
      indexRadius[1] = (unsigned int)MEDIANFILTERRADIUS;
      indexRadius[2] = (unsigned int)MEDIANFILTERRADIUS;
      smoothFilter->SetRadius( indexRadius );
      smoothFilter->Update();
      std::cout << "Median filtering completed." << std::endl;
      outputVolume = smoothFilter->GetOutput();
    }
  }
  
  // Return smoothed image.
  return outputVolume;
}



// Compute non-suppressed image from Dixon images by
// adding the fat and water images.
FloatImageType::Pointer
AddDixonImages( FloatImageType::Pointer watersuppressedVolume,
		FloatImageType::Pointer fatsuppressedVolume ) 
{
  // Simply add the voxel values from the fat and water images.
  FloatImageType::Pointer outputVolume = NULL;
  typedef itk::AddImageFilter<FloatImageType, FloatImageType>
    AddImageFilterType;
  AddImageFilterType::Pointer sumFilter = 
    AddImageFilterType::New();
  sumFilter->SetInput1( watersuppressedVolume );
  sumFilter->SetInput2( fatsuppressedVolume );  
  sumFilter->Update();
  outputVolume = sumFilter->GetOutput();
  return outputVolume;
}



// String manipulations to separate the path from filename.

std::string ExtractDirectory( const std::string& path )
  {
    return path.substr( 0, path.find_last_of( PathSeparator ) +1 );
  }

std::string ExtractFilename( const std::string& path )
  {
  return path.substr( path.find_last_of( PathSeparator ) +1 );
  }

std::string ChangeExtension( const std::string& path, const std::string& ext )
  {
  std::string filename = ExtractFilename( path );
  return ExtractDirectory( path ) +filename.substr( 0, filename.find_last_of( '.' ) ) +ext;
  }


// ToDo
// Convert itk image to itk oriented image type.

// template<class T>
// ConvertITKImagetoITKOrientedImage(const T inputVolume)
// {
// 	// Have to convert itkimage to itkorientedimage type.
// 	this->m_outputMaskVolume = MaskImageType::New();
// 	this->m_outputMaskVolume->CopyInformation( tempVolume );
// 	this->m_outputMaskVolume->SetRegions( tempVolume->GetBufferedRegion() );
// 	this->m_outputMaskVolume->Allocate();
// 	typedef itk::ImageRegionConstIterator< ScalarImageKmeansImageFilterType::OutputImageType >
// 	  ImageConstIteratorType;
// 	typedef itk::ImageRegionIterator< MaskImageType >
// 	  ImageIteratorType;
// 	ImageConstIteratorType itSRC( tempVolume, 
// 				      tempVolume->GetBufferedRegion() );
// 	ImageIteratorType itDST( this->m_outputMaskVolume, 
// 				 this->m_outputMaskVolume->GetBufferedRegion() );
// 	while( !itSRC.IsAtEnd() )
// 	  { itDST.Set( itSRC.Get() );
// 	    ++itSRC;
// 	    ++itDST; }
// 	tempVolume = 0;

// }


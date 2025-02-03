/*===========================================================================

Program:   Segmentation/Classification of thymus tissues from MRI.
Module:    $RCSfile: ThymusSegmentationAlgorithms.cxx,v $
Language:  C++
Date:      $Date: 2011/04/11 10:42:32 $
Version:   $Revision: 1.0 $
Author:    S. K. Makrogiannis
3T MRI Facility National Institute on Aging/National Institutes of Health.

=============================================================================*/

#include "itkObject.h"
#include "itkImage.h"
#include "itkOrientedImage.h"
#include "itkVectorImage.h"
#include "itkScalarToArrayCastImageFilter.h"
#include "itkImageRegionIterator.h"
#include "itkImageRegionIteratorWithIndex.h"

#include "itkImageFileReader.h" 
#include "itkImageFileWriter.h" 

// Preprocessing.
//#include "itkSubtractImageFilter.h"
#include "itkCurvatureAnisotropicDiffusionImageFilter.h"
//#include "itkAdaptiveHistogramEqualizationImageFilter.h"
// #include <itkLaplacianSharpeningImageFilter.h>
#include "itkMinimumMaximumImageFilter.h"

// itk Clustering or classification.
#include "itkScalarImageKmeansImageFilter.h"
// #include "itkImageKmeansImageFilter.h"
#include "itkImageToVectorImageFilter.h"
#include "itkRescaleIntensityImageFilter.h"
#include "itkMRFImageFilter.h"
#include "itkDistanceToCentroidMembershipFunction.h"
#include "itkMinimumDecisionRule.h"
#include "itkImageClassifierBase.h"
#include "itkOtsuMultipleThresholdsImageFilter.h"

// itk Segmentation.
#include "itkVectorConfidenceConnectedImageFilter.h"
#include "itkSignedMaurerDistanceMapImageFilter.h" 
#include "itkGradientMagnitudeRecursiveGaussianImageFilter.h"
#include "itkVectorGradientMagnitudeImageFilter.h"
#include "itkFastMarchingImageFilter.h"
#include "itkSigmoidImageFilter.h"
#include "itkShapeDetectionLevelSetImageFilter.h"
#include "itkBinaryThresholdImageFilter.h"
#include "itkCastImageFilter.h"
#include "itkGeodesicActiveContourLevelSetImageFilter.h"

// itk Morphological operations.
// #include "itkBinaryErodeImageFilter.h"
// #include "itkBinaryDilateImageFilter.h"
// #include "itkBinaryBallStructuringElement.h"

// itk Connected component processing.
#include "itkConnectedComponentImageFilter.h"
#include "itkRelabelComponentImageFilter.h"
#include "itkBinaryThresholdImageFilter.h"
// Miscellaneous.
#include "CommonFunctions.h"
#include "ParzenEdgeEstimation.h"
#include "VectorParzenEdgeEstimation.h"
#include "ThymusSegmentationAlgorithms.h"

// Resampling
#include "itkResampleImageFilter.h"
#include "itkIdentityTransform.h"
#include "itkLinearInterpolateImageFunction.h"
#include "itkRecursiveGaussianImageFilter.h"

// ThymusSegmentation class members.

void 
ThymusSegmentationObject::
SetInputs( FloatImageType::Pointer nonSuppressedVolume,
	   FloatImageType::Pointer waterSuppressedVolume,
	   FloatImageType::Pointer fatSuppressedVolume )
{

  this->m_nonsuppressedVolume = nonSuppressedVolume;
  this->m_fatsuppressedVolume = fatSuppressedVolume;
  this->m_watersuppressedVolume = waterSuppressedVolume;
}



//! Smooth out NS, WS and NS images and form a vector image.

void
ThymusSegmentationObject::
FormNSWSFSVectorVolume()
{
  // Smooth out input volumes.
  FloatImageType::Pointer smoothedNonSuppressedVolume = 
    SmoothInputVolume( this->m_nonsuppressedVolume, MEDIANFILTER );
  FloatImageType::Pointer smoothedWaterSuppressedVolume = 
    SmoothInputVolume( this->m_watersuppressedVolume, MEDIANFILTER );
  FloatImageType::Pointer smoothedFatSuppressedVolume = 
    SmoothInputVolume( this->m_fatsuppressedVolume, MEDIANFILTER );

  // Generate vector image from NS, WS and FS volumes.
  // Instantiate data types for analysis.
  typedef itk::ScalarToArrayCastImageFilter< FloatImageType, VectorFloatImageType > VectorImageGeneratorType;
  VectorImageGeneratorType::Pointer inputVectorImageGenerator =  
    VectorImageGeneratorType::New() ;
  inputVectorImageGenerator->SetInput(0, smoothedNonSuppressedVolume) ;
  inputVectorImageGenerator->SetInput(1, smoothedWaterSuppressedVolume) ;
  inputVectorImageGenerator->SetInput(2, smoothedFatSuppressedVolume) ;
  inputVectorImageGenerator->Update() ;
  this->m_VectorVolume = inputVectorImageGenerator->GetOutput();
}


// //! Resample input images in a multi-scale level-set segmentation loop.
// void
// ThymusSegmentationObject::
// ResampleVolume()
// {
//   // Set the sampling factors.
//   const double factorX = 2;
//   const double factorY = 2;
//   const double factorZ = 2;

//   // Apply Gaussian filter before resampling
//   typedef itk::RecursiveGaussianImageFilter< 
//                                   InternalImageType,
//                                   InternalImageType > GaussianFilterType;

//   GaussianFilterType::Pointer smootherX = GaussianFilterType::New();
//   GaussianFilterType::Pointer smootherY = GaussianFilterType::New();
//   GaussianFilterType::Pointer smootherZ = GaussianFilterType::New();
//   smootherX->SetInput( caster->GetOutput() );
//   smootherY->SetInput( smootherX->GetOutput() );
//   smootherZ->SetInput( smootherY->GetOutput() );
//   const InputImageType::SpacingType& inputSpacing = inputImage->GetSpacing();

//   const double sigmaX = inputSpacing[0] * factorX;
//   const double sigmaY = inputSpacing[1] * factorY;
//   const double sigmaZ = inputSpacing[2] * factorZ;

//   smootherX->SetSigma( sigmaX );
//   smootherY->SetSigma( sigmaY );
//   smootherZ->SetSigma( sigmaZ );

//   smootherX->SetDirection( 0 );
//   smootherY->SetDirection( 1 );
//   smootherZ->SetDirection( 2 );

//   smootherX->SetNormalizeAcrossScale( false );
//   smootherY->SetNormalizeAcrossScale( false );
//   smootherZ->SetNormalizeAcrossScale( false );

//   // Set the transformation to identity.
//   typedef itk::IdentityTransform< double, Dimension >  TransformType;
//   TransformType::Pointer transform = TransformType::New();
//   transform->SetIdentity();
//   resampler->SetTransform( transform );

//   // Instantiate the interpolator.
//   typedef itk::LinearInterpolateImageFunction< 
//                                    InternalImageType, double >  InterpolatorType;
  
//   InterpolatorType::Pointer interpolator = InterpolatorType::New();
//   resampler->SetInterpolator( interpolator );
//   resampler->SetDefaultPixelValue( 0 ); // value for regions without source

//   // Set the spacing.
//     OutputImageType::SpacingType spacing;

//   spacing[0] = inputSpacing[0] * factorX;
//   spacing[1] = inputSpacing[1] * factorY;
//   spacing[2] = inputSpacing[2] * factorZ;

//   resampler->SetOutputSpacing( spacing );

//   // Preserve the origin and direction while undersampling.
//   resampler->SetOutputOrigin( inputImage->GetOrigin() );
//   resampler->SetOutputDirection( inputImage->GetDirection() );

//   // Set the output size.
//   InputImageType::SizeType   inputSize = 
//               inputImage->GetLargestPossibleRegion().GetSize();

//   typedef InputImageType::SizeType::SizeValueType SizeValueType;

//   InputImageType::SizeType   size;

//   size[0] = static_cast< SizeValueType >( inputSize[0] / factorX );
//   size[1] = static_cast< SizeValueType >( inputSize[1] / factorY );
//   size[2] = static_cast< SizeValueType >( inputSize[2] / factorZ );

//   resampler->SetSize( size );

//   // Set the resampler input and trigger execution.
//   resampler->SetInput( smootherZ->GetOutput() );
//   resample->Update();

// }



//! Separate the air (smallest intensity) from other tissue types
//! using the NS image.
void
ThymusSegmentationObject::
SeparateAirFromOtherTissues() {
  const int nClusters = 4;
  
  // Apply multiple Otsu thresholding to the NS or F+W image
  // with N=4.
  typedef itk::OtsuMultipleThresholdsImageFilter<FloatImageType, MaskImageType>
    TissueClustererType;
  TissueClustererType::Pointer TissueClusterer = TissueClustererType::New();
  TissueClusterer->SetNumberOfHistogramBins( 128 );
  TissueClusterer->SetNumberOfThresholds( nClusters  );
  TissueClusterer->SetInput(this->m_nonsuppressedVolume);
  TissueClusterer->Update();

  // Compute average intensities of each label.
  MaskImageType::Pointer clusterLabelImage = 
    TissueClusterer->GetOutput();

  typedef itk::ImageRegionIteratorWithIndex<MaskImageType> MaskImageIteratorWithIndexType;
  MaskImageIteratorWithIndexType itMask(clusterLabelImage, 
					clusterLabelImage->GetBufferedRegion());
  vnl_vector<float> averagesVector(nClusters, 0.0), 
    cardinalitiesVector(nClusters, 0.0);

  // Label the group of smallest intensity as air and generate a
  // foreground/background mask.
  for(itMask.GoToBegin();!itMask.IsAtEnd();++itMask) {
    averagesVector[itMask.Get()] += 
      this->m_nonsuppressedVolume->GetPixel(itMask.GetIndex());
    cardinalitiesVector[itMask.Get()]++;
  }

  // Display cluster averages and cardinalities.
  std::cout << "Otsu multithreshold clustering." << std::endl;
  for(int i=0;i<nClusters;i++) {
    std::cout << "Cluster # " << i
	      << ", average intensity: "
	      << averagesVector[i] 
	      << ", voxel count "
	      << cardinalitiesVector[i]
	      << std::endl;
  }

  
  // Find minimum.
  MaskImageType::PixelType airLabel = 
    averagesVector.min_value();

  // Apply thresholding to generate the air mask.
  typedef itk::BinaryThresholdImageFilter<MaskImageType, 
    MaskImageType> 
    SelectAirClusterFilterType;
  SelectAirClusterFilterType::Pointer airRegionThresholdFilter = 
    SelectAirClusterFilterType::New();
  airRegionThresholdFilter->SetInput( clusterLabelImage );
  airRegionThresholdFilter->SetInsideValue( FOREGROUND );
  airRegionThresholdFilter->SetOutsideValue( BACKGROUND );
  airRegionThresholdFilter->SetLowerThreshold( airLabel );
  airRegionThresholdFilter->SetUpperThreshold( airLabel );  // number of regions we want detected.
  airRegionThresholdFilter->Update();
  std::cout << "Thresholding, done." 
	    << std::endl;

  // Copy pointer to the air mask image.
  this->m_AirMaskImage = airRegionThresholdFilter->GetOutput();

}



//! Generate a probabilistic edge map using Parzen kernels.
//  Scalar input version.
FloatImageType::Pointer
ThymusSegmentationObject::
ParzenEdgeMap(FloatImageType::Pointer inputImage)
{
  int halfwindowSize = (int) this->m_ParzenRadius;  // 2;
  double sigma = (double) this->m_ParzenBandwidth;  //30;

  int nSamples = (2*halfwindowSize+1)*
    (2*halfwindowSize+1)*
    (2*halfwindowSize+1);

  std::cout << "Parzen kernel based edge detection." 
	    << std::endl
	    << "Parzen radius = "
	    << this->m_ParzenRadius
	    << ", Parzen badwidth = "
	    << this->m_ParzenBandwidth
	    << std::endl;

  typedef itk::ConstNeighborhoodIterator< FloatImageType > NeighborhoodIteratorType;
  typedef itk::ImageRegionIterator< FloatImageType>        IteratorType;

  // Allocate memory for output image.
  FloatImageType::Pointer edgeMap = FloatImageType::New();
  edgeMap->CopyInformation( inputImage );
  edgeMap->SetRegions( inputImage->GetLargestPossibleRegion() );
  edgeMap->Allocate();
  edgeMap->FillBuffer( 0.0 );

  // Neighborhood iterator.
  NeighborhoodIteratorType::RadiusType radius;
  radius.Fill(halfwindowSize);
  NeighborhoodIteratorType it( radius, inputImage,
                               inputImage->GetRequestedRegion() );

  // Region iterator.
  IteratorType out(edgeMap, inputImage->GetRequestedRegion());
  typedef std::vector<double> DoubleVectorType;
  DoubleVectorType localVector(nSamples, 0.0);
  double totalMax = 0.0;
  for (it.GoToBegin(), out.GoToBegin(); !it.IsAtEnd(); ++it, ++out)
    {
      // For each voxel of the input volume,
      // generate a vector with the values within this 3D window.
	    
      for(int k=0;k<nSamples;k++) localVector[k] = it.GetPixel(k);
      // Estimate the probabilistic edge and
      // assign it to the output volume.
      double pmv = ParzenEdgeEstimation<DoubleVectorType>(localVector, 
							  (2*halfwindowSize+1), 
							  sigma);
      out.Set(pmv);
      if (totalMax<pmv) totalMax = pmv;
    }

  std::cout << "Overall pmv maximum = " 
	    << totalMax 
	    << std::endl;

  // Post-process values (normalize and complement to 0 if needed).
  for (out.GoToBegin(); !out.IsAtEnd(); ++out)
    {
      double edgeProbability = 255 * (1 - ( out.Get() / totalMax ) );  // for edge map.
      // double edgeProbability = 255 * (out.Get() / totalMax);  // for speed image.
      out.Set(edgeProbability);
    }

  return edgeMap;
}



//! Generate a probabilistic edge map using Parzen kernels.
// Works for vector images.
FloatImageType::Pointer
ThymusSegmentationObject::
VectorParzenEdgeMap(VectorFloatImageType::Pointer inputImage)
{
  int halfwindowSize = (int) this->m_ParzenRadius;  // 2;
  double sigma = (double) this->m_ParzenBandwidth;  //30;

  int nSamples = (2*halfwindowSize+1)*
    (2*halfwindowSize+1)*
    (2*halfwindowSize+1);

  std::cout << "Parzen kernel based edge detection." 
	    << std::endl
	    << "Parzen radius = "
	    << this->m_ParzenRadius
	    << ", Parzen badwidth = "
	    << this->m_ParzenBandwidth
	    << std::endl;

  typedef itk::ConstNeighborhoodIterator< VectorFloatImageType > NeighborhoodIteratorType;
  typedef itk::ImageRegionIterator< FloatImageType>        IteratorType;

  // Allocate memory for output image.

  FloatImageType::Pointer edgeMap = FloatImageType::New();
  edgeMap->CopyInformation( inputImage );
  edgeMap->SetRegions( inputImage->GetLargestPossibleRegion() );
  edgeMap->Allocate();
  edgeMap->FillBuffer( 0.0 );

  // Neighborhood iterator.
  NeighborhoodIteratorType::RadiusType radius;
  radius.Fill(halfwindowSize);
  NeighborhoodIteratorType it( radius, inputImage,
                               inputImage->GetRequestedRegion() );

  // Region iterator.
  IteratorType out(edgeMap, inputImage->GetRequestedRegion());
  // typedef std::vector<double> DoubleVectorType;
  // DoubleVectorType localVector(nSamples, 0.0);
  typedef vnl_matrix<double> DoubleArrayType;
  DoubleArrayType localSamples(nSamples, VoxelComponents, 0.0);

  double totalMax = 0.0;

  for (it.GoToBegin(), out.GoToBegin(); !it.IsAtEnd(); ++it, ++out)
    {
      // For each voxel of the input volume,
      // generate a vector with the values within this 3D window.
      for(int k=0;k<nSamples;k++) 
	for(int m=0;m<(int)VoxelComponents;m++)
	  localSamples[k][m] = it.GetPixel(k)[m];

//       std::cout << "kernel =" << std::endl
// 		<< localSamples << std::endl;

      // Estimate the probabilistic edge and
      // assign it to the output volume.
      double pmv = VectorParzenEdgeEstimation<DoubleArrayType>(localSamples, 
							       (2*halfwindowSize+1), 
							       sigma);
      out.Set(pmv);
//       std::cout << "pmv = " 
// 		<< pmv 
// 		<< std::endl;

      if (totalMax<pmv) totalMax = pmv;
    }

  std::cout << "Overall pmv maximum = " 
	    << totalMax 
	    << std::endl;

  // Post-process values (normalize and complement to 0 if needed).
  for (out.GoToBegin(); !out.IsAtEnd(); ++out)
    {
      double edgeProbability = 255 * (1 - ( out.Get() / totalMax ) );  // for edge map.
      // double edgeProbability = 255 * ( out.Get() / totalMax );  // for speed image.
      out.Set(edgeProbability);
    }

  return edgeMap;
}


// Use fast marching to initialize the level-set segmentation process.
MaskImageType::Pointer
ThymusSegmentationObject::
InitializeROIbyFastMarching(FloatImageType::Pointer speedImage) 
{
  // Use fast marching to initialize the segmentation process.

  // Declare fast marching filter type.
  typedef itk::FastMarchingImageFilter< FloatImageType,
    FloatImageType >    
    FastMarchingFilterType;
  FastMarchingFilterType::Pointer fastMarching = 
    FastMarchingFilterType::New();

  // Declare the node containers for initialization.
  typedef FastMarchingFilterType::NodeContainer        NodeContainer;
  typedef FastMarchingFilterType::NodeType             NodeType;
  NodeContainer::Pointer seeds = NodeContainer::New();
  NodeType node;
  // Set the seed value and index.
//   const double seedValue = (double) (-1)* this->m_seedbasedROI.GetSize()[0];
  const double seedValue = 0.0;
  node.SetValue( seedValue );
  FloatImageType::IndexType seedIndex;
  this->m_watersuppressedVolume->TransformPhysicalPointToIndex(this->m_seedPoint, 
							       seedIndex);
  node.SetIndex( seedIndex );

  // Add index to container.
  seeds->Initialize();
  seeds->InsertElement( 0, node );

  // Set the fast marching filter attributes.
  std::cout << "FM stop time: " 
	    << this->m_fastmarchingStoppingTime
	    << std::endl;
  fastMarching->SetTrialPoints( seeds );
  fastMarching->SetInput( speedImage );
  // Set stopping time for fast marching.
  fastMarching->SetStoppingValue( this->m_fastmarchingStoppingTime );
  fastMarching->SetOutputSize( this->m_watersuppressedVolume->GetBufferedRegion().GetSize() );
  fastMarching->Update();

  // Apply thresholding.
 typedef itk::BinaryThresholdImageFilter< FloatImageType,
    MaskImageType >    ThresholdingFilterType;
  ThresholdingFilterType::Pointer thresholder = ThresholdingFilterType::New();
  thresholder->SetLowerThreshold( 0.0 );
  thresholder->SetUpperThreshold( this->m_fastmarchingStoppingTime - 5.0 );
  thresholder->SetOutsideValue( BACKGROUND );
  thresholder->SetInsideValue( FOREGROUND );
  thresholder->SetInput( fastMarching->GetOutput() );
  thresholder->Update();

  std::string fmImageFilename = this->m_outputPath + 
    "/" + this->m_subjectID + "_" + SegmentationMethodString[this->m_segmentationMethod] + "_" + 
    "FM_Mask.nii";
  itk::ImageFileWriter<MaskImageType>::Pointer maskWriter = 
    itk::ImageFileWriter<MaskImageType>::New();
  maskWriter->SetInput( thresholder->GetOutput() ); 
  maskWriter->SetFileName( fmImageFilename );
  maskWriter->Update();
  maskWriter = 0;

  return thresholder->GetOutput();
}

// Apply k-means to vector image and rescale output according to the data type.

MaskImageType::Pointer 
ThymusSegmentationObject::
KMeansonVectorVolumes()
{
  // // Define and istantiate the K-Means filter.
  MaskImageType::Pointer outputMaskVolume = 0;
//   int useNonContiguousLabels = 0;
//   typedef itk::ImageKmeansImageFilter<VectorFloatImageType> KMeansFilterType;
//   KMeansFilterType::Pointer kMeansFilter = KMeansFilterType::New();
//   kMeansFilter->SetInput( this->m_VectorVolume );
//   kMeansFilter->SetUseNonContiguousLabels( useNonContiguousLabels );
//   KMeansFilterType::RealPixelType initialMean;
//   initialMean.Fill( 0.0 );
//   for( int cc = 0; cc < this->m_nClusters; ++cc )
//       kMeansFilter->AddClassWithInitialMean( initialMean );

//   kMeansFilter->Update();

//   // Store the centroids in our vector.
//   int nElements = VoxelComponents * this->m_nClusters;
//   this->m_kMeansCentroids.assign(nElements, 0.0);
//   for( int cc = 0; cc < nElements; ++cc ) {
//     this->m_kMeansCentroids[cc] =  kMeansFilter->GetFinalMeans()[cc];
//     //std::cout << this->m_kMeansCentroids[cc] << std::endl;
//   }

//   std::cout << "K-means clustering on NS, WS and FS volumes." 
//   	    << std::endl
//   	    << "# clusters = "
//   	    << this->m_nClusters
//   	    << std::endl;
//   std::cout << ", Final means are: " 
//   	    << std::endl;
//   std::cout << kMeansFilter->GetFinalMeans() 
//   	    << std::endl;

//   typedef itk::RescaleIntensityImageFilter<KMeansFilterType::OutputImageType, MaskImageType > RescalerType;
//   RescalerType::Pointer rescaler = RescalerType::New();
//   rescaler->SetInput( kMeansFilter->GetOutput() );
//   rescaler->SetOutputMinimum( BACKGROUND );
//   rescaler->SetOutputMaximum( FOREGROUND );
//   rescaler->Update();
//   outputMaskVolume = rescaler->GetOutput();

  return outputMaskVolume;
}



// Region growing on the vector image.

MaskImageType::Pointer 
ThymusSegmentationObject::
RegionGrowingOnVectorImage()
{
  // Set up the segmentation process.
  MaskImageType::Pointer outputMaskVolume;
  typedef  itk::VectorConfidenceConnectedImageFilter< VectorFloatImageType, 
    MaskImageType > ConnectedVectorFilterType;
  ConnectedVectorFilterType::Pointer vectorconfidenceConnected = 
    ConnectedVectorFilterType::New();
  
  std::cout << "Vector confidence connected region growing."
	    << std::endl;
  std::cout << "Region growing: multiplier= "
	    << this->m_rgMultiplier
	    << ", #Iterations= "
	    << this->m_rgIterations
	    << std::endl;

  vectorconfidenceConnected->SetInput( this->m_VectorVolume );
  vectorconfidenceConnected->SetMultiplier( this->m_rgMultiplier );
  vectorconfidenceConnected->SetNumberOfIterations ( this->m_rgIterations );
  vectorconfidenceConnected->SetReplaceValue( FOREGROUND );

  // Add seeds from ROI.
  typedef itk::ImageRegionIteratorWithIndex<VectorFloatImageType> VectorImageIteratorType;
  VectorImageIteratorType itVectorFloat(this->m_VectorVolume, this->m_seedbasedROI);
  VectorFloatImageType::IndexType idx;
  
  for(itVectorFloat.GoToBegin();!itVectorFloat.IsAtEnd();++itVectorFloat) 
    {  
      idx = itVectorFloat.GetIndex();
      vectorconfidenceConnected->AddSeed(idx);
    }

  // Fire up the segmentor.
  vectorconfidenceConnected->Update();
  outputMaskVolume = vectorconfidenceConnected->GetOutput();
  std::cout << "Region Growing on NS, WS and FS volumes, done." 
	    << std::endl;

  return outputMaskVolume;
}



// Level-set segmentation.

void
ThymusSegmentationObject::
LevelSetBasedSegmentation()
{
  // Set up the segmentation process.

  // Compute gradient magnitude.
  typedef itk::VectorGradientMagnitudeImageFilter<VectorFloatImageType> 
    GradientMagnitudeFilterType;
  typedef   itk::GradientMagnitudeRecursiveGaussianImageFilter<FloatImageType, 
    FloatImageType >  GradientFilterType;

  FloatImageType::Pointer gradientImage; 
  typedef itk::ImageRegionConstIterator< GradientMagnitudeFilterType::OutputImageType >
    VectorImageConstIteratorType;
  typedef itk::ImageRegionIterator< FloatImageType >
    FloatImageIteratorType;

  // Choose gradient method based on the image intensity dimensionality and
  // the desired algorithm (gradient magnitude or stochastic method).
  switch(this->m_segmentationMethod) {
  // Scalar pixel intesities.
  case WS_LEVEL_SETS: case WS_GEOD_ACT_CONT:
    {
      // Non-linear smoothing.
      FloatImageType::Pointer smoothedWaterSuppressedVolume = 
	SmoothInputVolume( this->m_watersuppressedVolume, MEDIANFILTER );

      switch(this->m_gradientMethod) {
      case CONVENTIONAL:
	{ 
	  // Compute the gradient magnitude from the WS image only.
	  GradientFilterType::Pointer  gradientMagnitude = GradientFilterType::New();
	  gradientMagnitude->SetSigma( this->m_gradientSigma ); // 0.0625
	  gradientMagnitude->SetInput( smoothedWaterSuppressedVolume );
	  gradientMagnitude->Update();
	  gradientImage = gradientMagnitude->GetOutput();
	}
	break;
      case PARZEN:
	{
	  // Edge map using Parzen estimators.
	  gradientImage = ParzenEdgeMap( smoothedWaterSuppressedVolume );
	}
	break;
      }
    }
    break;
    // Vector pixel intensities.
  case LEVEL_SETS_3D: case GEOD_ACT_CONT_3D:
    {
      switch(this->m_gradientMethod) {
      case CONVENTIONAL:
	{
	  GradientMagnitudeFilterType::Pointer 
	    vectorgradientMagnitude = 
	    GradientMagnitudeFilterType::New();
	  //     First off, calculate the gradient map on the previously smoothed volume.
	  vectorgradientMagnitude->SetInput( this->m_VectorVolume );
	  vectorgradientMagnitude->SetUseImageSpacing( true );
	  vectorgradientMagnitude->SetUsePrincipleComponents( true );
	  vectorgradientMagnitude->Update();
     
// 	  // Have to convert itkimage to itkorientedimage type.
// 	  gradientImage = FloatImageType::New();
// 	  gradientImage->CopyInformation( vectorgradientMagnitude->GetOutput() );
// 	  gradientImage->SetRegions( vectorgradientMagnitude->GetOutput()->GetBufferedRegion() );
// 	  gradientImage->Allocate();
// 	  VectorImageConstIteratorType itSRC( vectorgradientMagnitude->GetOutput(), 
// 					      vectorgradientMagnitude->GetOutput()->GetBufferedRegion() );
// 	  FloatImageIteratorType itDST( gradientImage, 
// 					gradientImage->GetBufferedRegion() );
// 	  while( !itSRC.IsAtEnd() )
// 	    { itDST.Set( itSRC.Get() );
// 	      ++itSRC;
// 	      ++itDST; }
	  gradientImage = vectorgradientMagnitude->GetOutput();
	}
	break;
      case PARZEN:
	{
	  // Edge map using Parzen estimators on all three images.
	  gradientImage = VectorParzenEdgeMap( this->m_VectorVolume );
	}
	break;
      }
    }

  }
  // Write gradient image to file.
  std::string gradientImageFilename = this->m_outputPath + "/" +
    this->m_subjectID + "_" + GradientMethodString[this->m_gradientMethod] + "_" + 
    SegmentationMethodString[this->m_segmentationMethod] + "_" + 
    "gradientmag.nii";
  itk::ImageFileWriter<FloatImageType>::Pointer floatWriter = 
    itk::ImageFileWriter<FloatImageType>::New();
  floatWriter->SetInput( gradientImage ); 
  floatWriter->SetFileName( gradientImageFilename );
  floatWriter->Update();
  floatWriter = 0;

  // Then apply sigmoid filter to produce edge potential image.
  std::cout << "Sigmoid Filter-- Alpha = "
	    << this->m_sigmoidAlpha
	    << ", Beta = "
	    << this->m_sigmoidBeta
	    << std::endl;

  typedef   itk::SigmoidImageFilter<FloatImageType, FloatImageType >  SigmoidFilterType;
  SigmoidFilterType::Pointer sigmoid = SigmoidFilterType::New();
  sigmoid->SetAlpha( this->m_sigmoidAlpha );   
  sigmoid->SetBeta( this->m_sigmoidBeta  );  
  sigmoid->SetOutputMinimum(  0.0  );
  sigmoid->SetOutputMaximum(  1.0  );
  sigmoid->SetInput( gradientImage ); 
  sigmoid->Update();

  // Initialize level set function.
  MaskImageType::Pointer roiVolume; 
  switch(this->m_levelsetInitializationMethod) {
  case SEED_BASED_ROI:
    {
      // If user asked for multi-scale analysis then use previous segmentation
      // to initialize current, else start from scratch.
      if( (this->m_nScales > 1) && (this->m_currentScale > 0) )
	{
	  roiVolume = this->m_outputMaskVolume;
	  std::cout << "Using previous output for initialization." << std::endl;
	}
      else
	{
	  // Use ROI to compute the distance map.

	  roiVolume = MaskImageType::New();
	  roiVolume->CopyInformation( this->m_nonsuppressedVolume );
	  roiVolume->SetRegions( this->m_nonsuppressedVolume->GetLargestPossibleRegion() );
	  roiVolume->Allocate();
	  roiVolume->FillBuffer( BACKGROUND );

	  itk::ImageRegionIteratorWithIndex<MaskImageType> 
	    itROI( roiVolume,   
		   this->m_seedbasedROI );
	  
	  // Generate an initialization of the level set from 
	  // the distance map values obtained from lumen mask.
	  for (itROI.GoToBegin(); 
	       !itROI.IsAtEnd(); 
	       ++itROI) 
	    itROI.Set( FOREGROUND );
	  std::cout << "Seed-based ROI generation." << std::endl;
	}
    }
      break;
  case FAST_MARCHING:
    {
      
      // If user asked for multi-scale analysis then use previous segmentation
      // to initialize current, else start from scratch.
      if( (this->m_nScales > 1) && (this->m_currentScale > 0) )
	{
	  roiVolume = this->m_outputMaskVolume;
	  std::cout << "Using previous output for initialization." << std::endl;
	}
      else
	{
	  // Use fast marching to initialize the segmentation process.
	  roiVolume = 
	    InitializeROIbyFastMarching( sigmoid->GetOutput() );
	  std::cout << "FM-based ROI generation." << std::endl;
	}
    }
    break;
  }


  // Calculate distance map from initial ROI.
  // Distance map filter type definition.
  FloatImageType::Pointer ROIDistanceOutput;
  typedef itk::SignedMaurerDistanceMapImageFilter<MaskImageType, FloatImageType> LabelDistanceTransformType;
  LabelDistanceTransformType::Pointer ROIDistFilter = LabelDistanceTransformType::New();
  ROIDistFilter->SetInput( roiVolume );
  ROIDistFilter->SetUseImageSpacing( true );
  // ROIDistFilter->SetInsideIsPositive( true );
  ROIDistFilter->SetSquaredDistance( true );  // false
  ROIDistFilter->Update();
  ROIDistanceOutput = ROIDistFilter->GetOutput();
  // ROIDistFilter = 0;

  // Instantiate caster because level set filters don't run on oriented images.
  typedef itk::CastImageFilter< itk::Image<float, 3>, FloatImageType > CastFilterType;
  CastFilterType::Pointer castImage = CastFilterType::New();
  FloatImageType::Pointer levelSetOutputImage = FloatImageType::New();

  // Run segmentation using shape detection or geodesic active contours.
  switch(this->m_segmentationMethod) {
  case WS_LEVEL_SETS: case LEVEL_SETS_3D:
  {
    typedef  itk::ShapeDetectionLevelSetImageFilter< FloatImageType, FloatImageType >    
      ShapeDetectionFilterType;

    // Run level set filter.
    std::cout << "Shape detection level sets." 
	      << std::endl
	      << "Propagation factor = "
	      << this->m_levelsetPropagationScalingFactor
	      << ", Curvature factor = "
	      << this->m_levelsetCurvatureScalingFactor
	      << std::endl;

    ShapeDetectionFilterType::Pointer shapeDetection = ShapeDetectionFilterType::New();                           
    shapeDetection->SetMaximumRMSError( this->m_levelsetMaximumRMSError );
    shapeDetection->SetNumberOfIterations( (int) this->m_levelsetMaximumIterations );
    shapeDetection->SetPropagationScaling(  this->m_levelsetPropagationScalingFactor );
    shapeDetection->SetCurvatureScaling( this->m_levelsetCurvatureScalingFactor ); 
    shapeDetection->SetIsoSurfaceValue( 0.0 );
    shapeDetection->SetInput( ROIDistanceOutput );
    shapeDetection->SetFeatureImage( sigmoid->GetOutput() );  // may use other speed volume modified by the distance.
    shapeDetection->Update();
    levelSetOutputImage = 
      shapeDetection->GetOutput();

    // Print out some useful information.
    std::cout << "Level set segmentation of thymus completed." << std::endl;
    std::cout << "Max. no. iterations: " << shapeDetection->GetNumberOfIterations() << "," <<
      "Max. RMS error: " << shapeDetection->GetMaximumRMSError() << std::endl;
    std::cout << "# of iterations: " << shapeDetection->GetElapsedIterations() << ","
              << "RMS change: " << shapeDetection->GetRMSChange() << std::endl;
  }
    break;
  case WS_GEOD_ACT_CONT: case GEOD_ACT_CONT_3D:
  {

    typedef itk::GeodesicActiveContourLevelSetImageFilter< FloatImageType, 
      FloatImageType > 
      GeodesicActiveContourFilterType;

    GeodesicActiveContourFilterType::Pointer geodesicActiveContours = 
      GeodesicActiveContourFilterType::New();

    geodesicActiveContours->SetMaximumRMSError( this->m_levelsetMaximumRMSError );
    geodesicActiveContours->SetNumberOfIterations( (int) this->m_levelsetMaximumIterations );
    geodesicActiveContours->SetPropagationScaling(  this->m_levelsetPropagationScalingFactor );
    geodesicActiveContours->SetCurvatureScaling( this->m_levelsetCurvatureScalingFactor ); 
    geodesicActiveContours->SetAdvectionScaling( this->m_levelsetAdvectionScalingFactor );
    geodesicActiveContours->SetIsoSurfaceValue( 0.0 );
    geodesicActiveContours->SetInput( ROIDistanceOutput );
    geodesicActiveContours->SetFeatureImage( sigmoid->GetOutput() );  // may use other speed volume modified by the distance.
    geodesicActiveContours->Update();

    std::cout << "Geodesic active contours." 
	      << std::endl
	      << "Propagation factor = "
	      << this->m_levelsetPropagationScalingFactor
	      << ", Curvature factor = "
	      << this->m_levelsetCurvatureScalingFactor
	      << ", Advection factor = "
	      << this->m_levelsetAdvectionScalingFactor
	      << std::endl;

    levelSetOutputImage = 
      geodesicActiveContours->GetOutput();

    // Print out some useful information.
    std::cout << "Level set segmentation of thymus completed." << std::endl;
    std::cout << "Max. no. iterations: " << geodesicActiveContours->GetNumberOfIterations() << "," <<
      "Max. RMS error: " << geodesicActiveContours->GetMaximumRMSError() << std::endl;
    std::cout << "# of iterations: " << geodesicActiveContours->GetElapsedIterations() << ","
              << "RMS change: " << geodesicActiveContours->GetRMSChange() << std::endl;
  }
  break;
  }

  // Type casting before thresholding.
  castImage->SetInput( levelSetOutputImage );
  castImage->Update();

  // Apply thresholding to the level sets
  // to get the segmentation result.
  typedef itk::BinaryThresholdImageFilter<FloatImageType, MaskImageType> ThresholdingFilterType;
  ThresholdingFilterType::Pointer thresholder = ThresholdingFilterType::New();
  thresholder->SetInput( castImage->GetOutput() );
  thresholder->SetOutsideValue(  BACKGROUND  );
  thresholder->SetInsideValue(  FOREGROUND );
  thresholder->SetLowerThreshold( -1000.0 );
  thresholder->SetUpperThreshold(     0.0 );
  thresholder->Update();
  // Free possibly already allocated volume.
  this->m_outputMaskVolume = 0;
  this->m_outputMaskVolume = thresholder->GetOutput();

  std::cout << "Level set-based segmentation, done." 
	    << std::endl;
}



// Run level sets in multiple edge map scales.
void 
ThymusSegmentationObject::
MultiscaleLevelSetBasedSegmentation() {
  //! Loop over the number of scales.
  std::cout << "Using multi-scale level set-based segmentation scheme." 
	    << std::endl;

  //! Record the original parameters.
  double originalParzenSigma = this->m_ParzenSigma;
  int originalLevelSetMaximumIterations = this->m_levelsetMaximumIterations;
  float originalGradientSigma = this->m_gradientSigma;
  const double parzensigmaBase = 2.0;
  const float gradientsigmaBase = 2.0;

  for(this->m_currentScale=0;
      this->m_currentScale<this->m_nScales;
      this->m_currentScale++) {
    
    //! Set the Parzen variables.
    this->m_ParzenSigma = 
      originalParzenSigma * 
      (double) std::pow(parzensigmaBase, (double)( (-1)*this->m_currentScale) );

    // //! Set the level-set iteration number.
    // this->m_levelsetMaximumIterations = 
    //   originalLevelSetMaximumIterations *
    //   std::pow(2, (this->m_currentScale - this->m_nScales + 1) );

    //! Set gradient magnitude sigma (works for scalar inputs only). 
    this->m_gradientSigma = 
      originalGradientSigma *
      (float)std::pow(gradientsigmaBase, (float)( (-1)*this->m_currentScale) );

    //! Run level sets.
    this->LevelSetBasedSegmentation();
  }
}



// Gaussian Mixture Modeling + Bayesian maximum ratio.

MaskImageType::Pointer 
ThymusSegmentationObject::
EM_GMM_BayesianClassification()
{
  // Set up the segmentation process.
  MaskImageType::Pointer outputMaskVolume;

  std::cout << "GMM using E-M optimization and Bayesian clustering on NS, WS and FS volumes, done." 
	    << std::endl;

  return outputMaskVolume;
}



// MRF Segmentation after clustering.

MaskImageType::Pointer 
ThymusSegmentationObject::
MRFSegmentationAfterClustering(MaskImageType::Pointer LabelImage)
{
  // Set up the segmentation process.
	MaskImageType::Pointer outputMaskVolume = 0;

//   // MRF filter.
//   // Set inputs and parameters.
//   typedef itk::MRFImageFilter< VectorFloatImageType, MaskImageType > MRFFilterType;
//   MRFFilterType::Pointer mrfFilter = MRFFilterType::New();
//   mrfFilter->SetInput( this->m_VectorVolume );
//   mrfFilter->SetNumberOfClasses( this->m_nClusters );
//   mrfFilter->SetMaximumNumberOfIterations( this->m_mrfIterations );
//   mrfFilter->SetErrorTolerance( 1e-7 );
//   mrfFilter->SetSmoothingFactor( this->m_mrfsmoothingFactor );

//   std::cout << "Markov Random Field filter." 
// 	    << std::endl
// 	    << "# clusters = "
// 	    << this->m_nClusters
// 	    << ", MRF smoothing factor = "
// 	    << this->m_mrfsmoothingFactor
// 	    << ", MRF iterations = "
// 	    << this->m_mrfIterations
// 	    << std::endl;

//   // Declare and instantiate classifier.
//   typedef itk::ImageClassifierBase<VectorFloatImageType, MaskImageType > 
//     SupervisedClassifierType;
//   SupervisedClassifierType::Pointer mrfClassifier = SupervisedClassifierType::New();

// //   // Check progress.
// //   typedef ShowProgressObject ProgressType;
// //   ProgressType progressWatch(classifierPointer);
// //   SimpleMemberCommand<ProgressType>::Pointer command;
// //   command = SimpleMemberCommand<ProgressType>::New();
// //   command->SetCallbackFunction(&progressWatch,
// //                                &ProgressType::ShowProgress);
// //   classifierPointer->AddObserver(itk::ProgressEvent(), command);

//   // Set decision rule.
//   typedef itk::MinimumDecisionRule DecisionRuleType;
//   DecisionRuleType::Pointer  classifierDecisionRule = DecisionRuleType::New();
//   mrfClassifier->SetDecisionRule( classifierDecisionRule.GetPointer() );
  
//   // Set membership functions.
//   typedef itk::Statistics::DistanceToCentroidMembershipFunction< VectorFloatPixelType > 
//     MembershipFunctionType;

//   typedef MembershipFunctionType::Pointer MembershipFunctionPointer;

//   //double meanDistance = 0.0;
//   vnl_vector<double> centroid( VoxelComponents, 0.0 ); 
//   int itCentroid;

//   for( int i=0; i < (int)this->m_nClusters; i++ )
//     {
//       MembershipFunctionPointer membershipFunction = 
// 	MembershipFunctionType::New();

//     // centroid[0] = atof( argv[i+numberOfArgumentsBeforeMeans] ); 
//       for(int j=0; j < (int)VoxelComponents; j++) {
// 	itCentroid = j + i * (int)VoxelComponents;
// 	centroid[j] = this->m_kMeansCentroids[ itCentroid ];
// 	//std::cout << itCentroid << ", " << this->m_kMeansCentroids[ itCentroid ] << std::endl;
//       }

//       //std::cout << centroid << std::endl;
//       membershipFunction->SetCentroid( centroid );

//       mrfClassifier->AddMembershipFunction( membershipFunction );
//       //meanDistance += static_cast< double > (centroid[0]);
//     }
//   //meanDistance /= this->m_nClusters;

//   mrfFilter->SetSmoothingFactor( this->m_mrfsmoothingFactor );

//   // Set the neighbors' number and weights.
//   mrfFilter->SetNeighborhoodRadius( this->m_mrfNeighborhoodRadius );
// //   std::vector< double > weights;
// //   weights.push_back(1.5);
// //   weights.push_back(2.0);
// //   weights.push_back(1.5);
// //   weights.push_back(2.0);
// //   weights.push_back(0.0); // This is the central pixel
// //   weights.push_back(2.0);
// //   weights.push_back(1.5);
// //   weights.push_back(2.0);
// //   weights.push_back(1.5);

// //   // Scale weights.
// //   double totalWeight = 0;
// //   for(std::vector< double >::const_iterator wcIt = weights.begin(); 
// //       wcIt != weights.end(); ++wcIt )
// //     {
// //     totalWeight += *wcIt;
// //     }
// //   for(std::vector< double >::iterator wIt = weights.begin(); 
// //       wIt != weights.end(); wIt++ )
// //     {
// //     *wIt = static_cast< double > ( (*wIt) * meanDistance / (2 * totalWeight));
// //     }

// //   mrfFilter->SetMRFNeighborhoodWeight( weights );

//   // Set the classifier and run.
//   mrfFilter->SetClassifier( mrfClassifier );
//   mrfFilter->Update();

// //   // Rescale MRF outputs according to data type.
// //   typedef MRFFilterType::OutputImageType  MRFOutputImageType;
// //   typedef itk::RescaleIntensityImageFilter<MRFOutputImageType, 
// //     MaskImageType >   RescalerType;

// //   RescalerType::Pointer intensityRescaler = RescalerType::New();
// //   intensityRescaler->SetOutputMinimum(   0 );
// //   intensityRescaler->SetOutputMaximum( 255 );
// //   intensityRescaler->SetInput( mrfFilter->GetOutput() );  
// //   intensityRescaler->Update();

// //   outputMaskVolume = intensityRescaler->GetOutput();

//   outputMaskVolume = mrfFilter->GetOutput();
//   outputMaskVolume->CopyInformation(this->m_nonsuppressedVolume);

//   std::cout << "Multi-parametric clustering and MRF on NS, WS and FS volumes, done." 
// 	    << std::endl;

   return outputMaskVolume;
}



// Connected components and/or moprhological operations
// after clustering to pick the thymus region.
void
ThymusSegmentationObject::
PostprocessClustering(
MaskImageType::Pointer outputMaskVolume
)
{

  // Morphological operations.

  // //   MaskImageType::Pointer morphologyClosedVolume =
  // //     BinaryMorphologicalClosing( outputMaskVolume,
  // //   				structureElementRadius );

  //   MaskImageType::Pointer morphologyOpenedVolume = 
  //     BinaryMorphologicalOpening( outputMaskVolume,
  //   				structureElementRadius);


  // Write clustering output to file.
  std::string clusterImageFilename = this->m_outputPath + 
    "/" + this->m_subjectID + "_" + SegmentationMethodString[this->m_segmentationMethod] + "_" + 
    "ClusteringOutput.nii";
  itk::ImageFileWriter<MaskImageType>::Pointer maskWriter = 
    itk::ImageFileWriter<MaskImageType>::New();
  maskWriter->SetInput( outputMaskVolume ); 
  maskWriter->SetFileName( clusterImageFilename );
  maskWriter->Update();
  maskWriter = 0;

  // Thymus region.
  LabelImageType::IndexType seedIdx; 
  outputMaskVolume->TransformPhysicalPointToIndex(this->m_seedPoint, 
						  seedIdx);
  std::cout << "Seed point: "
	    << this->m_seedPoint
	    << std::endl;
  std::cout << "Seed index: "
	    << seedIdx
	    << std::endl;


  // Determine Thymus label using ROI.
  int thymusRegionLabel = DetermineMedianROILabel( outputMaskVolume,
						   this->m_seedbasedROI);

  std::cout << "Thymus cluster label: "
	    << thymusRegionLabel
	    << std::endl;


  // Apply threshold to pick the thymus cluster.

  typedef itk::BinaryThresholdImageFilter<MaskImageType, 
    MaskImageType> 
    SelectThymusRegionFilterType;
  SelectThymusRegionFilterType::Pointer thymusRegionThresholdFilter = 
    SelectThymusRegionFilterType::New();
  thymusRegionThresholdFilter->SetInput( outputMaskVolume );
  thymusRegionThresholdFilter->SetInsideValue( FOREGROUND );
  thymusRegionThresholdFilter->SetOutsideValue( BACKGROUND );
  thymusRegionThresholdFilter->SetLowerThreshold( thymusRegionLabel );
  thymusRegionThresholdFilter->SetUpperThreshold( thymusRegionLabel );  // number of regions we want detected.
  thymusRegionThresholdFilter->Update();
  std::cout << "Thresholding, done." 
	    << std::endl;

  // Connected component labeling.

  typedef itk::ConnectedComponentImageFilter<MaskImageType, 
    MaskImageType, MaskImageType> ConnectedComponentLabelFilterType;
  ConnectedComponentLabelFilterType::Pointer labelMaskFilter = 
    ConnectedComponentLabelFilterType::New();
  labelMaskFilter->SetInput( thymusRegionThresholdFilter->GetOutput() );
  labelMaskFilter->SetMaskImage( thymusRegionThresholdFilter->GetOutput() );
  labelMaskFilter->SetFullyConnected( true );
  labelMaskFilter->Update();
  std::cout << "Connected component labeling, done." 
	    << std::endl;
  thymusRegionThresholdFilter = 0;

  thymusRegionLabel = DetermineMedianROILabel( labelMaskFilter->GetOutput(),
					       this->m_seedbasedROI);

  std::cout << "Thymus connected component label: "
	    << thymusRegionLabel
	    << std::endl;
  
  std::string connectedcomponentsFilename = this->m_outputPath + "/" + 
    this->m_subjectID + "_" + SegmentationMethodString[this->m_segmentationMethod] + "_" + 
    "ConnectedComponents.nii";
  itk::ImageFileWriter<MaskImageType>::Pointer labelWriter = 
    itk::ImageFileWriter<MaskImageType>::New();
  labelWriter = itk::ImageFileWriter<MaskImageType>::New();
  labelWriter->SetInput( labelMaskFilter->GetOutput() ); 
  labelWriter->SetFileName( connectedcomponentsFilename );
  labelWriter->Update();
  labelWriter = 0;

  // Second thresholding to pick the thymus component.

  thymusRegionThresholdFilter = SelectThymusRegionFilterType::New();
  thymusRegionThresholdFilter->SetInput( labelMaskFilter->GetOutput() );
  thymusRegionThresholdFilter->SetInsideValue( FOREGROUND );
  thymusRegionThresholdFilter->SetOutsideValue( BACKGROUND );
  thymusRegionThresholdFilter->SetLowerThreshold( thymusRegionLabel );
  thymusRegionThresholdFilter->SetUpperThreshold( thymusRegionLabel );  // number of regions we want detected.
  thymusRegionThresholdFilter->Update();

  this->m_outputMaskVolume = thymusRegionThresholdFilter->GetOutput();

  std::cout << "Second thresholding, done." 
	    << std::endl;

}


// This function calls all the segmentation algorithms.

void
ThymusSegmentationObject::
RunThymusSegmentationAlgorithms()
{

  //MaskImageType::Pointer outputMaskVolume;

  switch ( this->m_segmentationMethod )
    {
    case WS_KMEANS: // K-Means
      {
	// Non-linear smoothing.
	FloatImageType::Pointer SmoothedWaterSuppressedVolume = 
	  SmoothInputVolume( this->m_watersuppressedVolume, MEDIANFILTER );

	// Clustering (fcm or k-means).
	// Use 2 classes.
	typedef itk::ScalarImageKmeansImageFilter<FloatImageType> 
	  ScalarImageKmeansImageFilterType;

	ScalarImageKmeansImageFilterType::Pointer scalarImageKmeansImageFilter = 
	  ScalarImageKmeansImageFilterType::New();

	scalarImageKmeansImageFilter->SetInput( SmoothedWaterSuppressedVolume );
	//scalarImageKmeansImageFilter->SetInput( this->m_nonsuppressedVolume );
	scalarImageKmeansImageFilter->SetDebug( true );

	for (int i=1; i < this->m_nClusters; i++) {
	  scalarImageKmeansImageFilter->AddClassWithInitialMean( (float) 0.0 );
	}
	scalarImageKmeansImageFilter->Update();

	//this->m_outputMaskVolume = scalarImageKmeansImageFilter->GetOutput();
	ScalarImageKmeansImageFilterType::OutputImageType::Pointer tempVolume = 
	  (scalarImageKmeansImageFilter->GetOutput());

	// Have to convert itkimage to itkorientedimage type.
	MaskImageType::Pointer outputMaskVolume = MaskImageType::New();
	outputMaskVolume->CopyInformation( tempVolume );
	outputMaskVolume->SetRegions( tempVolume->GetBufferedRegion() );
	outputMaskVolume->Allocate();
	typedef itk::ImageRegionConstIterator< ScalarImageKmeansImageFilterType::OutputImageType >
	  ImageConstIteratorType;
	typedef itk::ImageRegionIterator< MaskImageType >
	  ImageIteratorType;
	ImageConstIteratorType itSRC( tempVolume, 
				      tempVolume->GetBufferedRegion() );
	ImageIteratorType itDST( outputMaskVolume, 
				 outputMaskVolume->GetBufferedRegion() );
	while( !itSRC.IsAtEnd() )
	  { itDST.Set( itSRC.Get() );
	    ++itSRC;
	    ++itDST; }
	tempVolume = 0;

	std::cout << "K-means clustering on water-suppressed volume, done." << std::endl;
	std::cout << "Final means are: " << std::endl;
	std::cout << scalarImageKmeansImageFilter->GetFinalMeans() << std::endl;

	this->PostprocessClustering(outputMaskVolume);

      }
      break;
    case KMEANS_3D: // multiparametric K-means.
      {

	// Generate vector image from NS, WS and FS volumes.
	this->FormNSWSFSVectorVolume();

	// 	// Cast to unsigned char vector.
	// 	typedef itk::Vector<MaskPixelType, 3> VectorMaskPixelType;
	// 	typedef itk::Image<VectorMaskPixelType, ImageDimension> VectorMaskImageType;

	// 	typedef itk::CastImageFilter<VectorFloatImageType, VectorMaskImageType> MaskCastImageFilterType;
	// 	MaskCastImageFilterType::Pointer convertToVectorMaskImage = MaskCastImageFilterType::New();
	// 	//		convertToVectorMaskImage->SetInput( formVectorImageFilter->GetOutput() );
	// 	//	convertToVectorMaskImage->Update();

	// Run K-Means.
	MaskImageType::Pointer outputMaskVolume = 
	  KMeansonVectorVolumes();
	this->PostprocessClustering(outputMaskVolume);

      }
      break;
    case REGIONGROWING_3D: // Region growing using vector images.
      {
	// Generate vector image from NS, WS and FS volumes.
	this->FormNSWSFSVectorVolume();

	// Run region growing.
	this->m_outputMaskVolume = RegionGrowingOnVectorImage();
      }
      break;
    case WS_LEVEL_SETS: case WS_GEOD_ACT_CONT:  
      {
	// If the user asked for multiple scale segmentation.
	if (this->m_nScales > 1) 
	  this->MultiscaleLevelSetBasedSegmentation();
	else
	  // Set up ROI->Mask image->Distance Maps
	  // Smoothing->Edge Maps (Parzen)->Speed Image
	  // Set up and run level sets.
	  this->LevelSetBasedSegmentation();
	  
      }
      break;
    case LEVEL_SETS_3D: case GEOD_ACT_CONT_3D:
      {
	// Generate vector image from NS, WS and FS volumes.
	this->FormNSWSFSVectorVolume();

	// If the user asked for multiple scale segmentation.
	if (this->m_nScales > 1) 
	  this->MultiscaleLevelSetBasedSegmentation();
	else
	  // Set up ROI->Mask image->Distance Maps
	  // Smoothing->Edge Maps (Parzen)->Speed Image
	  // Set up and run level sets.
	  this->LevelSetBasedSegmentation();
      }
      break;
    case EM_GMM_BAYES_3D: // Gaussian Mixture Modeling + Bayesian maximum ratio.
      {
	// Generate vector image from NS, WS and FS volumes.
	// this->FormNSWSFSVectorVolume();

	//       typedef itk::BayesianClassifierInitializationImageFilter< VectorFloatImageType > BayesianInitializerType;
	//       BayesianInitializerType::Pointer bayesianInitializer = BayesianInitializerType::New();
      }
      break;
    case KMEANS_MRF_3D: // Multi-parametric K-means followed by MRF.
      {
	// Generate vector image from NS, WS and FS volumes.
	this->FormNSWSFSVectorVolume();

	// Run K-Means.
	MaskImageType::Pointer kmeansOutputMaskImage = 
	  KMeansonVectorVolumes();

	// MRF segmentation.
 	MaskImageType::Pointer outputMaskVolume = 
 	  MRFSegmentationAfterClustering( kmeansOutputMaskImage );

      	this->PostprocessClustering(kmeansOutputMaskImage);
      }
      break;
    }

}



infoStruct
ThymusSegmentationObject::
DisplayVolumetrics()
{

  // Rank components wrt to size and relabel.

  typedef itk::RelabelComponentImageFilter<MaskImageType, 
    MaskImageType> 
    RelabelFilterType;
  RelabelFilterType::Pointer  sortLabelsImageFilter = 
    RelabelFilterType::New();
  sortLabelsImageFilter->SetInput( this->m_outputMaskVolume );
  sortLabelsImageFilter->Update();

//   // Determine the connected component label of thymus.
//   thymusRegionLabel = DetermineMedianROILabel( sortLabelsImageFilter->GetOutput(),
// 					       seedbasedROI);

  unsigned int thymusSizeinPixels = 
    sortLabelsImageFilter->GetSizeOfObjectsInPixels()[ 0 ];
  float thymusSizeinPhysicalUnits = 
    sortLabelsImageFilter->GetSizeOfObjectsInPhysicalUnits()[ 0 ];
  unsigned int nComponents = 
    sortLabelsImageFilter->GetNumberOfObjects();


  // Display volumetric results.

  std::cout << std::endl
	    << "=======Thymus Quantification results======" 
	    << std::endl
	    << "Number of spatially connected regions: " 
	    << nComponents
	    << std::endl
	    << "Thymus component (in pixels) " 
	    << thymusSizeinPixels
	    << std::endl
	    << "Thymus component in physical units (mm^3) "
	    << thymusSizeinPhysicalUnits
	    << std::endl;
  
  infoStruct thymusInfo;
  thymusInfo.SizeinPixels = thymusSizeinPixels;
  thymusInfo.SizeinPhysicalUnits = thymusSizeinPhysicalUnits;

  return( thymusInfo );
}


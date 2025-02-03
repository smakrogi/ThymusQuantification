/*===========================================================================

  Program:   Segmentation/Classification of thymus tissues from MRI.
  Module:    $RCSfile: ThymusSegmentationAlgorithms.cxx,v $
  Language:  C++
  Date:      $Date: 2013/03/04 14:42:00 $
  Version:   $Revision: 1.0 $
  Author:    S. K. Makrogiannis
  3T MRI Facility National Institute on Aging/National Institutes of Health.

  =============================================================================*/

#include "itkObject.h"
#include "itkImage.h"
// #include "itkOrientedImage.h"
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
#include "itkScalarChanAndVeseDenseLevelSetImageFilter.h"
// #include "itkScalarChanAndVeseSparseLevelSetImageFilter.h"
#include "itkScalarChanAndVeseLevelSetFunctionData.h"
#include "itkConstrainedRegionBasedLevelSetFunctionSharedData.h"
#include "itkAtanRegularizedHeavisideStepFunction.h"

// Set-up thymus shape priors for GAC segmentation.
#include "itkGeodesicActiveContourShapePriorLevelSetImageFilter.h"
#include "itkChangeInformationImageFilter.h"
// #include "itkBoundedReciprocalImageFilter.h"
#include "itkPCAShapeSignedDistanceFunction.h"
#include "itkEuler3DTransform.h"
#include "itkOnePlusOneEvolutionaryOptimizer.h"
#include "itkNormalVariateGenerator.h"
#include "vnl/vnl_sample.h"
#include "itkNumericSeriesFileNames.h"
#include "itkSpatialFunctionImageEvaluatorFilter.h"

#include "itkCommand.h"

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
// #include "itkResampleImageFilter.h"
// #include "itkIdentityTransform.h"
// #include "itkLinearInterpolateImageFunction.h"
// #include "itkRecursiveGaussianImageFilter.h"

// Pixel-wise arithmetic.
#include "itkMinimumMaximumImageCalculator.h"
#include "itkMultiplyByConstantImageFilter.h"



// ThymusSegmentation class members.

//! Method for setting the image inputs.
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
  // Check if previous vector volume has been allocated and
  // free memory if needed.
  // if( this->m_VectorVolume  )
  //   this->m_VectorVolume = 0;

  
  // Generate vector image from NS, WS and FS volumes.
  // Instantiate data types for analysis.
  typedef itk::ScalarToArrayCastImageFilter< FloatImageType, 
    VectorFloatImageType > 
    VectorImageGeneratorType;
  VectorImageGeneratorType::Pointer inputVectorImageGenerator =  
    VectorImageGeneratorType::New() ;
  inputVectorImageGenerator->SetInput(0, this->m_smoothedNonSuppressedVolume);
  inputVectorImageGenerator->SetInput(1, this->m_smoothedWaterSuppressedVolume);
  inputVectorImageGenerator->SetInput(2, this->m_smoothedFatSuppressedVolume);
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

  // Post-process values (normalize and complement to 0 if needed).
  for (out.GoToBegin(); !out.IsAtEnd(); ++out)
    {
      double edgeProbability = 255 * (1 - ( out.Get() / totalMax ) );  // for edge map.
      // double edgeProbability = out.Get() / totalMax;  // for speed image.
      out.Set(edgeProbability);
    }

  return edgeMap;
}



//! Generate a probabilistic edge map using Parzen kernels.
//! Works for vector images.
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
  DoubleArrayType localSamples(nSamples, ImageDimension, 0.0);

  double totalMax = 0.0;

  for (it.GoToBegin(), out.GoToBegin(); !it.IsAtEnd(); ++it, ++out)
    {
      // For each voxel of the input volume,
      // generate a vector with the values within this 3D window.
      for(int k=0;k<nSamples;k++) 
	for(int m=0;m<(int)ImageDimension;m++)
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
      // double edgeProbability = out.Get() / totalMax;  // for speed image.
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
    "/" + this->m_subjectID + "_" + SegmentationMethodStrings[this->m_segmentationMethod] + "_" + 
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


//! Over-ride command type to report convergence progress and
//! parameter values.
template<class TFilter>
class CommandIterationUpdate : public itk::Command
{
public:
  typedef CommandIterationUpdate   Self;
  typedef itk::Command             Superclass;
  typedef itk::SmartPointer<Self>  Pointer;
  itkNewMacro( Self );
protected:
  CommandIterationUpdate() {};
public:

  void Execute(itk::Object *caller, const itk::EventObject & event)
  {
    Execute( (const itk::Object *) caller, event);
  }

  void Execute(const itk::Object * object, const itk::EventObject & event)
  {
    const TFilter * filter =
      dynamic_cast< const TFilter * >( object );
    if( typeid( event ) != typeid( itk::IterationEvent ) )
      { return; }
      
    int remainder = filter->GetElapsedIterations() % 50;
    // std::cout << remainder << std::endl;
    if(remainder == 0) {
      std::cout << filter->GetElapsedIterations() << ": ";
      std::cout << filter->GetRMSChange() << " ";
      std::cout << filter->GetCurrentParameters() << std::endl;
    }
  }

};


//! Geodesic active contours with shape prior model.
FloatImageType::Pointer
ThymusSegmentationObject::
GeodesicActiveContoursWithShapePriors(FloatImageType::Pointer ROIDistanceOutput,
				      FloatImageType::Pointer speedImage)
{
  //! Set path to shape model files.
  char* tempString  = NULL;
  std::string shapemodelPath ;
  tempString = getenv(SHAPE_MODEL_PATH_ENV_VAR);
  if ( tempString == NULL ) {
    std::cerr << "Can not find shape model path" 
	      << std::endl;
    exit(1);
  }
  else
    shapemodelPath = tempString;
 
  //! GAC with shape priors filter.
  typedef  itk::GeodesicActiveContourShapePriorLevelSetImageFilter<
  FloatImageType,
    FloatImageType >   GeodesicActiveContourWithShapePriorFilterType;
  GeodesicActiveContourWithShapePriorFilterType::Pointer geodesicActiveContourWithShapePrior =
    GeodesicActiveContourWithShapePriorFilterType::New();

  // // The \doxygen{ChangeInformationImageFilter} is the first filter in the preprocessing
  // // stage and is used to force the image origin to the center of the image.
  // //
  // //  \index{itk::ChangeInformationImageFilter!CenterImageOn()}

  typedef itk::ChangeInformationImageFilter<
  FloatImageType >  CenterFilterType;

  CenterFilterType::Pointer center = CenterFilterType::New();
  center->CenterImageOn();
  center->SetInput( ROIDistanceOutput );
  center->Update();

  CenterFilterType::Pointer center2 = CenterFilterType::New();
  center2->CenterImageOn();
  center2->SetInput( speedImage );
  center2->Update();

  geodesicActiveContourWithShapePrior->SetPropagationScaling(  this->m_levelsetPropagationScalingFactor );
  geodesicActiveContourWithShapePrior->SetShapePriorScaling( this->m_levelsetShapePriorScalingFactor );
  geodesicActiveContourWithShapePrior->SetCurvatureScaling( this->m_levelsetCurvatureScalingFactor ); // 1.0 
  geodesicActiveContourWithShapePrior->SetAdvectionScaling( this->m_levelsetAdvectionScalingFactor ); // 1.0 
  geodesicActiveContourWithShapePrior->SetMaximumRMSError( this->m_levelsetMaximumRMSError ); // 0.005
  geodesicActiveContourWithShapePrior->SetNumberOfIterations( (int) this->m_levelsetMaximumIterations ); // 400
  geodesicActiveContourWithShapePrior->SetIsoSurfaceValue( 0.0 );
  geodesicActiveContourWithShapePrior->SetInput( center->GetOutput() ); // ROIDistanceOutput
  geodesicActiveContourWithShapePrior->SetFeatureImage( center2->GetOutput() ); // speedImage // may use other speed volume modified by the distance.
  // geodesicActiveContourWithShapePrior->UseImageSpacingOn();

  //  Each iteration, the current ``best-fit'' shape is estimated from the
  //  edge potential image and the current contour. To increase speed, only
  //  information within the sparse field layers of the current contour is used
  //  in the estimation. The default number of sparse field layers is
  //  the same as
  //  the ImageDimension which does not contain enough information to get
  //  a reliable best-fit shape estimate. Thus, we override the default and

  geodesicActiveContourWithShapePrior->SetNumberOfLayers( 6 ); // 4,6,8,12


  //  Next, we define the shape model. In this example,
  //  we use an implicit shape model based on the principal components
  //  such that:
  //
  //  \begin{equation}
  //  \psi^{*}(\mathbf{x}) = \mu(\mathbf{x}) + \sum_k \alpha_k u_k(\mathbf{x})
  //  \end{equation}
  //
  //  where $\mu(\mathbf{x})$ is the mean signed distance computed from training
  //  set of segmented objects and $u_k(\mathbf{x})$ are the first $K$ principal
  //  components of the offset (signed distance - mean).
  //  The coefficients $\{\alpha_k\}$ form the
  //  set of \emph{shape} parameters.
  //
  //  Given a set of training data, the \doxygen{ImagePCAShapeModelEstimator}
  //  can be used to obtain
  //  the mean and principal mode shape images required by PCAShapeSignedDistanceFunction.
  //
  //  \index{itk::PCAShapeSignedDistanceFunction!New()}
  //  \index{itk::PCAShapeSignedDistanceFunction!SetNumberOfPrincipalComponents()}

  const unsigned int numberOfPCAModes =  this->m_levelsetNumberOfPCAModes;

  // Software Guide : BeginCodeSnippet
  typedef itk::PCAShapeSignedDistanceFunction<
  double,
    ImageDimension,
    FloatImageType >     ShapeFunctionType;

  ShapeFunctionType::Pointer shape = ShapeFunctionType::New();

  shape->SetNumberOfPrincipalComponents( this->m_levelsetNumberOfPCAModes );

  //  In this example, we will read the mean shape and
  //  principal mode images from file. We will assume that
  //  the filenames of the mode images form a numeric series starting from index 0.
  //
  //  \index{itk::PCAShapeSignedDistanceFunction!SetMeanImage()}
  //  \index{itk::PCAShapeSignedDistanceFunction!SetPrincipalComponentsImages()}
  //
  typedef  itk::ImageFileReader< FloatImageType > ReaderType;
  ReaderType::Pointer meanShapeReader = ReaderType::New();
  meanShapeReader->SetFileName( shapemodelPath+DefaultMeanShapeFilename );
  std::cout << "Reading " << DefaultMeanShapeFilename << std::endl;
  meanShapeReader->Update();
  
  CenterFilterType::Pointer center3 = CenterFilterType::New();
  center3->CenterImageOn();
  center3->SetInput( meanShapeReader->GetOutput() );
  center3->Update();


  std::vector<FloatImageType::Pointer> shapeModeImages( this->m_levelsetNumberOfPCAModes );

  itk::NumericSeriesFileNames::Pointer fileNamesCreator =
    itk::NumericSeriesFileNames::New();

  fileNamesCreator->SetStartIndex( 0 );
  fileNamesCreator->SetEndIndex( this->m_levelsetNumberOfPCAModes - 1 );
  fileNamesCreator->SetSeriesFormat( DefaultShapeModeFilenamePattern );
  const std::vector<std::string> & shapeModeFileNames =
    fileNamesCreator->GetFileNames();

  std::cout << "Reading " << std::endl;
  for ( unsigned int k = 0; k < numberOfPCAModes; k++ )
    {
      ReaderType::Pointer shapeModeReader = ReaderType::New();
      shapeModeReader->SetFileName( shapemodelPath+shapeModeFileNames[k] );
      std::cout << shapeModeFileNames[k] << std::endl;
      shapeModeReader->Update();

      CenterFilterType::Pointer center4 = CenterFilterType::New();
      center4->CenterImageOn();
      center4->SetInput( shapeModeReader->GetOutput() );
      center4->Update();

      shapeModeImages[k] = center4->GetOutput(); // shapeModeReader->GetOutput()
    }

  shape->SetMeanImage( center3->GetOutput() ); // meanShapeReader->GetOutput()
  shape->SetPrincipalComponentImages( shapeModeImages );

  // Further we assume that the shape modes have been normalized
  // by multiplying with the corresponding singular value. Hence,
  // we can set the principal component standard deviations to all
  // ones.
  //
  //  \index{itk::PCAShapeSignedDistanceFunction!Set\-Principal\-Component\-Standard\-Deviations()}
 
  ShapeFunctionType::ParametersType pcaStandardDeviations( this->m_levelsetNumberOfPCAModes );
  pcaStandardDeviations.Fill( 1.0 );

  // Use following code to read-in eigenvalues if PCs 
  // have not been normalized before.
  // std::ifstream infile("eigenvalues.txt");
  // std::string line;
  // double value;
  // unsigned int k = 0;

  // if (std::getline(infile, line))
  //   {
  //     std::istringstream iss(line);
  //     while (iss >> value) { /* process first line */ }
  //     pcaStandardDeviations[k] = value;
  //     k++;
  //   }

  std::cout << "PC standard deviations: "
	    << pcaStandardDeviations
	    << std::endl;

  shape->SetPrincipalComponentStandardDeviations(pcaStandardDeviations);


  // Next, we instantiate a \doxygen{Euler3DTransform} and connect it to the
  // PCASignedDistanceFunction. The transform represent
  // the pose of the shape. The parameters of the transform
  // forms the set of \emph{pose} parameters.
  //
  //  \index{itk::PCAShapeSignedDistanceFunction!SetTransform()}
  //  \index{itk::ShapeSignedDistanceFunction!SetTransform()}

  typedef itk::Euler3DTransform<double>    TransformType;
  TransformType::Pointer transform = TransformType::New();

  shape->SetTransform( transform );

  // Before updating the level set at each iteration, the parameters
  // of the current best-fit shape is estimated by minimizing the
  // \doxygen{ShapePriorMAPCostFunction}. The cost function is composed of
  // four terms: contour fit, image fit, shape prior and pose prior.
  // The user can specify the weights applied to each term.
  //
  //  \index{itk::ShapePriorMAPCostFunction!SetWeights()}

  typedef itk::ShapePriorMAPCostFunction<
  FloatImageType,
    FloatPixelType >     CostFunctionType;

  CostFunctionType::Pointer costFunction = CostFunctionType::New();

  CostFunctionType::WeightsType weights;
  weights[0] =  this->m_levelsetMAPContourWeight; // 1.0;  // weight for contour fit term
  weights[1] =  this->m_levelsetMAPImageGradientWeight; // 10.0; // weight for image fit term  // 20.0
  weights[2] =  this->m_levelsetMAPShapeWeight; // 1.0;  // weight for shape prior term
  weights[3] =  1.0;  // weight for pose prior term
  costFunction->SetWeights(weights);


  // Contour fit measures the likelihood of seeing the current
  // evolving contour for a given set of shape/pose parameters.
  // This is computed by counting the number of pixels inside
  // the current contour but outside the current shape.
  //
  // Image fit measures the likelihood of seeing certain image
  // features for a given set of shape/pose parameters. This is
  // computed by assuming that ( 1 - edge potential ) approximates
  // a zero-mean, unit variance Gaussian along the normal of
  // the evolving contour. Image fit is then computed by computing
  // the Laplacian goodness of fit of the Gaussian:
  //
  // \begin{equation}
  // \sum \left( G(\psi(\mathbf{x})) - |1 - g(\mathbf{x})| \right)^2
  // \end{equation}
  //
  // where $G$ is a zero-mean, unit variance Gaussian and $g$ is
  // the edge potential feature image.
  //
  // The pose parameters are assumed to have a uniform distribution
  // and hence do not contribute to the cost function.
  // The shape parameters are assumed to have a Gaussian distribution.
  // The parameters of the distribution are user-specified. Since we
  // assumed the principal modes have already been normalized,
  // we set the distribution to zero mean and unit variance.
  //
  //  \index{itk::ShapePriorMAPCostFunction!SetShapeParameterMeans()}
  //  \index{itk::ShapePriorMAPCostFunction!SetShapeParameterStandardDeviations()}

  CostFunctionType::ArrayType mean( shape->GetNumberOfShapeParameters() );
  CostFunctionType::ArrayType stddev( shape->GetNumberOfShapeParameters() );

  mean.Fill( 0.0 );
  stddev.Fill( 1.0 );
  costFunction->SetShapeParameterMeans( mean );
  costFunction->SetShapeParameterStandardDeviations( stddev );

  // In this example, we will use the \doxygen{OnePlusOneEvolutionaryOptimizer}
  // to optimize the cost function.

  typedef itk::OnePlusOneEvolutionaryOptimizer    OptimizerType;
  OptimizerType::Pointer optimizer = OptimizerType::New();

  // The evolutionary optimization algorithm is based on testing
  // random permutations of the parameters. As such, we need to provide
  // the optimizer with a random number generator. In the following lines,
  // we create a \doxygen{NormalVariateGenerator}, seed it, and
  // connect it to the optimizer.
  //
  //  \index{itk::Statistics::NormalVariateGenerator!Initialize()}
  //  \index{itk::OnePlusOneEvolutionaryOptimizer!SetNormalVariateGenerator()}

  typedef itk::Statistics::NormalVariateGenerator GeneratorType;
  GeneratorType::Pointer generator = GeneratorType::New();

  generator->Initialize( 20020702 );

  optimizer->SetNormalVariateGenerator( generator );

  // The cost function has $K+3$ parameters. The first $K$
  // parameters are the principal component multipliers, followed
  // by three 3D rotation parameters (in radians) and the (x,y,z)
  // translation parameters (in mm).  We need to carefully
  // scale the different types of parameters to compensate
  // for the differences in the dynamic ranges of the parameters.
  //
  //  \index{itk::OnePlusOneEvolutionaryOptimizer!SetScales()}
  //
  OptimizerType::ScalesType scales( shape->GetNumberOfParameters() );
  scales.Fill( 1.0 );
  for( unsigned int k = 0; k < (unsigned int)this->m_levelsetNumberOfPCAModes; k++ )
    {
      scales[k] = 20.0;  // scales for the pca mode multiplier
    }
  scales[this->m_levelsetNumberOfPCAModes] = 350.0;  // scale for 3D rotation
  scales[this->m_levelsetNumberOfPCAModes+1] = 350.0;  // scale for 3D rotation
  scales[this->m_levelsetNumberOfPCAModes+2] = 350.0;  // scale for 3D rotation

  optimizer->SetScales( scales );

  // Next, we specify the initial radius, the shrink and
  // grow mutation factors and termination criteria of the optimizer.
  // Since the best-fit shape is re-estimated each iteration of
  // the curve evolution, we do not need to spend too much time finding the true
  // minimizing solution each time; we only need to head towards it. As such,
  // we only require a small number of optimizer iterations.
  //
  //  \index{itk::OnePlusOneEvolutionaryOptimizer!Initialize()}
  //  \index{itk::OnePlusOneEvolutionaryOptimizer!SetEpsilon()}
  //  \index{itk::OnePlusOneEvolutionaryOptimizer!SetMaximumIteration()}

  double initRadius = 1.05;
  double grow = 1.1;
  double shrink = pow(grow, -0.25);
  optimizer->Initialize(initRadius, grow, shrink);

  optimizer->SetEpsilon(1.0e-6); // minimal search radius

  optimizer->SetMaximumIteration(15);

  // Before starting the segmentation process we need to also supply the initial
  // best-fit shape estimate. In this example, we start with the unrotated mean shape
  // with the initial x- and y- translation specified through command-line
  // arguments.

  ShapeFunctionType::ParametersType parameters( shape->GetNumberOfParameters() );
  parameters.Fill( 0.0 );
  // parameters[numberOfPCAModes + 1] = atof( argv[16] ); // startX
  // parameters[numberOfPCAModes + 2] = atof( argv[17] ); // startY

  // Finally, we connect all the components to the filter and add our
  // observer.

  geodesicActiveContourWithShapePrior->SetShapeFunction( shape );
  geodesicActiveContourWithShapePrior->SetCostFunction( costFunction );
  geodesicActiveContourWithShapePrior->SetOptimizer( optimizer );
  geodesicActiveContourWithShapePrior->SetInitialParameters( parameters );

  //! Add observer through our custom command type.
  typedef CommandIterationUpdate<GeodesicActiveContourWithShapePriorFilterType> CommandType;
  CommandType::Pointer observer = CommandType::New();
  geodesicActiveContourWithShapePrior->AddObserver( itk::IterationEvent(), observer );

  //  The invocation of the \code{Update()} method triggers the
  //  execution of the pipeline.  As usual, the call is placed in a
  //  \code{try/catch} block to handle exceptions should errors occur.
  try{
    geodesicActiveContourWithShapePrior->Update();
  }
  catch( itk::ExceptionObject & excep )
    {
      std::cerr << "Exception caught !" << std::endl;
      std::cerr << excep << std::endl;
    }


  FloatImageType::Pointer levelsetOutput = 0;
  levelsetOutput =
    geodesicActiveContourWithShapePrior->GetOutput();
  //! Set origin to original.
  levelsetOutput->SetOrigin(this->m_watersuppressedVolume->GetOrigin());


  // Also write out the initial and final best fit shape
  typedef itk::SpatialFunctionImageEvaluatorFilter<
  ShapeFunctionType,
    FloatImageType,
    FloatImageType >  EvaluatorFilterType;

  EvaluatorFilterType::Pointer evaluator = EvaluatorFilterType::New();
  evaluator->SetInput( levelsetOutput );
  evaluator->SetFunction( shape );
  shape->SetParameters( geodesicActiveContourWithShapePrior->GetInitialParameters() );

  typedef itk::BinaryThresholdImageFilter<
  FloatImageType,
    MaskImageType    >       ThresholdingFilterType;

  ThresholdingFilterType::Pointer thresholder = ThresholdingFilterType::New();
  
  thresholder->SetLowerThreshold( -1000.0 );
  thresholder->SetUpperThreshold(     0.0 );

  thresholder->SetOutsideValue(  BACKGROUND  );
  thresholder->SetInsideValue(  FOREGROUND );

  thresholder->SetInput( evaluator->GetOutput() );
  
  // Write intermediate images for debugging purposes.
  typedef  itk::ImageFileWriter<  MaskImageType  > WriterType;
  WriterType::Pointer writer = WriterType::New();
  writer->SetInput(thresholder->GetOutput() );
  writer->SetFileName( "GACShapePriorInitialShape.mha" );
  writer->Update();

  shape->SetParameters( geodesicActiveContourWithShapePrior->GetCurrentParameters() );
  evaluator->Modified();
  writer->SetFileName( "GACShapePriorFinalShape.mha" );
  writer->Update();

  typedef  itk::ImageFileWriter<  FloatImageType  > FloatWriterType;
  FloatWriterType::Pointer writer2 = FloatWriterType::New();
  writer2->SetInput(center2->GetOutput());
  writer2->SetFileName( "GACShapePriorSpeedImage.mha" );
  writer2->Update();
  
  std::cout << "Geodesic active contours with shape prior." 
	    << std::endl
	    << "Propagation factor = "
	    << this->m_levelsetPropagationScalingFactor
	    << ", Curvature factor = "
	    << this->m_levelsetCurvatureScalingFactor
	    << ", Advection factor = "
	    << this->m_levelsetAdvectionScalingFactor
	    << ", Shape prior factor = "
	    << this->m_levelsetShapePriorScalingFactor
	    << ", MAP Contour weight = "
	    << this->m_levelsetMAPContourWeight
	    << ", MAP Im. Grad. weight = "
	    << this->m_levelsetMAPImageGradientWeight
	    << ". MAP Shape weight = "
	    << this->m_levelsetMAPShapeWeight
	    << std::endl;
  // geodesicActiveContourWithShapePrior->Print(std::cout);

  // Print out some useful information.
  std::cout << "Level set segmentation of thymus completed." << std::endl;
  std::cout << "Max. no. iterations: " << geodesicActiveContourWithShapePrior->GetNumberOfIterations() << "," <<
    "Max. RMS error: " << geodesicActiveContourWithShapePrior->GetMaximumRMSError() << std::endl;
  std::cout << "# of iterations: " << geodesicActiveContourWithShapePrior->GetElapsedIterations() << ","
	    << "RMS change: " << geodesicActiveContourWithShapePrior->GetRMSChange() << std::endl;

  // // For debugging purposes write out the level set output.
  // FloatWriterType::Pointer writer3 = FloatWriterType::New();
  // writer3->SetInput(levelsetOutput);
  // writer3->SetFileName( "GACShapePriorLevelSetOutput.mha" );
  // writer3->Update();

  return levelsetOutput;
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

  FloatImageType::Pointer gradientImage = 0; 

  typedef itk::ImageRegionConstIterator< GradientMagnitudeFilterType::OutputImageType >
    VectorImageConstIteratorType;
  typedef itk::ImageRegionIterator< FloatImageType >
    FloatImageIteratorType;

  // Choose gradient method based on the image intensity dimensionality and
  // the desired algorithm (gradient magnitude or stochastic method).
  switch(this->m_segmentationMethod) {
    // Scalar pixel intesities.
  case WS_LEV_SETS: case WS_GAC: case WS_CHAN_VESE: case WS_GAC_SHP:
    {
      switch(this->m_gradientMethod) {
      case CONVENTIONAL:
	{ 
	  // Compute the gradient magnitude from the WS image only.
	  GradientFilterType::Pointer  gradientMagnitude = GradientFilterType::New();
	  gradientMagnitude->SetSigma( this->m_gradientSigma ); // 0.0625
	  gradientMagnitude->SetInput( this->m_smoothedWaterSuppressedVolume );
	  gradientMagnitude->Update();
	  gradientImage = gradientMagnitude->GetOutput();
	}
	break;
      case PARZEN:
	{
	  // Edge map using Parzen estimators.
	  gradientImage = ParzenEdgeMap( this->m_smoothedWaterSuppressedVolume );
	}
	break;
      }
    }
    break;
    // Vector pixel intensities.
  case LEV_SETS_3D: case GAC_3D: case GAC_SHP_3D:
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
    this->m_subjectID + "_" + GradientMethodStrings[this->m_gradientMethod] + "_" + 
    SegmentationMethodStrings[this->m_segmentationMethod] + "_" + 
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
  typedef itk::SignedMaurerDistanceMapImageFilter<MaskImageType, 
    FloatImageType> 
    LabelDistanceTransformType;
  LabelDistanceTransformType::Pointer ROIDistFilter = 
    LabelDistanceTransformType::New();
  ROIDistFilter->SetInput( roiVolume );
  ROIDistFilter->SetUseImageSpacing( true );

  // In chan-vese filter, region interior is 
  // denoted by positive values.
  if ( this->m_segmentationMethod == WS_CHAN_VESE )
    ROIDistFilter->SetInsideIsPositive( true );

  ROIDistFilter->SetSquaredDistance( true );  // false
  ROIDistFilter->Update();
  ROIDistanceOutput = ROIDistFilter->GetOutput();
  // ROIDistFilter = 0;

  // // Instantiate caster because level set filters don't run on oriented images.
  // typedef itk::CastImageFilter< itk::Image<float, 3>, 
  //   FloatImageType > 
  //   CastFilterType;
  // CastFilterType::Pointer castImage = 
  //   CastFilterType::New();
  FloatImageType::Pointer levelSetOutputImage = 
    FloatImageType::New();

  // Run segmentation using shape detection or geodesic active contours.
  switch(this->m_segmentationMethod) {
  case WS_LEV_SETS: case LEV_SETS_3D:
    {
      typedef  itk::ShapeDetectionLevelSetImageFilter< FloatImageType, 
	FloatImageType >    
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
  case WS_GAC: case GAC_3D:
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
  case WS_GAC_SHP: case GAC_SHP_3D:

    // Call function for GAC with shape prior.
    levelSetOutputImage = this->GeodesicActiveContoursWithShapePriors( ROIDistanceOutput,
								       sigmoid->GetOutput() );
    break;
  case WS_CHAN_VESE:
    {

      // First, check the minimum and maximum intensities.
      // If min=0 and max=1, multiply pixel intensities with 255.
      typedef itk::MinimumMaximumImageCalculator <FloatImageType>
	ImageCalculatorFilterType;
      ImageCalculatorFilterType::Pointer imageCalculatorFilter
	= ImageCalculatorFilterType::New ();
      imageCalculatorFilter->SetImage(this->m_smoothedWaterSuppressedVolume);
      imageCalculatorFilter->Compute();
      float intensityMinimum  = imageCalculatorFilter->GetMinimum();
      float intensityMaximum  = imageCalculatorFilter->GetMaximum();

      // Check if we are processing the fat ratio image.
      // In that case multiply by 2^6-1=63 to facilitate chan-vese segmentation.
      FloatImageType::Pointer rescaledSmoothedWaterSuppressedVolume = 0;
      typedef itk::MultiplyByConstantImageFilter <FloatImageType, FloatPixelType, 
	FloatImageType> MultiplyByConstantImageFilterType;
      if ( (intensityMaximum<=1.0) & (intensityMinimum>=0) ) {
	// Multiply value by constant.
	float scalingConstant = (float)63.0;
	std::cout << "Multiplying ratio values by "
		  << scalingConstant
		  << " before Chan-Vese execution."
		  << std::endl; 
	MultiplyByConstantImageFilterType::Pointer 
	  multiplierFilter = MultiplyByConstantImageFilterType::New();
	multiplierFilter->SetInput( this->m_smoothedWaterSuppressedVolume );
	multiplierFilter->SetConstant( scalingConstant );
	multiplierFilter->Update();
	rescaledSmoothedWaterSuppressedVolume = multiplierFilter->GetOutput();
      } 
      else
	rescaledSmoothedWaterSuppressedVolume = this->m_watersuppressedVolume;
     

      // Algorithm parameters.
      double epsilon = 1.;
      // double reinitialization_weight = 0.;
      double l1 = 1.;
      double l2 = 1.;
      //      double curvature_weight = 0.;
      //      double area_weight = 0.;
      //      double volume_weight = 0.;
      //      double volume = 0.;
 
      // Typedef and instantiate filter and helper objects.
      //
      //  We now define the image type using a particular pixel type and
      //  dimension. In this case the \code{float} type is used for the pixels
      //  due to the requirements of the smoothing filter.
      typedef itk::ScalarChanAndVeseLevelSetFunctionData< FloatImageType,
	FloatImageType > DataHelperType;
      
      typedef itk::ConstrainedRegionBasedLevelSetFunctionSharedData<
	FloatImageType, FloatImageType, DataHelperType > SharedDataHelperType;
 
      // typedef itk::ScalarChanAndVeseLevelSetFunction< FloatImageType,
      // 	FloatImageType, SharedDataHelperType > LevelSetFunctionType;

      typedef itk::ScalarChanAndVeseLevelSetFunction< FloatImageType,
	FloatImageType, SharedDataHelperType > LevelSetFunctionType;

      typedef itk::ScalarChanAndVeseDenseLevelSetImageFilter< FloatImageType,
	FloatImageType, FloatImageType, LevelSetFunctionType,
	SharedDataHelperType > MultiLevelSetType;

      // typedef itk::ScalarChanAndVeseSparseLevelSetImageFilter< FloatImageType,
      // 	FloatImageType, FloatImageType, LevelSetFunctionType,
      // 	SharedDataHelperType > MultiLevelSetType;


      //  We declare now the type of the numerically discretized Step and Delta functions that
      //  will be used in the level-set computations for foreground and background regions
      //
      typedef itk::AtanRegularizedHeavisideStepFunction< FloatPixelType,
	FloatPixelType >  DomainFunctionType;
      
      DomainFunctionType::Pointer domainFunction = 
	DomainFunctionType::New();
      domainFunction->SetEpsilon( epsilon );

      MultiLevelSetType::Pointer chanandveseLevelSetFilter = 
	MultiLevelSetType::New();

      // First set the number of level-sets that evolve.
      chanandveseLevelSetFilter->SetFunctionCount( 1 );
 
      //  Set the feature image and initial level-set image as output of the
      //  fast marching image filter.
      chanandveseLevelSetFilter->SetFeatureImage( rescaledSmoothedWaterSuppressedVolume );
      chanandveseLevelSetFilter->SetLevelSet( 0, ROIDistanceOutput );
      chanandveseLevelSetFilter->SetNumberOfIterations( (int) this->m_levelsetMaximumIterations);
      chanandveseLevelSetFilter->SetMaximumRMSError( this->m_levelsetMaximumRMSError );
      chanandveseLevelSetFilter->SetUseImageSpacing( 1 );
      chanandveseLevelSetFilter->SetInPlace( false );

      //  For the level set with phase 0, set different parameters and weights. This may
      //  to be set in a loop for the case of multiple level-sets evolving simultaneously.
      //
      chanandveseLevelSetFilter->
	GetDifferenceFunction(0)->
	SetDomainFunction( domainFunction );
      chanandveseLevelSetFilter->
	GetDifferenceFunction(0)->
	SetCurvatureWeight( this->m_chanveseCurvatureWeight );
      chanandveseLevelSetFilter->
	GetDifferenceFunction(0)->
	SetAreaWeight( this->m_chanveseAreaWeight );
      // chanandveseLevelSetFilter->
      // 	GetDifferenceFunction(0)->
      // 	SetReinitializationSmoothingWeight( reinitialization_weight );
      chanandveseLevelSetFilter->GetDifferenceFunction(0)->
	SetVolumeMatchingWeight( this->m_chanveseVolumeWeight );
      chanandveseLevelSetFilter->
	GetDifferenceFunction(0)->
	SetVolume( this->m_chanveseVolume );
      chanandveseLevelSetFilter->
	GetDifferenceFunction(0)->
	SetLambda1( l1 );
      chanandveseLevelSetFilter->
	GetDifferenceFunction(0)->
	SetLambda2( l2 );
      chanandveseLevelSetFilter->
	GetDifferenceFunction(0)->
	SetOverlapPenaltyWeight(this->m_chanveseOverlapWeight);

      chanandveseLevelSetFilter->Update();

      std::cout << "Chan-Vese level sets." 
		<< std::endl
		<< "curvature weight = "
		<< this->m_chanveseCurvatureWeight
		<< ", area weight  = "
		<< this->m_chanveseAreaWeight
		<< ", overlap weight  = "
		<< this->m_chanveseOverlapWeight
		<< std::endl;

      levelSetOutputImage = 
	chanandveseLevelSetFilter->GetOutput();

      // Print out some useful information.
      std::cout << "Level set segmentation of thymus completed." << std::endl;
      std::cout << "Max. no. iterations: " << chanandveseLevelSetFilter->GetNumberOfIterations() << "," <<
	"Max. RMS error: " << chanandveseLevelSetFilter->GetMaximumRMSError() << std::endl;
      std::cout << "# of iterations: " << chanandveseLevelSetFilter->GetElapsedIterations() << ","
		<< "RMS change: " << chanandveseLevelSetFilter->GetRMSChange() << std::endl;
 
    }
  }

  // // Type casting before thresholding.
  // castImage->SetInput( levelSetOutputImage );
  // castImage->Update();

  // Apply thresholding to the level sets
  // to get the segmentation result.
  typedef itk::BinaryThresholdImageFilter<FloatImageType, MaskImageType> 
    ThresholdingFilterType;
  ThresholdingFilterType::Pointer thresholder = ThresholdingFilterType::New();
  thresholder->SetInput( levelSetOutputImage ); // previously: castImage->GetOutput() );
  thresholder->SetOutsideValue(  BACKGROUND  );
  thresholder->SetInsideValue(  FOREGROUND );
  thresholder->SetLowerThreshold( -1000.0 );
  thresholder->SetUpperThreshold(     0.0 );
  thresholder->Update();
  // Free possibly already allocated volume.
  this->m_outputMaskVolume = 0;
  this->m_outputMaskVolume = thresholder->GetOutput();
  this->m_outputMaskVolume->SetOrigin(this->m_watersuppressedVolume->GetOrigin());

  std::cout << "Level set-based segmentation, done." 
	    << std::endl;
}


// Run level sets in multiple scales generated in the Parzen domain and
// the scale-space derived from gaussian, diffusion and median filtering.
void 
ThymusSegmentationObject::
MultiscaleLevelSetBasedSegmentation() {
  //! Loop over the number of scales.
  std::cout << "Using multi-scale level set-based segmentation scheme." 
	    << std::endl;

  //! Record the original parameters.
  double originalParzenBandwidth = this->m_ParzenBandwidth;
  // int originalLevelSetMaximumIterations = this->m_levelsetMaximumIterations;
  // float originalGradientSigma = this->m_gradientSigma;
  const double parzenbandwidthBase = (double) this->m_ParzenBandwidthBase;
  // const float gradientsigmaBase = 2.0;

  //! Iterate segmentation over several scales.
  for( this->m_currentScale = 0;
       this->m_currentScale < this->m_nScales;
       this->m_currentScale++ ) {
    
    //! Switch to selected multi-scale generator.
    switch(this->m_multiscaleGenerator) {
    case PARZEN_EDGES:
      {
	//! If multiple scales are generated by Parzen bandwidth.
	//! Set the Parzen variables.
	this->m_ParzenBandwidth = 
	  originalParzenBandwidth * 
	  (double) std::pow(parzenbandwidthBase, (double)( (-1)*this->m_currentScale) );


	// //! Set the level-set iteration number.
	// this->m_levelsetMaximumIterations = 
	//   originalLevelSetMaximumIterations *
	//   std::pow(2, (this->m_currentScale - this->m_nScales + 1) );
      
	break;
      case GAUSSIAN: case DIFFUSION: case MEDIAN:
	{
	  // //! Set gradient magnitude sigma (works for scalar inputs only). 
	  // this->m_gradientSigma = 
	  //   originalGradientSigma *
	  //   (float)std::pow(gradientsigmaBase, (float)( (-1)*this->m_currentScale) );

	  //! Set as starting point the original images.
	  this->m_smoothedNonSuppressedVolume = 0;
	  this->m_smoothedNonSuppressedVolume = this->m_nonsuppressedVolume;
	  this->m_smoothedWaterSuppressedVolume = 0;
	  this->m_smoothedWaterSuppressedVolume = this->m_watersuppressedVolume;
	  this->m_smoothedFatSuppressedVolume = 0;
	  this->m_smoothedFatSuppressedVolume = this->m_fatsuppressedVolume;

	  //! Generate multiple scales and run filtering.
	  int scalespaceIterations = 
	    (int) MY_ROUND(std::pow((float)2.0, (float)(this->m_nScales - this->m_currentScale - 1) ) );
	
	  for (int i=0; 
	       i<scalespaceIterations; 
	       i++) {

	    // Choose smoothing based on the image intensity dimensionality and
	    // the desired algorithm (gradient magnitude or stochastic method).
	    switch(this->m_segmentationMethod) {
	      // Scalar pixel intesities.
	    case WS_LEV_SETS: case WS_GAC: case WS_CHAN_VESE:
	      {
		// Non-linear smoothing.
		this->m_smoothedWaterSuppressedVolume = 
		  SmoothInputVolume( this->m_smoothedWaterSuppressedVolume, 
				     this->m_multiscaleGenerator );
	      }
	      break;
	    case LEV_SETS_3D: case GAC_3D:
	      {
		// Smooth out all three input volumes and
		// generate the vector volume.
		this->m_smoothedNonSuppressedVolume = 
		  SmoothInputVolume( this->m_smoothedNonSuppressedVolume, 
				     this->m_multiscaleGenerator );

		this->m_smoothedWaterSuppressedVolume = 
		  SmoothInputVolume( this->m_smoothedWaterSuppressedVolume, 
				     this->m_multiscaleGenerator );

		this->m_smoothedFatSuppressedVolume = 
		  SmoothInputVolume( this->m_smoothedFatSuppressedVolume, 
				     this->m_multiscaleGenerator );

		// Construct the vector volume from the three components.
		this->FormNSWSFSVectorVolume();
	      }
	    }
	  }
	  break;
	}
      }
      break;
    }

    //! Run level sets.
    this->LevelSetBasedSegmentation();
  }
}



// Gaussian Mixture Modeling + Bayesian maximum ratio.

// MaskImageType::Pointer 
// ThymusSegmentationObject::
// GMM_BayesianClassification()
// {
// // Set up the segmentation process.
// MaskImageType::Pointer outputMaskVolume;

// std::cout << "GMM using E-M optimization and Bayesian clustering on NS, WS and FS volumes, done." 
// 	    << std::endl;

// return outputMaskVolume;
// }



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
PostprocessClustering(MaskImageType::Pointer outputMaskVolume)
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
    "/" + this->m_subjectID + "_" + 
    SegmentationMethodStrings[this->m_segmentationMethod] + "_" + 
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
    this->m_subjectID + "_" + 
    SegmentationMethodStrings[this->m_segmentationMethod] + "_" + 
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
	this->m_smoothedWaterSuppressedVolume = 
	  SmoothInputVolume( this->m_watersuppressedVolume, MEDIAN );

	// Clustering (fcm or k-means).
	// Use 2 classes.
	typedef itk::ScalarImageKmeansImageFilter<FloatImageType> 
	  ScalarImageKmeansImageFilterType;

	ScalarImageKmeansImageFilterType::Pointer scalarImageKmeansImageFilter = 
	  ScalarImageKmeansImageFilterType::New();

	scalarImageKmeansImageFilter->SetInput( this->m_smoothedWaterSuppressedVolume );
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
	// Smooth out input volumes.
	this->m_smoothedNonSuppressedVolume = 
	  SmoothInputVolume( this->m_nonsuppressedVolume, 
			     MEDIAN );
	this->m_smoothedWaterSuppressedVolume = 
	  SmoothInputVolume( this->m_watersuppressedVolume, 
			     MEDIAN );
	this->m_smoothedFatSuppressedVolume = 
	  SmoothInputVolume( this->m_fatsuppressedVolume, 
			     MEDIAN );

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
    case REG_GROW_3D: // Region growing using vector images.
      {

	// Smooth out input volumes.
	this->m_smoothedNonSuppressedVolume = 
	  SmoothInputVolume( this->m_nonsuppressedVolume, 
			     MEDIAN );
	this->m_smoothedWaterSuppressedVolume = 
	  SmoothInputVolume( this->m_watersuppressedVolume, 
			     MEDIAN );
	this->m_smoothedFatSuppressedVolume = 
	  SmoothInputVolume( this->m_fatsuppressedVolume, 
			     MEDIAN );

	// Generate vector image from NS, WS and FS volumes.
	this->FormNSWSFSVectorVolume();

	// Run region growing.
	this->m_outputMaskVolume = RegionGrowingOnVectorImage();
      }
      break;
    case WS_LEV_SETS: case WS_GAC: case WS_GAC_SHP: case WS_CHAN_VESE: 
      {
	// Non-linear smoothing.
	this->m_smoothedWaterSuppressedVolume = 
	  SmoothInputVolume( this->m_watersuppressedVolume, MEDIAN );

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
    case LEV_SETS_3D: case GAC_3D: case GAC_SHP_3D:
      {

	// Smooth out input volumes.
	this->m_smoothedNonSuppressedVolume = 
	  SmoothInputVolume( this->m_nonsuppressedVolume, 
			     MEDIAN );
	this->m_smoothedWaterSuppressedVolume = 
	  SmoothInputVolume( this->m_watersuppressedVolume, 
			     MEDIAN );
	this->m_smoothedFatSuppressedVolume = 
	  SmoothInputVolume( this->m_fatsuppressedVolume, 
			     MEDIAN );

	// Generate vector image from NS, WS and FS volumes
	// and apply median filtering to remove noise.
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
    case GMM_BAYES_3D: // Gaussian Mixture Modeling + Bayesian maximum ratio.
      {
	// // Smooth out input volumes.
	// this->m_smoothedNonSuppressedVolume = 
	//   SmoothInputVolume( this->m_nonsuppressedVolume, 
	// 		     MEDIAN );
	// this->m_smoothedWaterSuppressedVolume = 
	//   SmoothInputVolume( this->m_watersuppressedVolume, 
	// 		     MEDIAN );
	// this->m_smoothedFatSuppressedVolume = 
	//   SmoothInputVolume( this->m_fatsuppressedVolume, 
	// 		     MEDIAN );

	// Generate vector image from NS, WS and FS volumes.
	// this->FormNSWSFSVectorVolume();

	//       typedef itk::BayesianClassifierInitializationImageFilter< VectorFloatImageType > BayesianInitializerType;
	//       BayesianInitializerType::Pointer bayesianInitializer = BayesianInitializerType::New();
      }
      break;
    case KMEANS_MRF_3D: // Multi-parametric K-means followed by MRF.
      {

	// Smooth out input volumes.
	this->m_smoothedNonSuppressedVolume = 
	  SmoothInputVolume( this->m_nonsuppressedVolume, 
			     MEDIAN );
	this->m_smoothedWaterSuppressedVolume = 
	  SmoothInputVolume( this->m_watersuppressedVolume, 
			     MEDIAN );
	this->m_smoothedFatSuppressedVolume = 
	  SmoothInputVolume( this->m_fatsuppressedVolume, 
			     MEDIAN );

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


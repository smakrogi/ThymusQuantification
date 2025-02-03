/*===========================================================================

Program:   Quantification of thymus tissues from MRI.
Module:    $RCSfile: ThymusQuantification.cxx,v $
Language:  C++
Date:      $Date: 2011/10/21 11:04:00 $
Version:   $Revision: 1.0 $
Author:    S. K. Makrogiannis
3T MRI Facility National Institute on Aging/National Institutes of Health.

=============================================================================*/

#if defined(_MSC_VER)
#pragma warning ( disable : 4786 )
#endif

// Miscellaneous.
#include <cmath>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cstdio>
#include <ctime>

// Image and parameter I/O.
#include "ThymusQuantificationCLP.h"
#include "itkImageFileReader.h" 
#include "itkImageFileWriter.h" 
#include "itkImage.h"
//#include "itkOrientedImage.h"
// #include "itkGDCMImageIO.h"
// #include "itkGDCMSeriesFileNames.h"
// #include "itkImageSeriesReader.h"
#include "itkVectorImage.h"
//#include "itkImageToVectorImageFilter.h"
#include "itkScalarToArrayCastImageFilter.h"
//#include "itkCastImageFilter.h"
#include "itkImageRegionIterator.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "IO_Utils_ITK.h"

// Connected component processing.
#include "OptionList.h"
#include "CommonFunctions.h"
#include "ThymusSegmentationAlgorithms.h"


// Doxygen comments.
/*! \mainpage Thymus Quantification Index Page
 *
 */

/////////////////////////////////////////////////////////////
// Main function.
/////////////////////////////////////////////////////////////
int main(int argc, char *argv[])
{
  // Inputs will be: Non-suppressed, Fat-suppressed image, Water-Suppressed volumes
  // and thymus seed coordinates.


  // Algorithm outline.
  // 1. Non-linear smoothing by curvature anisotropic diffusion.
  // 2. Segmentation/Clustering of Thymus.
  // 3. Connected component labeling.
  // 4. Pick the label of the seed for fat.
  // 5. Compute volume of thymus.
  // 6. Compute fat and water ratios.


  // Switch off multi-threading.
  itk::MultiThreader::SetGlobalDefaultNumberOfThreads(1);
  itk::MultiThreader::SetGlobalMaximumNumberOfThreads(1);


  // Parse command line or Slicer arguments.
  PARSE_ARGS;

  // Main variables.
  int nClusters;
  float halfROIEdgeLength; // in mm.(1.0)
  float structureelementRadius = (float) STRUCTUREELEMENTRADIUS;
  bool debug = DEBUG;
  // std::string subjectID, 
  std::string  nonsuppressedfile,
    watersuppressedfile,
    fatsuppressedfile,
    outputPath,
    segmentationAlgorithm,
    gradientAlgorithm,
    levelsetInitializationAlgorithm,
    multiscaleGeneratorString;
  float sigmoidAlpha,
    sigmoidBeta,
    rgMultiplier,
    mrfSmoothFactor,
    levelsetPropagationScalingFactor,
    levelsetCurvatureScalingFactor,
    levelsetAdvectionScalingFactor,
    parzenBandwidth,
    parzenBandwidthBase;
  int segmentationMethod=GAC_3D,
    gradientMethod=PARZEN,
    levelsetInitializationMethod=SEED_BASED_ROI;
  int levelsetMaximumIterations,
    levelsetNumberOfPCAModes,
    rgIterations,
    mrfIterations,
    mrfNeighborhoodRadius,
    parzenRadius;
  float levelsetMaximumRMSError;
  double levelsetShapePriorScalingFactor,
    levelsetMAPContourWeight,
    levelsetMAPImageGradientWeight,
    levelsetMAPShapeWeight;
  double fastmarchingStoppingTime;
  double chanveseCurvatureWeight,
    chanveseAreaWeight,
    chanveseVolumeWeight,
    chanveseVolume,
    chanveseOverlapWeight;
  int levelsetNScales, multiscaleGenerator=PARZEN_EDGES;

  std::string nonsuppressedCompleteFilename;
  std::string fatsuppressedCompleteFilename;
  std::string watersuppressedCompleteFilename;
  FloatImageType::Pointer nonsuppressedVolume;
  FloatImageType::Pointer fatsuppressedVolume;
  FloatImageType::Pointer watersuppressedVolume;
  FloatImageType::PointType seedPoint; // in RAS
  std::vector<double> seedPointDouble;
  std::string outputFileName;
  clock_t startProcess, endProcess;


  // Set default values.
  outputPath = "./";
  segmentationAlgorithm = "GAC_3D";
  nClusters = CLUSTERS;
  gradientAlgorithm = "PARZEN";
  levelsetInitializationAlgorithm = "SEED_BASED_ROI";
  parzenRadius = PARZENRADIUS;
  parzenBandwidth = PARZENBANDWIDTH;
  parzenBandwidthBase = PARZENBANDWIDTHBASE;
  halfROIEdgeLength = HALFROIEDGELENGTH;
  rgMultiplier = RGMULTIPLIER; 
  rgIterations = RGITERATIONS;
  mrfIterations = MRFITERATIONS;
  mrfSmoothFactor = MRFSMOOTHINGFACTOR; 
  mrfNeighborhoodRadius = MRFNEIGHBORHOODRADIUS;
  sigmoidBeta = LEVELSETSIGMOIDBETA;
  sigmoidAlpha = (-1.0) * (LEVELSETSIGMOIDBETA / float(SIGMOIDBETAALPHARATIO));
  levelsetInitializationAlgorithm = "FAST_MARCHING";
  fastmarchingStoppingTime = FASTMARCHINGSTOPPINGTIME;
  levelsetPropagationScalingFactor = LEVELSETPROPAGATIONSCALINGFACTOR;
  levelsetCurvatureScalingFactor = LEVELSETCURVATURESCALINGFACTOR;
  levelsetAdvectionScalingFactor = LEVELSETADVECTIONSCALINGFACTOR;
  levelsetMaximumIterations = LEVELSETMAXITERATIONS;
  levelsetMaximumRMSError = LEVELSETMAXRMSERROR;
  levelsetShapePriorScalingFactor = LEVELSETSHAPEPRIORSCALINGFACTOR;
  levelsetNumberOfPCAModes =  LEVELSETPCAMODES;
  levelsetMAPContourWeight = LEVELSETMAPCONTOURWEIGHT;
  levelsetMAPImageGradientWeight = LEVELSETMAPIMAGEGRADIENTWEIGHT;
  levelsetMAPShapeWeight = LEVELSETMAPSHAPEWEIGHT;
  chanveseCurvatureWeight = CHANVESECURVATUREWEIGHT;    
  chanveseAreaWeight = CHANVESEAREAWEIGHT;
  chanveseVolumeWeight = CHANVESEVOLUMEWEIGHT;
  chanveseVolume = CHANVESEVOLUME;
  chanveseOverlapWeight = CHANVESEOVERLAPWEIGHT;
  levelsetNScales = NSCALES;
  multiscaleGeneratorString = "PARZEN_EDGES";


  // Pass values from xml file.
  nClusters = nClustersGUI;
  parzenRadius = parzenRadiusGUI;
  parzenBandwidth = parzenBandwidthGUI;
  parzenBandwidthBase = parzenBandwidthBaseGUI;
  halfROIEdgeLength = halfROIEdgeLengthGUI;
  rgMultiplier = rgMultiplierGUI;
  rgIterations = rgIterationsGUI;
  sigmoidBeta = sigmoidBetaGUI;
  sigmoidAlpha = (-1.0) * (sigmoidBetaGUI / sigmoidBetaAlphaRatioGUI);
  fastmarchingStoppingTime = fastmarchingStoppingTimeGUI;
  levelsetPropagationScalingFactor = levelsetPropagationScalingFactorGUI;
  levelsetCurvatureScalingFactor = levelsetCurvatureScalingFactorGUI;
  levelsetAdvectionScalingFactor = levelsetAdvectionScalingFactorGUI;
  levelsetMaximumIterations = levelsetMaximumIterationsGUI;
  levelsetMaximumRMSError = levelsetMaximumRMSErrorGUI;
  levelsetShapePriorScalingFactor = levelsetShapePriorScalingFactorGUI;
  levelsetNumberOfPCAModes =  levelsetNumberOfPCAModesGUI;
  levelsetMAPContourWeight = levelsetMAPContourWeightGUI;
  levelsetMAPImageGradientWeight = levelsetMAPImageGradientWeightGUI;
  levelsetMAPShapeWeight = levelsetMAPShapeWeightGUI;
  chanveseCurvatureWeight = chanveseCurvatureWeightGUI;
  chanveseAreaWeight = chanveseAreaWeightGUI;
  chanveseVolumeWeight = chanveseVolumeWeightGUI;
  chanveseVolume = chanveseVolumeGUI;
  chanveseOverlapWeight = chanveseOverlapWeightGUI;
  levelsetNScales = levelsetNScalesGUI;


  // Parse GUI inputs.
  nonsuppressedCompleteFilename = inputVolumeNS; 
  watersuppressedCompleteFilename = inputVolumeWS;
  fatsuppressedCompleteFilename = inputVolumeFS;
  seedPoint[0] = (float) (-1) * seed[0];
  seedPoint[1] = (float) (-1) * seed[1];
  seedPoint[2] = (float) seed[2];
  outputFileName = outputVolume;
  outputPath = outputDirectory;

    
  // Set the segmentation algorithm using string matching.
  int i = 0; 
  bool flag = false;
  segmentationAlgorithm = segmentationalgorithmGUI;
  while( i < (int) sizeof(SegmentationMethodStrings) && flag == false ) {
    if ( segmentationAlgorithm.compare(SegmentationMethodStrings[i]) == 0 ) {
      segmentationMethod = i;
      flag = true;
    }
    i++;
  }


  // Set the gradient algorithm using string matching.
  i = 0; flag = false;
  gradientAlgorithm = gradientalgorithmGUI;
  while( i < (int) sizeof(GradientMethodStrings) && flag == false ) {
    if ( gradientAlgorithm.compare(GradientMethodStrings[i]) == 0 ) {
      gradientMethod = i;
      flag = true;
    }
    i++;
  }


  // Set the level-set initialization algorithm using string matching.
  i = 0; flag = false;
  levelsetInitializationAlgorithm = levelsetinitializationalgorithmGUI;
  while( i < (int) sizeof(LevelsetInitializationMethodStrings) && flag == false ) {
    if ( levelsetInitializationAlgorithm.compare(LevelsetInitializationMethodStrings[i]) == 0 ) {
      levelsetInitializationMethod = i;
      flag = true;
    }
    i++;
  }


  // Set the scale-space generation algorithm using string matching.
   i = 0; flag = false;
   multiscaleGeneratorString = multiscaleGeneratorGUI;
    while( i < (int) sizeof(MultiscaleGeneratorMethodStrings) && flag == false ) {
    if ( multiscaleGeneratorString.compare(MultiscaleGeneratorMethodStrings[i]) == 0 ) {
      multiscaleGenerator = i;
      flag = true;
    }
    i++;
  }
 


  if ( flag == false ) {
    std::cerr << "Unknown segmentation or gradient method. "
  	      << std::endl;
    return EXIT_FAILURE;
  }

  std::cout << std::endl;
  std::cout << "Segmentation Algorithm: " 
	    << segmentationAlgorithm
	    << std::endl;
  std::cout << "Gradient Method: " 
	    << gradientAlgorithm
	    << std::endl;
  std::cout << "ROI initialization method: " 
	    << levelsetInitializationAlgorithm
	    << std::endl;
  std::cout << "Multi-scale generation method: " 
	    << multiscaleGeneratorString
	    << std::endl;
  std::cout << "ROI radius: "
	    << halfROIEdgeLength
	    << std::endl;


  // Read-in the volumes.

  typedef itk::ImageFileReader<FloatImageType> FloatReaderType;

  FloatReaderType::Pointer reader = FloatReaderType::New();

  reader->SetFileName( watersuppressedCompleteFilename.c_str() ); 
  try
    {
      reader->Update();
    }
  catch (itk::ExceptionObject &ex)
    {
      std::cout << ex << std::endl;
      return EXIT_FAILURE;
    }

  watersuppressedVolume = 
    ReorientImage<FloatImageType>( reader->GetOutput(), 
				   itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_RAI 
				   );
  reader = 0;

  std::cout << std::endl
	    << watersuppressedCompleteFilename 
	    << " Water-suppressed image read in." 
	    << std::endl;
  std::cout << "Matrix properties--"  
	    << "Size:" << watersuppressedVolume->GetBufferedRegion().GetSize() << ", "
	    << "Spacing:" <<  watersuppressedVolume->GetSpacing() 
	    << std::endl;


  reader = FloatReaderType::New();
  reader->SetFileName( fatsuppressedCompleteFilename.c_str() ); 
  try
    {
      reader->Update();
    }
  catch (itk::ExceptionObject &ex)
    {
      std::cout << ex << std::endl;
      return EXIT_FAILURE;
    }
  fatsuppressedVolume = 
    ReorientImage<FloatImageType>( reader->GetOutput(), 
				   itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_RAI 
				   );
  reader = 0;

  std::cout << std::endl
	    << fatsuppressedCompleteFilename 
	    << " Fat-suppressed image read in." 
	    << std::endl;
  std::cout << "Matrix properties--"  
	    << "Size:" << fatsuppressedVolume->GetBufferedRegion().GetSize() << ", "
	    << "Spacing:" <<  fatsuppressedVolume->GetSpacing() 
	    << std::endl;


  // If non-suppressed filename is empty,
  // add the two suppressed images.


  if( nonsuppressedCompleteFilename.compare( DefaultNSFilename ) == 0 ) {
    nonsuppressedVolume = AddDixonImages( watersuppressedVolume, 
					  fatsuppressedVolume );
    std::cout << std::endl
	      << " Added water and fat images, assuming Dixon scan types." 
	      << std::endl;
  }
  else {
    reader = FloatReaderType::New();
    reader->SetFileName( nonsuppressedCompleteFilename.c_str() ); 
    try
      {
	reader->Update();
      }
    catch (itk::ExceptionObject &ex)
      {
	std::cout << ex << std::endl;
	return EXIT_FAILURE;
      }
    nonsuppressedVolume = 
      ReorientImage<FloatImageType>( reader->GetOutput(), 
				     itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_RAI 
				     );
    reader = 0;

    std::cout << std::endl
	      << nonsuppressedCompleteFilename 
	      << " Non-suppressed image read in." 
	      << std::endl;
    std::cout << "Matrix properties--"  
	      << "Size:" << nonsuppressedVolume->GetBufferedRegion().GetSize() << ", "
	      << "Spacing:" <<  nonsuppressedVolume->GetSpacing() 
	      << std::endl;
  }

  // Start timing.
  startProcess = clock();


  // Generate a seed based neighborhood to be used for 
  // * determining region labels
  // * initializing region growing, level sets, etc :).

  std::vector<float> halfROISize(ImageDimension, halfROIEdgeLength);

  FloatImageType::RegionType seedbasedROI = GenerateROI ( nonsuppressedVolume,
							  seedPoint,
							  halfROISize );
  std::cout << "Seed ROI size  (px x px x px): " 
	    << seedbasedROI.GetSize()
	    << std::endl;


  // Instantiate Thymus segmentation object and run algorithms.

  ThymusSegmentationObject *thymusSegmentationModule = 
    new ThymusSegmentationObject;

  thymusSegmentationModule->SetInputs( nonsuppressedVolume,
				       watersuppressedVolume,
				       fatsuppressedVolume );
  

  thymusSegmentationModule->SetAlgorithmParameters( outputPath,
						    subjectID,
						    seedPoint,
						    seedbasedROI,
						    segmentationMethod,
						    nClusters,
						    gradientMethod,
						    parzenRadius,
						    parzenBandwidth,
						    parzenBandwidthBase,
						    structureelementRadius,
						    rgMultiplier,
						    rgIterations,
						    mrfIterations,
						    mrfSmoothFactor,
						    mrfNeighborhoodRadius,
						    sigmoidAlpha,
						    sigmoidBeta,
						    levelsetInitializationMethod,
						    fastmarchingStoppingTime,
						    levelsetPropagationScalingFactor,
						    levelsetCurvatureScalingFactor, 
						    levelsetAdvectionScalingFactor,
						    levelsetMaximumIterations,
						    levelsetMaximumRMSError,
						    levelsetShapePriorScalingFactor,
						    levelsetNumberOfPCAModes,
						    levelsetMAPContourWeight,
						    levelsetMAPImageGradientWeight,
						    levelsetMAPShapeWeight,
						    chanveseCurvatureWeight,
						    chanveseAreaWeight,
						    chanveseVolumeWeight,
						    chanveseVolume,
						    chanveseOverlapWeight,
						    multiscaleGenerator,
						    levelsetNScales); 
  thymusSegmentationModule->SetDebug ( debug );
  thymusSegmentationModule->RunThymusSegmentationAlgorithms();
  MaskImageType::Pointer outputMaskVolume = thymusSegmentationModule->GetOutput();
  infoStruct thymusInfoStruct = thymusSegmentationModule->DisplayVolumetrics();

  // Compute [water / (water + fat)] and [fat / (water + fat)] ratios over the segmented Thymus.
  // Allocate fat and water ratio volumes.
  FloatImageType::Pointer waterratioVolume = FloatImageType::New();
  waterratioVolume->CopyInformation( outputMaskVolume );
  waterratioVolume->SetRegions( outputMaskVolume->GetLargestPossibleRegion() );
  waterratioVolume->Allocate();
  waterratioVolume->FillBuffer( 0.0 );

  FloatImageType::Pointer fatratioVolume = FloatImageType::New();
  fatratioVolume->CopyInformation( outputMaskVolume );
  fatratioVolume->SetRegions( outputMaskVolume->GetLargestPossibleRegion() );
  fatratioVolume->Allocate();
  fatratioVolume->FillBuffer( 0.0 );

  // Define and instantiate image iterators.
  typedef itk::ImageRegionIteratorWithIndex<MaskImageType> MaskImageIteratorWithIndexType;
  MaskImageIteratorWithIndexType itMask(outputMaskVolume, 
					outputMaskVolume->GetBufferedRegion());
  vnl_vector<float> waterToSumRatioVector(thymusInfoStruct.SizeinPixels, 0.0), 
    fatToSumRatioVector(thymusInfoStruct.SizeinPixels, 0.0);
  // vnl_vector<float> waterToSumRatioVector, fatToSumRatioVector;
  float waterToSumRatio, fatToSumRatio, waterandfatSum;
  int count = 0;
  MaskImageType::IndexType idx, fsIdx;
  FloatImageType::PointType point;
  float minimalOffset = 10e-20;

  for(itMask.GoToBegin();!itMask.IsAtEnd();++itMask) {
    if (itMask.Get() != BACKGROUND ) {
      // Get index in the WS space.
      idx = itMask.GetIndex();
      // Convert index in WS space to physical coordinates.
      watersuppressedVolume->TransformIndexToPhysicalPoint(idx, point);
      // Convert physical coordinates to index in the FS space.
      fatsuppressedVolume->TransformPhysicalPointToIndex(point, fsIdx);
      // std::cout << itMask.Get() << idx << std::endl;

      waterandfatSum = ( watersuppressedVolume->GetPixel( idx ) + 
			 fatsuppressedVolume->GetPixel( fsIdx ) + 
			 minimalOffset );

      waterToSumRatio = fatsuppressedVolume->GetPixel(fsIdx) / waterandfatSum;

      fatToSumRatio = watersuppressedVolume->GetPixel( idx ) / waterandfatSum;

      //    waterToSumRatio = fatsuppressedVolume->GetPixel( idx ) / 
      //       ( fatsuppressedVolume->GetPixel( idx ) + watersuppressedVolume->GetPixel( idx ) + minimalOffset);
      //     fatToSumRatio = watersuppressedVolume->GetPixel( idx ) / 
      //       ( fatsuppressedVolume->GetPixel( idx ) + watersuppressedVolume->GetPixel( idx ) + minimalOffset);

      waterratioVolume->SetPixel(idx, waterToSumRatio);
      fatratioVolume->SetPixel(idx, fatToSumRatio);
      waterToSumRatioVector[count] = waterToSumRatio;
      fatToSumRatioVector[count] = fatToSumRatio;
      count++;
      }

    }

  // Compute rms error.
  vnl_vector<float> waterToSumRatioVectorCentered(thymusInfoStruct.SizeinPixels, 0.0), 
    fatToSumRatioVectorCentered(thymusInfoStruct.SizeinPixels, 0.0);
  waterToSumRatioVectorCentered = waterToSumRatioVector - waterToSumRatioVector.mean();
  fatToSumRatioVectorCentered = fatToSumRatioVector - fatToSumRatioVector.mean();

  // Display volumetric water and fat ratios on console and Slicer.
  std::stringstream outputString;
  outputString << std::endl
	       << "=======Water and Fat Ratios over Thymus mask======" << std::endl
	       << "[fat / (water + fat)] ratio stats" << std::endl
	       << "Minimum: " << fatToSumRatioVector.min_value() << std::endl
	       << "Maximum: " << fatToSumRatioVector.max_value() << std::endl
	       << "Average: " << fatToSumRatioVector.mean() << std::endl
	       << "RMS error: " << fatToSumRatioVectorCentered.rms() << std::endl
	       << "[water / (water + fat)] ratio stats" << std::endl
	       << "Minimum: " << waterToSumRatioVector.min_value() << std::endl
	       << "Maximum: " << waterToSumRatioVector.max_value() << std::endl
	       << "Average: " << waterToSumRatioVector.mean() << std::endl
	       << "RMS error: " << waterToSumRatioVectorCentered.rms() << std::endl
	       << "=======================================" 
	       << std::endl;
  std::cout << outputString.str();

  // Stop timing and display elapsed time.
  endProcess = clock();
  double elapsedTime = ( double(endProcess)-
			 double(startProcess) ) / CLOCKS_PER_SEC;
 
  std::cout << "Elapsed time for quantification: " 
            << elapsedTime
            << " sec"
            << std::endl;


  // Write to csv file.
  std::string csvFilename = outputPath + "/" + subjectID + 
    ".stats" + ".csv";
  std::ofstream out( csvFilename.c_str() ); 

  if(!out) { 
    std::cerr << "Cannot open file.\n"; 
  }
  out << "SubjectId,ThymusSize(px),ThymusVol(mm^3),WaterRatioMean,WaterRatioRMS,FatRatioMean,FatRatioRMS" << std::endl
      << subjectID << "," << thymusInfoStruct.SizeinPixels << "," 
      << thymusInfoStruct.SizeinPhysicalUnits << ","
      << waterToSumRatioVector.mean() << "," << waterToSumRatioVectorCentered.rms() << ","
      << fatToSumRatioVector.mean() << "," << fatToSumRatioVectorCentered.rms();
  out.close();


  // Write out the return parameters in "name = value" form
  std::ofstream rts;
  rts.open(returnParameterFile.c_str());
  rts 
    << "volumepxReturnIntGUI = " << thymusInfoStruct.SizeinPixels << std::endl
    << "volumemm3ReturnFloatGUI = " << thymusInfoStruct.SizeinPhysicalUnits << std::endl
    << "fatratiomeanReturnFloatGUI = " << fatToSumRatioVector.mean() << std::endl
    << "fatratiormsReturnFloatGUI = " << fatToSumRatioVectorCentered.rms() << std::endl
    << "waterratiomeanReturnFloatGUI = " << waterToSumRatioVector.mean() << std::endl
    << "waterratiormsReturnFloatGUI = " << waterToSumRatioVectorCentered.rms() << std::endl;
  rts.close();


  // Write water and fat ratio volumes to files.

  itk::ImageFileWriter<FloatImageType>::Pointer floatWriter = 
    itk::ImageFileWriter<FloatImageType>::New();
  
  std::string waterratioFilename = outputPath + "/" + subjectID + "_" + "segmentedWaterRatio.nii";
  floatWriter = itk::ImageFileWriter<FloatImageType>::New();
  floatWriter->SetInput( waterratioVolume ); 
  floatWriter->SetFileName( waterratioFilename );
  floatWriter->Update();
  floatWriter = 0;
  
  std::string fatratioFilename = outputPath + "/" + subjectID + "_" + "segmentedFatRatio.nii";
  floatWriter =  itk::ImageFileWriter<FloatImageType>::New();
  floatWriter->SetInput( fatratioVolume ); 
  floatWriter->SetFileName( fatratioFilename );
  floatWriter->Update();
  floatWriter = 0;


  // Save output mask to a volume.
  if( outputFileName.compare( "MaskImage" ) == 0 )
    outputFileName = outputPath + "/" + subjectID + "_" + outputFileName + ".nii";
  std::cout << std::endl;
  std::cout << "Saving mask image to: " 
	    << outputFileName.c_str()
	    << std::endl;
  itk::ImageFileWriter<MaskImageType>::Pointer labelWriter = 
    itk::ImageFileWriter<MaskImageType>::New();
  labelWriter->SetInput( outputMaskVolume ); 
  labelWriter->SetFileName( outputFileName.c_str() );
  labelWriter->Update();
  labelWriter = 0;

  delete( thymusSegmentationModule );

  return EXIT_SUCCESS;
}

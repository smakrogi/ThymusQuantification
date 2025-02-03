/*===========================================================================

Program:   Quantification of thymus tissues from MRI.
Module:    $RCSfile: ThymusQuantification.cxx,v $
Language:  C++
Date:      $Date: 2010/04/06 10:42:32 $
Version:   $Revision: 0.5 $
Author:    S. K. Makrogiannis
3T MRI Facility National Institute on Aging/National Institutes of Health.

=============================================================================*/

#if defined(_MSC_VER)
#pragma warning ( disable : 4786 )
#endif

// Image and parameter I/O.
#include "itkImageFileReader.h" 
#include "itkImageFileWriter.h" 
#include "OptionList.h"

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

// Connected component processing.
// #include "itkConnectedComponentImageFilter.h"
// #include "itkRelabelComponentImageFilter.h"
// #include "itkBinaryThresholdImageFilter.h"

// Miscellaneous.
#include <cmath>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cstdio>
#include <ctime>

#include "OptionList.h"
#include "CommonFunctions.h"
#include "ThymusSegmentationAlgorithms.h"


// Doxygen comments.
/*! \mainpage Thymus Quantification Index Page
 *
 */

void print_usage()
{
  std::cout << std::endl;
  std::cout << "3T NIA ThymusQuantification 0.6 (Apr. 27 2010 by SKM)" << std::endl;
  std::cout << std::endl;
  std::cout << "usage: ThymusQuantification --subjectid <ID>... " << std::endl;
  std::cout << "       --nonsuppressed <imagefile>... "  << std::endl ;
  std::cout << "       --watersuppressed <imagefile>... "  << std::endl ;
  std::cout << "       --fatsuppressed <imagefile>... " << std::endl ;
  std::cout << "       --seedpoint <x> <y> <z>... " << std::endl;
  std::cout << "       --outputpath <path>... " << std::endl;
  std::cout << "       --segmentalgo {WS_KMEANS,KMEANS_3D,REGIONGROWING_3D,WS_LEVEL_SETS," << std::endl;
  std::cout << "                      LEVEL_SETS_3D,WS_GEOD_ACT_CONT,GEOD_ACT_CONT_3D, " << std::endl;
  std::cout << "                      EM_GMM_BAYES_3D,KMEANS_MRF_3D}... (optional)" << std::endl;
  std::cout << "       --nclusters <int>..(optional). " << std::endl;
  std::cout << "       --gradalgo {CONVENTIONAL,PARZEN}...(optional). " << std::endl;
  std::cout << "       --parzrad  <int>...(optional). " << std::endl;
  std::cout << "       --parzband <float>...(optional). " << std::endl;
  std::cout << "       --roiradius <float>...(optional). " << std::endl;
  std::cout << "       --rgmult <int>...(optional). " << std::endl; 
  std::cout << "       --rgit <float>...(optional). " << std::endl; 
  std::cout << "       --mrfit <int>...(optional). " << std::endl; 
  std::cout << "       --mrfsmooth <float>...(optional). " << std::endl; 
  std::cout << "       --mrfnbrradius <int>...(optional). " << std::endl; 
  std::cout << "       --lssigmoidalpha <float>...(optional). " << std::endl;
  std::cout << "       --lssigmoidbeta <float>...(optional). " << std::endl;
  std::cout << "       --lspropscale <float>...(optional). " << std::endl;
  std::cout << "       --lscurvscale <float>...(optional). " << std::endl;
  std::cout << "       --lsadvscale <float>...(optional). " << std::endl;
  std::cout << "       --lsmaxit <int>...(optional). " << std::endl;
  std::cout << "       --lsmaxrmserror <float>...(optional). " << std::endl; }


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
  // 6. Compute water to fat ratios.

  // Main variables.
  int segmentationMethod = KMEANS_3D;
  // {WS_KMEANS=0,KMEANS_3D,REGIONGROWING_3D,WS_LEVEL_SETS,LEVEL_SETS_3D,EM_GMM_BAYES_3D,KMEANS_MRF_3D} 
  int gradientMethod;
  int nClusters;
  float halfROIEdgeLength; // in mm.(1.0)
  float structureelementRadius = (float) STRUCTUREELEMENTRADIUS;
  bool debug = DEBUG;
  std::string subjectID, 
    nonsuppressedfile,
    watersuppressedfile,
    fatsuppressedfile,
    outputPath,
    segmentationAlgorithm,
    gradientAlgorithm;
  float sigmoidAlpha,
    sigmoidBeta,
    rgMultiplier,
    mrfSmoothFactor,
    levelsetPropagationScalingFactor,
    levelsetCurvatureScalingFactor,
    levelsetAdvectionScalingFactor,
    parzenBandwidth;
  int levelsetMaximumIterations,
    rgIterations,
    mrfIterations,
    mrfNeighborhoodRadius,
    parzenRadius;
  float levelsetMaximumRMSError;

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

   if (argc <= 1)
    {
      print_usage() ;
      exit(EXIT_FAILURE) ;
    }

   OptionList options(argc, argv) ;

   try
    {
      // Reguired arguments.
      options.GetStringOption("subjectid", &subjectID, true) ;
      options.GetStringOption("nonsuppressed", &nonsuppressedfile, true) ;
      options.GetStringOption("watersuppressed", &watersuppressedfile, true) ;
      options.GetStringOption("fatsuppressed", &fatsuppressedfile, true) ;
      options.GetStringOption("outputpath", &outputPath, true) ;
      options.GetMultiDoubleOption("seedpoint", &seedPointDouble, true);

      // Optional arguments.
      // Segmentation algorithm.
      if ( options.GetStringOption("segmentalgo", &segmentationAlgorithm, false) == -1 )
	segmentationAlgorithm = "KMEANS_3D";
      // Clustering.
      nClusters = options.GetIntOption("nclusters", false);
      if ( nClusters == -1 ) nClusters = CLUSTERS;
      // Image gradient algorithm.
      if ( options.GetStringOption("gradalgo", &gradientAlgorithm, false) == -1 )
	gradientAlgorithm = "CONVENTIONAL";
      parzenRadius = options.GetIntOption("parzrad", false);
      if ( parzenRadius == -1 ) parzenRadius = PARZENRADIUS;
      parzenBandwidth = (float) options.GetDoubleOption("parzband", false);
      if ( parzenBandwidth == -1 ) parzenBandwidth = PARZENBANDWIDTH;
      // ROI size.
      halfROIEdgeLength = (float) options.GetDoubleOption("roiradius", false);
      if ( halfROIEdgeLength == -1 ) halfROIEdgeLength = HALFROIEDGELENGTH;
      // Region growing.
      rgMultiplier = (float) options.GetDoubleOption("rgmult", false);
      if ( rgMultiplier == -1 )rgMultiplier = RGMULTIPLIER; 
      rgIterations = options.GetIntOption("rgit", false);
      if ( rgIterations == -1 ) rgIterations = RGITERATIONS;
      // Markov Random Fields.
      mrfIterations = options.GetIntOption("mrfit", false);
      if ( mrfIterations == -1 ) mrfIterations = MRFITERATIONS;
      mrfSmoothFactor = (float) options.GetDoubleOption("mrfsmooth", false);
      if ( mrfSmoothFactor == -1 )mrfSmoothFactor = MRFSMOOTHINGFACTOR; 
      mrfNeighborhoodRadius = options.GetIntOption("mrfnbrradius", false);
      if ( mrfNeighborhoodRadius == -1 ) 
	mrfNeighborhoodRadius = MRFNEIGHBORHOODRADIUS;
      // Level sets.
      sigmoidAlpha = (float) options.GetDoubleOption("lssigmoidalpha", false);
      if ( sigmoidAlpha == -1 ) sigmoidAlpha = LEVELSETSIGMOIDALPHA;
      sigmoidBeta = (float) options.GetDoubleOption("lssigmoidbeta", false);
      if ( sigmoidBeta == -1 ) sigmoidBeta = LEVELSETSIGMOIDBETA;
      levelsetPropagationScalingFactor = (float) options.GetDoubleOption("lspropscale", false);
      if ( levelsetPropagationScalingFactor == -1 ) 
	levelsetPropagationScalingFactor = LEVELSETPROPAGATIONSCALINGFACTOR;
      levelsetCurvatureScalingFactor = (float) options.GetDoubleOption("lscurvscale", false);
      if ( levelsetCurvatureScalingFactor == -1 ) 
	levelsetCurvatureScalingFactor = LEVELSETCURVATURESCALINGFACTOR;
      levelsetAdvectionScalingFactor = (float) options.GetDoubleOption("lsadvscale", false);
      if ( levelsetAdvectionScalingFactor == -1 ) 
	levelsetAdvectionScalingFactor = LEVELSETADVECTIONSCALINGFACTOR;
      levelsetMaximumIterations = options.GetIntOption("lsmaxit", false);
      if ( levelsetMaximumIterations == -1 ) 
	levelsetMaximumIterations = LEVELSETMAXITERATIONS;
      levelsetMaximumRMSError = (float) options.GetDoubleOption("lsmaxrmserror", false);
      if ( levelsetMaximumRMSError == -1 ) 
	levelsetMaximumRMSError = LEVELSETMAXRMSERROR;
}
  catch(OptionList::RequiredOptionMissing e)
    {
      std::cout << "Error: The '" << e.OptionTag 
                << "' option is required but missing." 
                << std::endl ;
      print_usage() ;
      return EXIT_FAILURE ;
    }


  // Generate filenames with path and variables.
//   nonsuppressedCompleteFilename = subjectID + "/" + nonsuppressedfile; 
//   watersuppressedCompleteFilename = subjectID + "/" + watersuppressedfile;
//   fatsuppressedCompleteFilename = subjectID + "/" + fatsuppressedfile;
   nonsuppressedCompleteFilename = nonsuppressedfile; 
   watersuppressedCompleteFilename = watersuppressedfile;
   fatsuppressedCompleteFilename = fatsuppressedfile;
  seedPoint[0] = (float) seedPointDouble[0];  
  seedPoint[1] = (float) seedPointDouble[1]; 
  seedPoint[2] = (float) seedPointDouble[2]; 
  outputFileName = outputPath + "/" + 
    subjectID + "_" + gradientAlgorithm + "_" + 
    segmentationAlgorithm + "_" + "ThymusMask.nii";

    
  // Set the segmentation algorithm.
  int i = 0; 
  bool flag = false;
  while( i < (int) sizeof(SegmentationMethodString) && flag == false ) {
    if ( segmentationAlgorithm.compare(SegmentationMethodString[i]) == 0 ) {
      segmentationMethod = i;
      flag = true;
    }
    i++;
  }

  // Set the gradient algorithm.
  i = 0; flag = false;
  while( i < (int) sizeof(GradientMethodString) && flag == false ) {
    if ( gradientAlgorithm.compare(GradientMethodString[i]) == 0 ) {
      gradientMethod = i;
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
  std::cout << "ROI radius: "
	    << halfROIEdgeLength
	    << std::endl;

  // Read-in the volumes.

  typedef itk::ImageFileReader<FloatImageType> FloatReaderType;

  FloatReaderType::Pointer reader = FloatReaderType::New();

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
  nonsuppressedVolume = reader->GetOutput();
  reader = 0;

  std::cout << std::endl
	    << nonsuppressedCompleteFilename 
	    << " Non-suppressed image read in." 
	    << std::endl;
  std::cout << "Matrix properties--"  
	    << "Size:" << nonsuppressedVolume->GetBufferedRegion().GetSize() << ", "
	    << "Spacing:" <<  nonsuppressedVolume->GetSpacing() 
	    << std::endl;

  reader = FloatReaderType::New();
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
  watersuppressedVolume = reader->GetOutput();
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
  fatsuppressedVolume = reader->GetOutput();
  reader = 0;

  std::cout << std::endl
	    << fatsuppressedCompleteFilename 
	    << " Fat-suppressed image read in." 
	    << std::endl;
  std::cout << "Matrix properties--"  
	    << "Size:" << fatsuppressedVolume->GetBufferedRegion().GetSize() << ", "
	    << "Spacing:" <<  fatsuppressedVolume->GetSpacing() 
	    << std::endl;


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
						    structureelementRadius,
						    rgMultiplier,
						    rgIterations,
						    mrfIterations,
						    mrfSmoothFactor,
						    mrfNeighborhoodRadius,
						    sigmoidAlpha,
						    sigmoidBeta,
						    levelsetPropagationScalingFactor,
						    levelsetCurvatureScalingFactor, 
						    levelsetAdvectionScalingFactor,
						    levelsetMaximumIterations,
						    levelsetMaximumRMSError ); 
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
  float waterToSumRatio, fatToSumRatio;
  int count = 0;
  MaskImageType::IndexType idx;
  float minimalOffset = 10e-20;

  for(itMask.GoToBegin();!itMask.IsAtEnd();++itMask) {
    idx = itMask.GetIndex();
    // std::cout << itMask.Get() << idx << std::endl;

    waterToSumRatio = fatsuppressedVolume->GetPixel( idx ) / 
      ( fatsuppressedVolume->GetPixel( idx ) + watersuppressedVolume->GetPixel( idx ) + minimalOffset);
    fatToSumRatio = watersuppressedVolume->GetPixel( idx ) / 
      ( fatsuppressedVolume->GetPixel( idx ) + watersuppressedVolume->GetPixel( idx ) + minimalOffset);

    waterratioVolume->SetPixel(idx, waterToSumRatio);
    fatratioVolume->SetPixel(idx, fatToSumRatio);
    
    if (itMask.Get() != BACKGROUND ) {
      waterToSumRatioVector[count] = waterToSumRatio;
      fatToSumRatioVector[count] = fatToSumRatio;
      count++;
    }

  }


  // Display volumetric results.

  std::cout << "=======Water and Fat Ratios over Thymus mask======" 
	    << std::endl;
  std::cout << "[water / (water + fat)] ratio stats" << std::endl
	    << "Minimum: " << waterToSumRatioVector.min_value() << std::endl
	    << "Maximum: " << waterToSumRatioVector.max_value() << std::endl
	    << "Average: " << waterToSumRatioVector.mean() << std::endl
	    << "RMS: " << waterToSumRatioVector.rms() << std::endl;
  std::cout << "[fat / (water + fat)] ratio stats" << std::endl
	    << "Minimum: " << fatToSumRatioVector.min_value() << std::endl
	    << "Maximum: " << fatToSumRatioVector.max_value() << std::endl
	    << "Average: " << fatToSumRatioVector.mean() << std::endl
	    << "RMS: " << fatToSumRatioVector.rms() << std::endl;
  std::cout << "=======================================" 
	    << std::endl;



  // Stop timing and display elapsed time.
  endProcess = clock();
  double elapsedTime = ( double(endProcess)-
			 double(startProcess) ) / CLOCKS_PER_SEC;
 
  std::cout << "Elapsed time for quantification: " 
            << elapsedTime
            << " sec"
            << std::endl;


  // Write to csv file.
  std::string csvFilename = subjectID + "/" + subjectID + "_" +
    gradientAlgorithm + "_" + segmentationAlgorithm + 
    ".stats" + ".csv";
  std::ofstream out( csvFilename.c_str() ); 

  if(!out) { 
    std::cerr << "Cannot open file.\n"; 
  }
  out << "SubjectId,ThymusSize(px),ThymusVol(mm^3),WaterRatioMean,WaterRatioRMS,FatRatioMean,FatRatioRMS" << std::endl
      << subjectID << "," << thymusInfoStruct.SizeinPixels << "," << thymusInfoStruct.SizeinPhysicalUnits << ","
      << waterToSumRatioVector.mean() << "," << waterToSumRatioVector.rms() << ","
      << fatToSumRatioVector.mean() << "," << fatToSumRatioVector.rms();
  
  out.close();


  // Write water and fat ratio volumes to files.

  itk::ImageFileWriter<FloatImageType>::Pointer floatWriter = itk::ImageFileWriter<FloatImageType>::New();
  
  std::string waterratioFilename = outputPath + "/" + subjectID + "_" + "WaterRatio.nii";
  floatWriter = itk::ImageFileWriter<FloatImageType>::New();
  floatWriter->SetInput( waterratioVolume ); 
  floatWriter->SetFileName( waterratioFilename );
  floatWriter->Update();
  floatWriter = 0;
  
  std::string fatratioFilename = outputPath + "/" + subjectID + "_" + "FatRatio.nii";
  floatWriter =  itk::ImageFileWriter<FloatImageType>::New();
  floatWriter->SetInput( fatratioVolume ); 
  floatWriter->SetFileName( fatratioFilename );
  floatWriter->Update();
  floatWriter = 0;


  // Save output mask to a volume.
  itk::ImageFileWriter<MaskImageType>::Pointer labelWriter = 
    itk::ImageFileWriter<MaskImageType>::New();
  labelWriter->SetInput( outputMaskVolume ); 
  labelWriter->SetFileName( outputFileName.c_str() );
  labelWriter->Update();
  labelWriter = 0;

  delete( thymusSegmentationModule );

  return EXIT_SUCCESS;
}


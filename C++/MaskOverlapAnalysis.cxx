/*===========================================================================

  Program:   Volumetric Image Overlap Analysis.
  Module:    $RCSfile: MaskOverlapAnalysis.cxx,v $
  Language:  C++
  Date:      $Date: 2011/03/02 13:40:00 $
  Version:   $Revision: 0.1 $

  3T MRI Facility National Institute on Aging/National Institutes of Health.

  =============================================================================*/

#if defined(_MSC_VER)
#pragma warning ( disable : 4786 )
#endif

#include <iostream>
#include <fstream>

#include "MaskOverlapAnalysisCLP.h"
#include "itkImageFileReader.h"
#include "OrientedImageTypes.h"
#include "IO_Utils_ITK.h"
#include "ComputeOverlap.h"
#include "itkBinaryThresholdImageFilter.h"
#include "itkImageFileWriter.h" 

using namespace std;

int main(int argc, char* argv[])
{
  itk::MultiThreader::SetGlobalDefaultNumberOfThreads(1);
  itk::MultiThreader::SetGlobalMaximumNumberOfThreads(1);

//   std::string segmentationFile;
//   std::string referenceStandardFile;
//   string csvOutputFile;
//  itk::OStringStream msg;

  ofstream fd_csv;


	
  typedef itk::BinaryThresholdImageFilter<OrientedInputImageType, OrientedMaskImageType>
    ThresholdType;


  // Parse command line or Slicer arguments.
  PARSE_ARGS;


  // Read in segmented image.

  OrientedInputImageType::Pointer segmentedInput;

  std::string subjectName;

  if ( ReadArchetypeImage( segmentationFile, segmentedInput, subjectName ) ) {
    std::cout << segmentationFile << " Segmented image read in" << std::endl;
    std::cout << "Subject name = " << subjectName << std::endl;
  }
  else {
    std::cerr << "Error reading segmented image: " << segmentationFile << std::endl;
    return EXIT_FAILURE;
  }

  OrientedInputImageType::Pointer segmentedImage = 
    ReorientImage<OrientedInputImageType>( segmentedInput,
					   itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_RAI);
  segmentedInput = 0;


  // Construct a segmentation mask (find all pixels > 0) from the segmentation mask.
  //

  OrientedMaskImageType::Pointer segmentedmaskImage;

  ThresholdType::Pointer threshold = ThresholdType::New();

  threshold->SetInput( segmentedImage );
  threshold->SetInsideValue( FOREGROUND );
  threshold->SetOutsideValue( BACKGROUND );
  threshold->SetLowerThreshold( 1 );
  threshold->Update();

  segmentedmaskImage = threshold->GetOutput();
  // threshold = 0;

  // Read in mask/reference image using archetype reader.

  OrientedInputImageType::Pointer referenceInput;

  if ( ReadArchetypeImage( referenceStandardFile, referenceInput, subjectName ) ) {
    std::cout << referenceStandardFile << " Ground Truth image read in" << std::endl;
    std::cout << "Subject name = " << subjectName << std::endl;
  }
  else {
    std::cerr << "Error reading reference image: " << referenceStandardFile << std::endl;
    return EXIT_FAILURE;
  }

  OrientedInputImageType::Pointer referenceImage = 
    ReorientImage<OrientedInputImageType>( referenceInput,
					   itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_RAI);
  referenceInput = 0;

  // Construct an ROI mask (find all pixels > 0) from the ground truth.
  //

  OrientedMaskImageType::Pointer referencemaskImage;

  ThresholdType::Pointer threshold2 = ThresholdType::New();
  threshold2->SetInput( referenceImage );
  threshold2->SetInsideValue( FOREGROUND );
  threshold2->SetOutsideValue( BACKGROUND );
  threshold2->SetLowerThreshold( 1 );
  // threshold2->SetUpperThreshold( 1 );
  threshold2->Update();

  referencemaskImage = threshold2->GetOutput();
  // threshold2 = 0;


//   // Write image for debug purposes.
//   itk::ImageFileWriter<OrientedMaskImageType>::Pointer writer = 
//     itk::ImageFileWriter<OrientedMaskImageType>::New();
//   std::string writerfilename = "./debug.nii";
//   writer = itk::ImageFileWriter<OrientedMaskImageType>::New();
//   writer->SetInput( segmentedmaskImage ); 
//   writer->SetFileName( writerfilename );
//   writer->Update();
//   writer = 0; 

 

  // Open output csv file for writing.
  fd_csv.open(csvOutputFile.c_str());
  if (fd_csv.fail()) {
    cerr << "unable to open file for writing" << endl;
    exit(1);
  }
  fd_csv << "Subject_ID,Volume_groundtruth(mm^3),Volume_detection(mm^3),Volume_overlap(mm^3),Volume_union(mm^3),";
  fd_csv << "Overlap_score,DICE_score,Sensitivity(%),Specificity(%),Classification_Accuracy(%),Volume_Difference(%)" << endl;

  
  // Write out the return parameters in "name = value" form
  std::ofstream slicer_stream;
  slicer_stream.open(returnParameterFile.c_str());

  // Compute detection rates.
  ComputeOverlapScore( referencemaskImage, 
		       segmentedmaskImage, 
		       subjectID,
		       fd_csv, 
		       slicer_stream );

  fd_csv.close();
  slicer_stream.close();

 return 0;
}

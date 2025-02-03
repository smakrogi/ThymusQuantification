#if defined(_MSC_VER)
#pragma warning ( disable : 4786 )
#endif

#ifdef __BORLANDC__
#define ITK_LEAN_AND_MEAN
#endif

#include "itkImageRegionIteratorWithIndex.h"
#include "ComputeOverlap.h"
#include <cmath>

void
ComputeOverlapScore( OrientedMaskImageType::Pointer groundTruth, 
		     OrientedMaskImageType::Pointer segmentation, 
		     std::string subjectID,
		     ofstream &f_csv,
		     ofstream &slicer_stream)
{
  // compute the overlap score;
  float nPixelsInGroundTruth = 0;
  float nPixelsInDetection = 0;
  float nPixelsInOverlap = 0;
  float nPixelsInUnion = 0;
  float FP = 0;
  float FN = 0;
  float TP = 0;
  float TN = 0; // TP + TN = groundTruth

  itk::ImageRegionConstIteratorWithIndex< OrientedMaskImageType > 
    it( segmentation, 
	segmentation->GetLargestPossibleRegion() );
  OrientedMaskImageType::PointType point;
  OrientedMaskImageType::IndexType testIdx, referenceIdx;

  for (it.GoToBegin(); !it.IsAtEnd(); ++it) { 
    testIdx = it.GetIndex();
    segmentation->TransformIndexToPhysicalPoint( testIdx,
						 point );
    groundTruth->TransformPhysicalPointToIndex( point,
						referenceIdx );
    if( groundTruth->GetBufferedRegion().IsInside( referenceIdx ) ) {
	if ( groundTruth->GetPixel( referenceIdx ) == FOREGROUND ) {
	  nPixelsInGroundTruth++;
	  nPixelsInUnion++;
	  if ( it.Get() == FOREGROUND ) {
	    // got it
	    nPixelsInDetection++;
	    nPixelsInOverlap++;
	    TP++;
	  }
	  else { // missed
	    FN++;
	  }
	}
	else {
	  // could be false positive
	  if ( it.Get() == FOREGROUND ) {
	    nPixelsInDetection++;
	    nPixelsInUnion++;
	    FP++;
	  }
	  else {
	    // correct!
	    TN++;
	  }
	}
      }
  }

  float voxelVolume = groundTruth->GetSpacing()[0] *
    groundTruth->GetSpacing()[1] *
    groundTruth->GetSpacing()[2];
  float overlapScore = 
    static_cast<float>(nPixelsInOverlap)/static_cast<float>(nPixelsInUnion);
  float diceScore = 
    static_cast<float>(nPixelsInOverlap)*2/static_cast<float>(nPixelsInGroundTruth+nPixelsInDetection);
  float sensitivity = static_cast<float>(100*TP)/static_cast<float>(TP+FN);
  float specificity = static_cast<float>(100*TN)/static_cast<float>(TN+FP);
  float classificationAccuracy = static_cast<float>( 100 * ( TP + TN ) ) / 
    static_cast<float>( TP + FN + TN + FP  );
  float volumeDifference = static_cast<float>(nPixelsInDetection - nPixelsInGroundTruth) 
    / static_cast<float>(nPixelsInGroundTruth);
  volumeDifference =  100.0 * std::abs(volumeDifference);

  // To a file output.
  f_csv << subjectID << ","
	<< nPixelsInGroundTruth * voxelVolume << ","
        << nPixelsInDetection * voxelVolume << ","
        << nPixelsInOverlap * voxelVolume << ","
        << nPixelsInUnion * voxelVolume << "," 
        << overlapScore << "," 
        << diceScore << "," 
        << sensitivity << "," 
        << specificity << ","
	<< classificationAccuracy << ","
	<< volumeDifference 
        << std::endl;


  // Write out the return parameters in "name = value" form for 3D Slicer.
  slicer_stream 
    << "volumeInGroundTruthGUI = " << nPixelsInGroundTruth * voxelVolume << std::endl
    << "volumeInDetectionGUI = " << nPixelsInDetection * voxelVolume << std::endl
    << "volumeInOverlapGUI = " << nPixelsInOverlap * voxelVolume << std::endl
    << "volumeInUnionGUI = " << nPixelsInUnion * voxelVolume << std::endl
    << "overlapScoreGUI = " << overlapScore << std::endl
    << "diceScoreGUI = " << diceScore << std::endl
    << "classificationAccuracyGUI = " << classificationAccuracy << std::endl
    << "volumeDifferenceGUI = " << volumeDifference << std::endl;
 

  // To standard output.
  std::cout << "Subject ID: " << subjectID << std::endl;
std::cout << "NPixels in groundtruth: " << nPixelsInGroundTruth << std::endl;
  std::cout << "NPixels in detection: " << nPixelsInDetection << std::endl;
  std::cout << "NPixels in overlap: " << nPixelsInOverlap << std::endl;
  std::cout << "NPixels in union: " << nPixelsInUnion << std::endl;
  std::cout << "Overlap score: " 
	    << overlapScore
	    << std::endl;
  std::cout << "DICE score: " 
	    << diceScore
	    << std::endl;
  std::cout << "Sensitivity: " 
	    << sensitivity << "%" 
	    << std::endl;
  std::cout << "Specificity: " 
	    << specificity << "%" 
	    << std::endl;
  std::cout << "Classification Accuracy: " 
	    << classificationAccuracy << "%" 
	    << std::endl;
  std::cout << "Volume Difference: " 
	    << volumeDifference << "%" 
	    << std::endl;

  return;  
}

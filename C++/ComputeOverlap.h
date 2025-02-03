#ifndef __ComputeOverlap_h
#define __ComputeOverlap_h

#include <iostream>
#include <fstream>
#include "OrientedImageTypes.h"

#define FOREGROUND 255
#define BACKGROUND 0

using namespace std;

void
ComputeOverlapScore( OrientedMaskImageType::Pointer groundTruth,
		     OrientedMaskImageType::Pointer segmentation, 
		     std::string subjectID,
		     ofstream &f_csv,
		     ofstream &slicer_stream );

#endif

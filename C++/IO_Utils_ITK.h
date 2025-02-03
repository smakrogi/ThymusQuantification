#ifndef __IO_Utils_ITK_h
#define __IO_Utils_ITK_h

#include "OrientedImageTypes.h"
#include "itkSpatialOrientation.h"

// template<class ImageType> 
// itk::SmartPointer<ImageType> ReorientImage( itk::SmartPointer<ImageType> img,
//                                   itk::SpatialOrientation::ValidCoordinateOrientationFlags orientation );

// routine to enforce that the CT image and the mask have the same
// orientation. Our CT images are usually IS and the masks are SI.
// This routine makes them consistent.
void ZFlipIfNecessary( OrientedInputImageType::Pointer refImage,
                       OrientedMaskImageType::Pointer image );


// routine to read in an image.  The image may be a DICOM series, a
// single DICOM file, or a single file of another format. This
// capability (plus some more) will be supplanted by an
// itk::ImageArchetypeReader.
int ReadArchetypeImage(const std::string &filename,
                       OrientedInputImageType::Pointer &img,
                       std::string &subjectname);

// routine to write out an image. Image can be written as a single
// file or as a series of files. if "byslice" is true, then "filename"
// is expected to be a pattern such as foo-%d.mhd. This capability
// (plus some mroe) will be supplanted by an itk::ImageArchetypeWriter.
int WriteArchetypeImage(const std::string &filename,
                        OrientedOutputImageType *img,
                        bool byslice, bool compress=false); 

#ifndef ITK_MANUAL_INSTANTIATION
#include "IO_Utils_ITK.cxx"
#endif

#endif

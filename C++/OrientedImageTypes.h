#ifndef __OrientedImageTypes_h
#define __OrientedImageTypes_h


#include "itkImage.h"

typedef itk::Image<short, 3 >    OrientedInputImageType;
typedef itk::Image<short, 2 >    OrientedInputImage2DType;
typedef itk::Image<unsigned char, 3>  OrientedMaskImageType;
typedef itk::Image< unsigned char, 2 > OrientedMaskImage2DType;
typedef itk::Image<unsigned char, 3 > OrientedOutputImageType;
typedef itk::Image<float, 3>  OrientedResponseImageType;
typedef itk::Image<float, 3>  OrientedFloatImageType;
typedef itk::Image<float, 2>  OrientedFloatImage2DType;

typedef itk::Image< itk::Vector<double, 3> , 3 > OrientedFloatVectorImageType;


#endif

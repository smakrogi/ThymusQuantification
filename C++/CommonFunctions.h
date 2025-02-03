/*===========================================================================

Program:   Segmentation/Classification of thymus tissues from MRI.
Module:    $RCSfile: CommonFunctions.h,v $
Language:  C++
Date:      $Date: 2010/04/15 10:42:32 $
Version:   $Revision: 0.1 $
Author:    S. K. Makrogiannis
3T MRI Facility National Institute on Aging/National Institutes of Health.

=============================================================================*/

#ifndef _CommonFunctions_h
#define _CommonFunctions_h

#include <cmath>
#include "itkImage.h"
// #include "itkOrientedImage.h"

// Basic operations.
#ifndef MY_ROUND
#define MY_ROUND(X) ( ( X - floor(X) ) >= 0.5? ceil(X):floor(X) )
#endif 

#ifndef MY_SQUARE
#define MY_SQUARE(X)  (X * X)
#endif

#ifdef _WIN32
static const std::string PathSeparator = "\\";
#else
static const std::string PathSeparator = "/";
#endif


// Application parameters.
#define BACKGROUND 0
#define FOREGROUND 255

typedef enum {WS_KMEANS=0, 
	      KMEANS_3D, 
	      REG_GROW_3D, 
	      WS_LEV_SETS,
	      LEV_SETS_3D,
	      WS_GAC,
	      GAC_3D,
	      WS_GAC_SHP,
	      GAC_SHP_3D,
	      WS_CHAN_VESE,
	      GMM_BAYES_3D,
	      KMEANS_MRF_3D} 
  SegmentationMethod;

typedef enum {CONVENTIONAL=0, 
	      PARZEN}
  GradientMethod;

typedef enum {SEED_BASED_ROI=0,
	      FAST_MARCHING}
  LevelSetInitializationMethod;

typedef enum {PARZEN_EDGES=0,
	      GAUSSIAN,
	      DIFFUSION,
	      MEDIAN}
  MultiscaleGeneratorMethod;


#define DEBUG 1

#define HALFROIEDGELENGTH 2.0
#define SMOOTHINGSIGMA 0.0625
#define SMOOTHINGITERATIONS 8
#define MEDIANFILTERRADIUS 1
#define GRADIENTMETHOD 2
#define PARZENRADIUS 1
#define PARZENBANDWIDTH 0.1
#define PARZENBANDWIDTHBASE 1.2
#define CLUSTERS 3
#define LEVELSETSIGMOIDALPHA -4 //-250
#define LEVELSETSIGMOIDBETA 50 //1500
#define SIGMOIDBETAALPHARATIO 4.5
#define FASTMARCHINGSTOPPINGTIME 10.0
#define LEVELSETPROPAGATIONSCALINGFACTOR 1.0
#define LEVELSETCURVATURESCALINGFACTOR 0.1
#define LEVELSETADVECTIONSCALINGFACTOR 2.5
#define LEVELSETMAXITERATIONS 600
#define LEVELSETMAXRMSERROR 0.015
#define LEVELSETSHAPEPRIORSCALINGFACTOR 0.02
#define LEVELSETPCAMODES  3
#define LEVELSETMAPCONTOURWEIGHT 1.0
#define LEVELSETMAPIMAGEGRADIENTWEIGHT 10.0
#define LEVELSETMAPSHAPEWEIGHT 1.0
#define CHANVESECURVATUREWEIGHT 0.0
#define CHANVESEAREAWEIGHT 0.0
#define CHANVESEVOLUMEWEIGHT 0.0
#define CHANVESEVOLUME 0.0
#define CHANVESEOVERLAPWEIGHT 0.0
#define STRUCTUREELEMENTRADIUS 1
#define NSCALES 1
#define MRFITERATIONS 5
#define MRFSMOOTHINGFACTOR 2
#define MRFNEIGHBORHOODRADIUS 1
#define RGMULTIPLIER 0.5
#define RGITERATIONS 3

static const std::string 
GradientMethodStrings[] = { 
  "CONVENTIONAL",
  "PARZEN" };

static const std::string 
SegmentationMethodStrings[] = { 
  "WS_KMEANS",
  "KMEANS_3D",
  "REG_GROW_3D",
  "WS_LEV_SETS",
  "LEV_SETS_3D",
  "WS_GAC",
  "GAC_3D",
  "WS_GAC_SHP",
  "GAC_SHP_3D",
  "WS_CHAN_VESE",
  "GMM_BAYES_3D",
  "KMEANS_MRF_3D" };

static const std::string 
LevelsetInitializationMethodStrings[] = {
  "SEED_BASED_ROI",
  "FAST_MARCHING"};

static const std::string 
MultiscaleGeneratorMethodStrings[] = {"PARZEN_EDGES",
				      "GAUSSIAN",
				      "DIFFUSION",
				      "MEDIAN"};

static const std::string DefaultNSFilename = "Dixon";
static const std::string DefaultMeanShapeFilename = "meanImage.mha";
static const std::string DefaultShapeModeFilenamePattern = "pcImage%d.mha";
static const char* SHAPE_MODEL_PATH_ENV_VAR = "THYMUS_PCA_SHAPE_MODEL_PATH";

typedef struct data{	      
  int     SizeinPixels;      
  float   SizeinPhysicalUnits;
} infoStruct;

const unsigned int ImageDimension = 3;
const unsigned int VoxelComponents = 3;

typedef float FloatPixelType;
typedef unsigned char MaskPixelType;
typedef unsigned int LabelPixelType;
typedef short InputPixelType;
typedef short OutputPixelType;
typedef itk::Vector<FloatPixelType, VoxelComponents> VectorFloatPixelType;

typedef itk::Image<FloatPixelType, ImageDimension> FloatImageType;
typedef itk::Image<MaskPixelType, ImageDimension> MaskImageType;
typedef itk::Image<LabelPixelType, ImageDimension> LabelImageType;
typedef itk::Image<InputPixelType, ImageDimension> InputImageType;
typedef itk::Image<OutputPixelType, ImageDimension> OutputImageType;
typedef itk::Image<VectorFloatPixelType, ImageDimension> VectorFloatImageType;

FloatImageType::RegionType GenerateROI ( FloatImageType::Pointer inputVolume,
					 FloatImageType::PointType seedPoint,
					 std::vector<float> halfROISize );

int DetermineMedianROILabel( MaskImageType::Pointer InputVolume,
			     MaskImageType::RegionType seedbasedROI );

FloatImageType::Pointer SmoothInputVolume( FloatImageType::Pointer inputVolume,
					   int denoisingMethod );

FloatImageType::Pointer
AddDixonImages( FloatImageType::Pointer watersuppressedVolume,
		FloatImageType::Pointer fatsuppressedVolume );


std::string 
ExtractDirectory( const std::string& path );

std::string 
ExtractFilename( const std::string& path );

std::string 
ChangeExtension( const std::string& path, const std::string& ext );


#endif

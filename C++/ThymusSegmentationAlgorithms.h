/*===========================================================================

Program:   Quantification of thymus tissues from MRI.
Module:    $RCSfile: ThymusSegmentationAlgorithms.h,v $
Language:  C++
Date:      $Date: 2010/04/06 10:42:32 $
Version:   $Revision: 0.1 $
Author:    S. K. Makrogiannis
3T MRI Facility National Institute on Aging/National Institutes of Health.

=============================================================================*/
#ifndef _ThymusSegmentationAlgorithms_h
#define _ThymusSegmentationAlgorithms_h

#include "CommonFunctions.h"

/* namespace itk  */
/* {  */
/*   // class to support progress feeback */
/*   class ShowProgressObject */
/*   { */
/*   public: */
/*     ShowProgressObject(LightProcessObject * o) */
/*       {m_Process = o;} */
/*     void ShowProgress() */
/*     {std::cout << "Progress " << m_Process->GetProgress() << std::endl;} */
/*     LightProcessObject::Pointer m_Process; */
/*   }; */


class ThymusSegmentationObject
{

 public:

  // Functions.
  // Default constructor that sets all segmentation parameters.
  ThymusSegmentationObject() { 

    this->m_outputPath = "./";
    this->m_subjectID = "SubjectID";
    this->m_segmentationMethod = REG_GROW_3D;
    
    this->m_gradientMethod = CONVENTIONAL;
    this->m_gradientSigma = SMOOTHINGSIGMA;
    this->m_ParzenRadius = PARZENRADIUS;
    this->m_ParzenBandwidth = PARZENBANDWIDTH;
    this->m_ParzenBandwidthBase = PARZENBANDWIDTHBASE;
    this->m_nClusters = CLUSTERS;
    this->m_structureelementRadius = STRUCTUREELEMENTRADIUS;
    /*this->m_seedPoint = [41 41 27];*/

    this->m_rgMultiplier = RGMULTIPLIER;
    this->m_rgIterations = RGITERATIONS;

    this->m_mrfIterations = MRFITERATIONS;
    this->m_mrfsmoothingFactor = MRFSMOOTHINGFACTOR; // typically between 1 and 5.
    this->m_mrfNeighborhoodRadius = MRFNEIGHBORHOODRADIUS;

    this->m_sigmoidAlpha = LEVELSETSIGMOIDALPHA; // -16
    this->m_sigmoidBeta = LEVELSETSIGMOIDBETA; // 130
    this->m_levelsetInitializationMethod = SEED_BASED_ROI;
    this->m_fastmarchingStoppingTime = (double) FASTMARCHINGSTOPPINGTIME;
    /* this->m_levelsetMethod = "SHAPE_DETECTION"; // "SHAPE_DETECTION","GEOMETRIC_ACTIVE_CONTOURS" */
    this->m_levelsetMaximumRMSError = LEVELSETMAXRMSERROR;
    this->m_levelsetMaximumIterations = LEVELSETMAXITERATIONS;
    this->m_levelsetPropagationScalingFactor = LEVELSETPROPAGATIONSCALINGFACTOR;
    this->m_levelsetCurvatureScalingFactor = LEVELSETCURVATURESCALINGFACTOR;
    this->m_levelsetAdvectionScalingFactor = LEVELSETADVECTIONSCALINGFACTOR;
    this->m_levelsetShapePriorScalingFactor = LEVELSETSHAPEPRIORSCALINGFACTOR;
    this->m_levelsetNumberOfPCAModes =  LEVELSETPCAMODES;
    this->m_levelsetMAPContourWeight = LEVELSETMAPCONTOURWEIGHT;
    this->m_levelsetMAPImageGradientWeight = LEVELSETMAPIMAGEGRADIENTWEIGHT;
    this->m_levelsetMAPShapeWeight = LEVELSETMAPSHAPEWEIGHT;
    this->m_chanveseCurvatureWeight = CHANVESECURVATUREWEIGHT;
    this->m_chanveseAreaWeight = CHANVESEAREAWEIGHT;
    this->m_chanveseVolumeWeight = CHANVESEVOLUMEWEIGHT;
    this->m_chanveseVolume = CHANVESEVOLUME;
    this->m_chanveseOverlapWeight = CHANVESEOVERLAPWEIGHT;
    this->m_multiscaleGenerator = PARZEN_EDGES;
    this->m_nScales = NSCALES;
    this->m_currentScale = 0; 
}


  // Destructor, not much to add because of garbage collection
  // in ITK.
  ~ThymusSegmentationObject() {}


  // Setting the input images.
  void SetInputs( FloatImageType::Pointer nonSuppressedVolume,
		  FloatImageType::Pointer waterSuppressedVolume,
		  FloatImageType::Pointer fatSuppressedVolume );


  // Set algorithm parameters from user's input.
  void SetAlgorithmParameters ( std::string outputPath,
				std::string subjectID,
				FloatImageType::PointType seedPoint,
				FloatImageType::RegionType seedbasedROI,
				int segmentationMethod,
				int nClusters,
				int gradientMethod,
				int parzenRadius,
				float parzenBandwidth,
				float parzenBandwidthBase,
				float structureelementRadius,
				float rgMultiplier,
				int rgIterations,
				int mrfIterations,
				float mrfsmoothingFactor,
				int mrfNeighborhoodRadius,
				float sigmoidAlpha,
				float sigmoidBeta,
				int levelsetInitializationMethod,
				double fastmarchingStoppingTime,
				float levelsetPropagationScalingFactor,
				float levelsetCurvatureScalingFactor,
				float levelsetAdvectionScalingFactor,
				int levelsetMaximumIterations,
				float levelsetMaximumRMSError,
				double levelsetShapePriorScalingFactor,
				int levelsetNumberOfPCAModes,
				double levelsetMAPContourWeight,
				double levelsetMAPImageGradientWeight,
				double levelsetMAPShapeWeight,
				double chanveseCurvatureWeight,
				double chanveseAreaWeight,
				double chanveseVolumeWeight,
				double chanveseVolume,
				double chanveseOverlapWeight,
				int multiscaleGenerator,
				int nScales)
  { this->m_outputPath = outputPath;
    this->m_subjectID = subjectID;
    this->m_seedPoint = seedPoint;
    this->m_seedbasedROI = seedbasedROI; 
    this->m_segmentationMethod = segmentationMethod;
    this->m_nClusters = nClusters;
    this->m_gradientMethod = gradientMethod;
    this->m_ParzenRadius = parzenRadius;
    this->m_ParzenBandwidth = parzenBandwidth;
    this->m_ParzenBandwidthBase = parzenBandwidthBase;
    this->m_structureelementRadius = structureelementRadius; 
    this->m_rgMultiplier = rgMultiplier;
    this->m_rgIterations = rgIterations;
    this->m_mrfIterations = mrfIterations;
    this->m_mrfsmoothingFactor = mrfsmoothingFactor; 
    this->m_mrfNeighborhoodRadius = mrfNeighborhoodRadius;
    this->m_sigmoidAlpha = sigmoidAlpha;
    this->m_sigmoidBeta = sigmoidBeta;
    this->m_levelsetInitializationMethod = levelsetInitializationMethod;
    this->m_fastmarchingStoppingTime = fastmarchingStoppingTime;
    this->m_levelsetPropagationScalingFactor = 
      levelsetPropagationScalingFactor;
    this->m_levelsetCurvatureScalingFactor =
      levelsetCurvatureScalingFactor;
    this->m_levelsetAdvectionScalingFactor = 
      levelsetAdvectionScalingFactor;
    this->m_levelsetMaximumIterations = 
      levelsetMaximumIterations;
    this->m_levelsetMaximumRMSError = 
      levelsetMaximumRMSError;
    this->m_levelsetShapePriorScalingFactor = levelsetShapePriorScalingFactor;
    this->m_levelsetNumberOfPCAModes =  levelsetNumberOfPCAModes;
    this->m_levelsetMAPContourWeight = levelsetMAPContourWeight;
    this->m_levelsetMAPImageGradientWeight = levelsetMAPImageGradientWeight;
    this->m_levelsetMAPShapeWeight = levelsetMAPShapeWeight;
    this->m_chanveseCurvatureWeight = chanveseCurvatureWeight;
    this->m_chanveseAreaWeight = chanveseAreaWeight;
    this->m_chanveseVolumeWeight = chanveseVolumeWeight;
    this->m_chanveseVolume = chanveseVolume;
    this->m_chanveseOverlapWeight = chanveseOverlapWeight;
    this->m_multiscaleGenerator = multiscaleGenerator;
    this->m_nScales = nScales;
  }


  // Set debug mode.
  void SetDebug( bool value ) { m_debug = value; }

  // Return the segmentation mask after execution.
  MaskImageType::Pointer GetOutput() { return m_outputMaskVolume; }

  // Execute selected segmentation algorithm.
  void RunThymusSegmentationAlgorithms();

  // Return volumetric measures.
  infoStruct
    DisplayVolumetrics();

 protected:

  // Foreground/background separation.
  void
    SeparateAirFromOtherTissues();

  // Parzen edge detector.
  FloatImageType::Pointer
    ParzenEdgeMap(FloatImageType::Pointer inputImage);

  // Parzen edge detector on vector pixel values.
  FloatImageType::Pointer 
    VectorParzenEdgeMap(VectorFloatImageType::Pointer inputImage);

  // Form a vector volume from three input images.
  void
    FormNSWSFSVectorVolume();

  // Initialize ROI before level sets by fast marching.
  MaskImageType::Pointer
    InitializeROIbyFastMarching(FloatImageType::Pointer speedImage);

  // Apply post-processing after clustering to pick
  // thymus region.
  void
    PostprocessClustering(MaskImageType::Pointer outputMaskVolume);

  // Apply K-Means on the formed images with vector pixel values.
  MaskImageType::Pointer 
    KMeansonVectorVolumes();

  // Apply ITK's region growing algorithm.
  MaskImageType::Pointer 
    RegionGrowingOnVectorImage();
  
  // Geodesic active contours with shape prior model.
  FloatImageType::Pointer
    GeodesicActiveContoursWithShapePriors(FloatImageType::Pointer ROIDistanceOutput,
					  FloatImageType::Pointer speedImage);

  // All implemented level set pipelines.
  void
    LevelSetBasedSegmentation();

  // Multi-scale scheme for level set segmentation
  void 
    MultiscaleLevelSetBasedSegmentation();

  // Use EM-GMM followed by Bayesian decision making.
  /* MaskImageType::Pointer  */
  /*   EM_GMM_BayesianClassification(); */

  // Apply MRF segmentation to group previously detected clusters
  // using K-Means.
  MaskImageType::Pointer 
    MRFSegmentationAfterClustering(MaskImageType::Pointer LabelImage);

private:
  // Basic image variables.
  FloatImageType::PointType m_seedPoint;
  // Original volume data used
  // for computing volumetics and intensity statistics,
  // such as fat ratios and water ratios.
  FloatImageType::Pointer m_nonsuppressedVolume, 
    m_fatsuppressedVolume,
    m_watersuppressedVolume;
  // Smoothed volume data used in segmentation and
  // scale-space generation.
  FloatImageType::Pointer m_smoothedNonSuppressedVolume,
    m_smoothedWaterSuppressedVolume,
    m_smoothedFatSuppressedVolume;
  VectorFloatImageType::Pointer m_VectorVolume;
  MaskImageType::Pointer m_AirMaskImage;
  bool m_debug;

  std::string m_outputPath, m_subjectID;
  int m_segmentationMethod, m_nClusters;
  float m_structureelementRadius;
  FloatImageType::RegionType m_seedbasedROI;

  // Gradient-related.
  int m_gradientMethod;
  float m_gradientSigma;
  int m_ParzenRadius;
  double m_ParzenBandwidth;
  double m_ParzenBandwidthBase;

  // RG-related.
  float m_rgMultiplier;
  int m_rgIterations;

  // K-means-related.
  std::vector<double> m_kMeansCentroids;
  int m_mrfIterations;
  float m_mrfsmoothingFactor; // typically between 1 and 5.
  int m_mrfNeighborhoodRadius;

  // Level-set-related.
  // std::string m_levelsetMethod;
  int m_levelsetInitializationMethod;
  double m_fastmarchingStoppingTime;
  float m_sigmoidAlpha, m_sigmoidBeta;
  float m_levelsetMaximumRMSError;
  int m_levelsetMaximumIterations;
  float m_levelsetPropagationScalingFactor, 
    m_levelsetCurvatureScalingFactor,
    m_levelsetAdvectionScalingFactor;
  double m_levelsetShapePriorScalingFactor;
  int m_levelsetNumberOfPCAModes;
  double m_levelsetMAPContourWeight,
    m_levelsetMAPImageGradientWeight,
    m_levelsetMAPShapeWeight;

  double m_chanveseCurvatureWeight, 
    m_chanveseAreaWeight, 
    m_chanveseVolumeWeight,
    m_chanveseVolume,
    m_chanveseOverlapWeight;

  int m_nScales, m_currentScale;
  int m_multiscaleGenerator;

  // Output mask image.
  MaskImageType::Pointer m_outputMaskVolume;
};

/* } // namespace itk *\/ */

#endif

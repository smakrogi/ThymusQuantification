/*===========================================================================

Program:   Quantification of thymus tissues from MRI.
Module:    $RCSfile: ThymusQuantification.cxx,v $
Language:  C++
Date:      $Date: 2010/04/24 10:42:32 $
Version:   $Revision: 0.5 $
Author:    S. K. Makrogiannis
3T MRI Facility National Institute on Aging/National Institutes of Health.

=============================================================================*/

#ifndef __VectorParzenEdgeEstimation_h
#define __VectorParzenEdgeEstimation_h

#ifndef MY_SQUARE
#define MY_SQUARE(X)  (X * X)
#endif

// Estimates the probability of edge presence.
// Can be used to estimate "probabilistic" gradient magnitude.
// Arguments: vector image, window width, parzen 

template<class T>
double VectorParzenEdgeEstimation(const T inputArray, int window, double sigma)
{
int i;
double ds, pot, potentialSum=0.0, squaredSigma;
vnl_vector<int>  matrixSize(2, 0);
int nSamples = window*window*window; 
int vectorLength = inputArray.size() / nSamples;
vnl_vector<double> vectorMean(vectorLength, 0.0);
vnl_vector<double> inputVector(vectorLength, 0.0);


  // for( i=0; i<(window*window*window); i++)
  // totalSum += inputVector[i];
  // mean = totalSum / (float)(nSamples);

/*    std::cout << "matrix="   */
/* 	     << inputArray */
/* 	     << std::endl;  */

  for( i=0; i<vectorLength; i++)
   vectorMean[i] = (inputArray.get_column(i)).mean();

  //  Estimate pdf based on Parzen kernels.
  squaredSigma = MY_SQUARE(sigma);
  for(i=0;i<nSamples;i++)
  {
    inputVector = inputArray.get_row(i);
    ds = vnl_vector_ssd(inputVector, vectorMean);
    pot = (float)exp(-ds/(2*squaredSigma))/sigma;
    potentialSum += pot;
/*     std::cout << "input=" << inputVector */
/* 	      << ", mean=" << vectorMean  */
/* 	      << ", ssd=" << dx */
/* 	      << ", potential=" << potentialSum */
/* 	      << std::endl; */
  }
  potentialSum = potentialSum / (float)(nSamples);
  return potentialSum;
}


#endif

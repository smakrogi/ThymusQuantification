/*===========================================================================

Program:   Quantification of thymus tissues from MRI.
Module:    $RCSfile: ThymusQuantification.cxx,v $
Language:  C++
Date:      $Date: 2010/04/23 10:42:32 $
Version:   $Revision: 0.5 $
Author:    S. K. Makrogiannis
3T MRI Facility National Institute on Aging/National Institutes of Health.

=============================================================================*/

#ifndef __ParzenEdgeEstimation_h
#define __ParzenEdgeEstimation_h

#ifndef MY_SQUARE
#define MY_SQUARE(X)  (X * X)
#endif

// Estimates the probability of edge presence.
// Can be used to estimate "probabilistic" gradient magnitude.
// Arguments: vector image, window width, parzen 

template<class T>
double ParzenEdgeEstimation(const T inputVector, int window, double sigma)
{
int i;
double dx, ds, pot, totalSum=0.0, mean, potentialSum=0.0, squaredSigma;

  int nSamples = window*window*window;

  for( i=0; i<(window*window*window); i++)
 	  totalSum += inputVector[i];

  mean = totalSum / (float)(nSamples);

  //  Estimate pdf based on Parzen kernels.

  squaredSigma = MY_SQUARE(sigma);
  for(i=0;i<nSamples;i++)
  {
	  dx=inputVector[i]-mean;
	  ds= MY_SQUARE(dx);
	  pot = (float)exp(-ds/(2*squaredSigma))/sigma;
	  potentialSum+=pot;
  }


  potentialSum = potentialSum / (float)(nSamples);

  return potentialSum;
}


#endif

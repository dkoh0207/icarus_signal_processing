#ifndef __SIGPROC_TOOLS_MORPH1D_CXX__
#define __SIGPROC_TOOLS_MORPH1D_CXX__

#include "Morph1D.h"


void icarus_signal_processing::Morph1D::getWaveformParams(
  const std::vector<short>& waveform,
  float& mean,
  float& median,
  float& rms)
{
  getWaveformParams<short>(waveform, mean, median, rms);
  return;
}

void icarus_signal_processing::Morph1D::getWaveformParams(
  const std::vector<float>& waveform,
  float& mean,
  float& median,
  float& rms)
{
  getWaveformParams<float>(waveform, mean, median, rms);
  return;
}

void icarus_signal_processing::Morph1D::getWaveformParams(
  const std::vector<double>& waveform,
  float& mean,
  float& median,
  float& rms)
{
  getWaveformParams<double>(waveform, mean, median, rms);
  return;
}

template <typename T> 
void icarus_signal_processing::Morph1D::getWaveformParams(
  const std::vector<T>& waveform,
  float& mean,
  float& median,
  float& rms)
{
  /*
  Calculate waveform parameters for a given 1D waveform.

  INPUTS:
    - waveform: 1D Pedestal Corrected Waveform.
  
  MODIFIES:
    - mean: computes mean of input std::vector waveform.
    - median: computes median of input waveform.
    - rms: computes RMS value of median subtracted waveform. 
  */
  typename std::vector<T> localWaveform = waveform;

  float realMean(
    float(std::accumulate(localWaveform.begin(),
    localWaveform.end(),0)) / float(localWaveform.size()));

  if (localWaveform.size() % 2 == 0) {
    const auto m1 = localWaveform.begin() + localWaveform.size() / 2 - 1;
    const auto m2 = localWaveform.begin() + localWaveform.size() / 2;
    std::nth_element(localWaveform.begin(), m1, localWaveform.end());
    std::nth_element(localWaveform.begin(), m2, localWaveform.end());
    median = (*m1 + *m2) / 2.0;
  } else {
    median = localWaveform[localWaveform.size()/2];
  }

  mean = realMean;
  std::vector<float> adcLessPedVec;
  adcLessPedVec.resize(localWaveform.size());
  std::transform(
    localWaveform.begin(),localWaveform.end(),
    adcLessPedVec.begin(),std::bind(
      std::minus<float>(),std::placeholders::_1,median));
  rms = std::sqrt(std::inner_product(
    adcLessPedVec.begin(), adcLessPedVec.end(), adcLessPedVec.begin(), 0.)
    / float(adcLessPedVec.size()));
  return;
}


void icarus_signal_processing::Morph1D::getDilation(
  const Waveform<short>& waveform,
  const unsigned int structuringElement,
  Waveform<short>& dilationVec) const
{
  getDilation<short>(waveform, structuringElement, dilationVec);
  return;
}

void icarus_signal_processing::Morph1D::getDilation(
  const Waveform<float>& waveform,
  const unsigned int structuringElement,
  Waveform<float>& dilationVec) const
{
  getDilation<float>(waveform, structuringElement, dilationVec);
  return;
}

void icarus_signal_processing::Morph1D::getDilation(
  const Waveform<double>& waveform,
  const unsigned int structuringElement,
  Waveform<double>& dilationVec) const
{
  getDilation<double>(waveform, structuringElement, dilationVec);
  return;
}

template <typename T> 
void icarus_signal_processing::Morph1D::getDilation(
  const Waveform<T>& inputWaveform,
  const unsigned int structuringElement,
  Waveform<T>& dilationVec) const
{
  /*
  Module for 1D Dilation Filter.

  INPUTS:
    - waveform: 1D Pedestal Corrected Waveform.
    - structuringElement: Size of moving window
  
  MODIFIES:
    - dilationVec: Returned Dilation Vector.
  */
  // Set the window size
  int halfWindowSize(structuringElement/2);
  // Initialize min and max elements
  std::pair<typename Waveform<T>::const_iterator,
            typename Waveform<T>::const_iterator> minMaxItr =
            std::minmax_element(
              inputWaveform.begin(),inputWaveform.begin()+halfWindowSize);

  typename Waveform<T>::const_iterator minElementItr = minMaxItr.first;
  typename Waveform<T>::const_iterator maxElementItr = minMaxItr.second;

  // Initialize the erosion and dilation vectors
  dilationVec.resize(inputWaveform.size());
  // Now loop through remaining elements and complete the vectors
  typename Waveform<T>::iterator maxItr = dilationVec.begin();
  for (typename Waveform<T>::const_iterator inputItr = 
    inputWaveform.begin(); inputItr != inputWaveform.end(); inputItr++)
  {
    // There are two conditions to check:
    // 1) is the current min/max element outside the current window?
    // 2) is the new element smaller/larger than the current min/max?
    // Make sure we are not running off the end of the vector
    if (std::distance(inputItr,inputWaveform.end()) > halfWindowSize)
    {
      if (std::distance(minElementItr,inputItr) >= halfWindowSize)
          minElementItr = std::min_element(
            inputItr - halfWindowSize, inputItr + halfWindowSize + 1);
      else if (*(inputItr + halfWindowSize) < *minElementItr)
          minElementItr = inputItr + halfWindowSize;
      if (std::distance(maxElementItr,inputItr) >= halfWindowSize)
          maxElementItr = std::max_element(
            inputItr - halfWindowSize, inputItr + halfWindowSize + 1);
      else if (*(inputItr + halfWindowSize) > *maxElementItr)
          maxElementItr = inputItr + halfWindowSize;
    }
    // Update the vectors
    *maxItr++ = *maxElementItr;
  }
  return;
}


void icarus_signal_processing::Morph1D::getErosion(
  const Waveform<short>& waveform,
  const unsigned int structuringElement,
  Waveform<short>& erosionVec) const
{
  getErosion<short>(waveform, structuringElement, erosionVec);
  return;
}

void icarus_signal_processing::Morph1D::getErosion(
  const Waveform<float>& waveform,
  const unsigned int structuringElement,
  Waveform<float>& erosionVec) const
{
  getErosion<float>(waveform, structuringElement, erosionVec);
  return;
}

void icarus_signal_processing::Morph1D::getErosion(
  const Waveform<double>& waveform,
  const unsigned int structuringElement,
  Waveform<double>& erosionVec) const
{
  getErosion<double>(waveform, structuringElement, erosionVec);
  return;
}

template <typename T> 
void icarus_signal_processing::Morph1D::getErosion(
  const Waveform<T>& inputWaveform,
  const unsigned int structuringElement,
  Waveform<T>& erosionVec) const
{
  // Set the window size
  int halfWindowSize(structuringElement/2);
  // Initialize min and max elements
  std::pair<typename Waveform<T>::const_iterator,
            typename Waveform<T>::const_iterator> minMaxItr =
            std::minmax_element(
              inputWaveform.begin(),inputWaveform.begin()+halfWindowSize);

  typename Waveform<T>::const_iterator minElementItr = minMaxItr.first;
  typename Waveform<T>::const_iterator maxElementItr = minMaxItr.second;

  // Initialize the erosion and dilation vectors
  erosionVec.resize(inputWaveform.size());
  // Now loop through remaining elements and complete the vectors
  typename Waveform<T>::iterator minItr = erosionVec.begin();

  for (typename Waveform<T>::const_iterator inputItr = 
    inputWaveform.begin(); inputItr != inputWaveform.end(); inputItr++)
  {
    // There are two conditions to check:
    // 1) is the current min/max element outside the current window?
    // 2) is the new element smaller/larger than the current min/max?
    // Make sure we are not running off the end of the vector
    if (std::distance(inputItr,inputWaveform.end()) > halfWindowSize)
    {
      if (std::distance(minElementItr,inputItr) >= halfWindowSize)
          minElementItr = std::min_element(
            inputItr - halfWindowSize, inputItr + halfWindowSize + 1);
      else if (*(inputItr + halfWindowSize) < *minElementItr)
          minElementItr = inputItr + halfWindowSize;
      if (std::distance(maxElementItr,inputItr) >= halfWindowSize)
          maxElementItr = std::max_element(
            inputItr - halfWindowSize, inputItr + halfWindowSize + 1);
      else if (*(inputItr + halfWindowSize) > *maxElementItr)
          maxElementItr = inputItr + halfWindowSize;
    }
    // Update the vectors
    *minItr++ = *minElementItr;
  }
  return;
}


void icarus_signal_processing::Morph1D::getGradient(
  const Waveform<short>& waveform,
  const unsigned int structuringElement,
  Waveform<short>& gradientVec) const
{
  getGradient<short>(waveform, structuringElement, gradientVec);
  return;
}

void icarus_signal_processing::Morph1D::getGradient(
  const Waveform<float>& waveform,
  const unsigned int structuringElement,
  Waveform<float>& gradientVec) const
{
  getGradient<float>(waveform, structuringElement, gradientVec);
  return;
}

void icarus_signal_processing::Morph1D::getGradient(
  const Waveform<double>& waveform,
  const unsigned int structuringElement,
  Waveform<double>& gradientVec) const
{
  getGradient<double>(waveform, structuringElement, gradientVec);
  return;
}

template <typename T> 
void icarus_signal_processing::Morph1D::getGradient(
  const Waveform<T>& inputWaveform,
  const unsigned int structuringElement,
  Waveform<T>& gradientVec) const
{
  // Set the window size
  int halfWindowSize(structuringElement/2);
  // Initialize min and max elements
  std::pair<typename Waveform<T>::const_iterator,
            typename Waveform<T>::const_iterator> minMaxItr =
            std::minmax_element(
              inputWaveform.begin(),inputWaveform.begin()+halfWindowSize);

  typename Waveform<T>::const_iterator minElementItr = minMaxItr.first;
  typename Waveform<T>::const_iterator maxElementItr = minMaxItr.second;

  // Initialize the erosion and dilation vectors
  gradientVec.resize(inputWaveform.size());
  // Now loop through remaining elements and complete the vectors
  typename Waveform<T>::iterator difItr = gradientVec.begin();

  for (typename Waveform<T>::const_iterator inputItr = 
    inputWaveform.begin(); inputItr != inputWaveform.end(); inputItr++)
  {
    // There are two conditions to check:
    // 1) is the current min/max element outside the current window?
    // 2) is the new element smaller/larger than the current min/max?
    // Make sure we are not running off the end of the vector
    if (std::distance(inputItr,inputWaveform.end()) > halfWindowSize)
    {
      if (std::distance(minElementItr,inputItr) >= halfWindowSize)
          minElementItr = std::min_element(
            inputItr - halfWindowSize, inputItr + halfWindowSize + 1);
      else if (*(inputItr + halfWindowSize) < *minElementItr)
          minElementItr = inputItr + halfWindowSize;
      if (std::distance(maxElementItr,inputItr) >= halfWindowSize)
          maxElementItr = std::max_element(
            inputItr - halfWindowSize, inputItr + halfWindowSize + 1);
      else if (*(inputItr + halfWindowSize) > *maxElementItr)
          maxElementItr = inputItr + halfWindowSize;
    }
    // Update the vectors
    *difItr++ = *maxElementItr - *minElementItr;
  }
  return;
}


void icarus_signal_processing::Morph1D::getAverage(
  const Waveform<short>& waveform,
  const unsigned int structuringElement,
  Waveform<short>& averageVec) const
{
  getAverage<short>(waveform, structuringElement, averageVec);
  return;
}

void icarus_signal_processing::Morph1D::getAverage(
  const Waveform<float>& waveform,
  const unsigned int structuringElement,
  Waveform<float>& averageVec) const
{
  getAverage<float>(waveform, structuringElement, averageVec);
  return;
}

void icarus_signal_processing::Morph1D::getAverage(
  const Waveform<double>& waveform,
  const unsigned int structuringElement,
  Waveform<double>& averageVec) const
{
  getAverage<double>(waveform, structuringElement, averageVec);
  return;
}

template <typename T> 
void icarus_signal_processing::Morph1D::getAverage(
  const Waveform<T>& inputWaveform,
  const unsigned int structuringElement,
  Waveform<T>& averageVec) const
{
  // Set the window size
  int halfWindowSize(structuringElement/2);
  // Initialize min and max elements
  std::pair<typename Waveform<T>::const_iterator,
            typename Waveform<T>::const_iterator> minMaxItr =
            std::minmax_element(
              inputWaveform.begin(),inputWaveform.begin()+halfWindowSize);

  typename Waveform<T>::const_iterator minElementItr = minMaxItr.first;
  typename Waveform<T>::const_iterator maxElementItr = minMaxItr.second;

  // Initialize the erosion and dilation vectors
  averageVec.resize(inputWaveform.size());
  // Now loop through remaining elements and complete the vectors
  typename Waveform<T>::iterator avgItr = averageVec.begin();

  for (typename Waveform<T>::const_iterator inputItr = 
    inputWaveform.begin(); inputItr != inputWaveform.end(); inputItr++)
  {
    // There are two conditions to check:
    // 1) is the current min/max element outside the current window?
    // 2) is the new element smaller/larger than the current min/max?
    // Make sure we are not running off the end of the vector
    if (std::distance(inputItr,inputWaveform.end()) > halfWindowSize)
    {
      if (std::distance(minElementItr,inputItr) >= halfWindowSize)
          minElementItr = std::min_element(
            inputItr - halfWindowSize, inputItr + halfWindowSize + 1);
      else if (*(inputItr + halfWindowSize) < *minElementItr)
          minElementItr = inputItr + halfWindowSize;
      if (std::distance(maxElementItr,inputItr) >= halfWindowSize)
          maxElementItr = std::max_element(
            inputItr - halfWindowSize, inputItr + halfWindowSize + 1);
      else if (*(inputItr + halfWindowSize) > *maxElementItr)
          maxElementItr = inputItr + halfWindowSize;
    }
    // Update the vectors
    *avgItr++ = 0.5 * (*maxElementItr + *minElementItr);
  }
  return;
}


void icarus_signal_processing::Morph1D::getMedian(
  const Waveform<short>& waveform,
  const unsigned int structuringElement,
  Waveform<short>& medianVec) const
{
  getMedian<short>(waveform, structuringElement, medianVec);
  return;
}

void icarus_signal_processing::Morph1D::getMedian(
  const Waveform<float>& waveform,
  const unsigned int structuringElement,
  Waveform<float>& medianVec) const
{
  getMedian<float>(waveform, structuringElement, medianVec);
  return;
}

void icarus_signal_processing::Morph1D::getMedian(
  const Waveform<double>& waveform,
  const unsigned int structuringElement,
  Waveform<double>& medianVec) const
{
  getMedian<double>(waveform, structuringElement, medianVec);
  return;
}

template <typename T> 
void icarus_signal_processing::Morph1D::getMedian(
  const Waveform<T>& inputWaveform,
  const unsigned int structuringElement,
  Waveform<T>& medianVec) const
{
  // Set the window size
  int halfWindowSize(structuringElement/2);
  medianVec.resize(inputWaveform.size());
  size_t nTicks = inputWaveform.size();

  T median = 0;

  for (size_t i=0; i<nTicks; ++i) {
    typename std::vector<T> localVec;
    localVec.reserve(structuringElement);
    int lbx = i - (int) halfWindowSize;
    int ubx = i + (int) halfWindowSize;
    size_t lowerBoundx = std::max(lbx, 0);
    size_t upperBoundx = std::min(ubx, (int) nTicks);
    for (size_t j=lowerBoundx; j<upperBoundx+1; ++j) {
      localVec.push_back(inputWaveform[j]);
    }
    if (structuringElement % 2 == 0) {
      const auto m1 = localVec.begin() + localVec.size() / 2 - 1;
      const auto m2 = localVec.begin() + localVec.size() / 2;
      std::nth_element(localVec.begin(), m1, localVec.end());
      std::nth_element(localVec.begin(), m2, localVec.end());
      median = (*m1 + *m2) / 2.0;
    } else {
      median = localVec[localVec.size() / 2];
    }
    medianVec[i] = median;
  }
  return;
}


void icarus_signal_processing::Morph1D::getOpeningAndClosing(
  const Waveform<short>& inputWaveform,
  const unsigned int structuringElement,
  Waveform<short>& openingVec,
  Waveform<short>& closingVec) const
{
  getOpeningAndClosing<short>(
    inputWaveform, structuringElement, openingVec, closingVec);
  return;
}

void icarus_signal_processing::Morph1D::getOpeningAndClosing(
  const Waveform<float>& inputWaveform,
  const unsigned int structuringElement,
  Waveform<float>& openingVec,
  Waveform<float>& closingVec) const
{
  getOpeningAndClosing<float>(
    inputWaveform, structuringElement, openingVec, closingVec);
  return;
}

void icarus_signal_processing::Morph1D::getOpeningAndClosing(
  const Waveform<double>& inputWaveform,
  const unsigned int structuringElement,
  Waveform<double>& openingVec,
  Waveform<double>& closingVec) const
{
  getOpeningAndClosing<double>(
    inputWaveform, structuringElement, openingVec, closingVec);
  return;
}


template <typename T> 
void icarus_signal_processing::Morph1D::getOpeningAndClosing(
  const Waveform<T>& inputWaveform,
  const unsigned int structuringElement,
  Waveform<T>& openingVec,
  Waveform<T>& closingVec) const
{
  Waveform<T> dilationVec;
  Waveform<T> erosionVec;

  getDilation(inputWaveform, structuringElement, dilationVec);
  getErosion(inputWaveform, structuringElement, erosionVec);
  // Set the window size
  int halfWindowSize(structuringElement/2);
  // Start with the opening: get the max element in the input erosion vector
  typename Waveform<T>::iterator maxElementItr = 
    std::max_element(erosionVec.begin(),erosionVec.begin()+halfWindowSize);
  // Initialize the opening vector
  openingVec.resize(erosionVec.size());
  // Now loop through remaining elements and complete the vectors
  typename Waveform<T>::iterator maxItr = openingVec.begin();
  for (typename Waveform<T>::iterator inputItr = erosionVec.begin(); 
    inputItr != erosionVec.end(); inputItr++)
  {
    // There are two conditions to check:
    // 1) is the current min/max element outside the current window?
    // 2) is the new element smaller/larger than the current min/max?
    // Make sure we are not running off the end of the vector
    if (std::distance(inputItr,erosionVec.end()) > halfWindowSize)
    {
      if (std::distance(maxElementItr,inputItr) >= halfWindowSize)
          maxElementItr = std::max_element(
            inputItr - halfWindowSize, inputItr + halfWindowSize + 1);
      else if (*(inputItr + halfWindowSize) > *maxElementItr)
          maxElementItr = inputItr + halfWindowSize;
    }
    // Update the vectors
    *maxItr++ = *maxElementItr;
  }
  // Now go with the closing: get the min element in the input dilation vector
  typename Waveform<T>::iterator minElementItr = std::min_element(
    dilationVec.begin(),dilationVec.begin()+halfWindowSize);
  // Initialize the opening and closing vectors
  closingVec.resize(dilationVec.size());
  // Now loop through remaining elements and complete the vectors
  typename Waveform<T>::iterator minItr = closingVec.begin();
  for (typename Waveform<T>::iterator inputItr = dilationVec.begin(); 
    inputItr != dilationVec.end(); inputItr++)
  {
    // There are two conditions to check:
    // 1) is the current min/max element outside the current window?
    // 2) is the new element smaller/larger than the current min/max?
    // Make sure we are not running off the end of the vector
    if (std::distance(inputItr,dilationVec.end()) > halfWindowSize)
    {
      if (std::distance(minElementItr,inputItr) >= halfWindowSize)
          minElementItr = std::min_element(
            inputItr - halfWindowSize, inputItr + halfWindowSize + 1);
      else if (*(inputItr + halfWindowSize) < *minElementItr)
          minElementItr = inputItr + halfWindowSize;
    }
    // Update the vectors
    *minItr++ = *minElementItr;
  }
  return;
}


#endif

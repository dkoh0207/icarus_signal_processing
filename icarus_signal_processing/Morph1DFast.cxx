#ifndef __SIGPROC_TOOLS_MORPH1DFAST_CXX__
#define __SIGPROC_TOOLS_MORPH1DFAST_CXX__

#include "Morph1DFast.h"


void sigproc_tools::Morph1DFast::getDilation(
  const Waveform<short>& waveform,
  const unsigned int structuringElement,
  Waveform<short>& dilationVec) const
{
  getDilation<short>(waveform, structuringElement, dilationVec);
  return;
}

void sigproc_tools::Morph1DFast::getDilation(
  const Waveform<float>& waveform,
  const unsigned int structuringElement,
  Waveform<float>& dilationVec) const
{
  getDilation<float>(waveform, structuringElement, dilationVec);
  return;
}

void sigproc_tools::Morph1DFast::getDilation(
  const Waveform<double>& waveform,
  const unsigned int structuringElement,
  Waveform<double>& dilationVec) const
{
  getDilation<double>(waveform, structuringElement, dilationVec);
  return;
}

template <typename T>
void sigproc_tools::Morph1DFast::getDilation(
  const Waveform<T>& inputWaveform,
  const unsigned int structuringElement,
  Waveform<T>& dilationVec) const
{
  size_t N = inputWaveform.size();
  size_t k = (size_t) structuringElement;
  if (N <= k) {
    std::cout << "Input array size " << N << " must be greater than structuring element size " << k << std::endl;
    return;
  }
  size_t bufferSize = N + 2 * (k/2) + (k - (N % k));
  size_t windowSize = k/2;
  size_t paddingSize = k - (N % k);
  std::vector<T> suffixArr(bufferSize);
  std::vector<T> prefixArr(bufferSize);
  dilationVec.resize(N);

  // Padding Operations on Buffers
  for (size_t i=0; i<windowSize; ++i) {
    suffixArr[i] = std::numeric_limits<T>::min();
    prefixArr[i] = std::numeric_limits<T>::min();
  }

  for (size_t i=N+windowSize; i<bufferSize; ++i) {
    suffixArr[i] = std::numeric_limits<T>::min();
    prefixArr[i] = std::numeric_limits<T>::min();
  }

  // Compute Prefix and Suffix Buffers
  for (size_t i=0; i<N+paddingSize; ++i) {
    if (i % k == 0) {
      prefixArr[i+windowSize] = inputWaveform[i];
    } else {
      prefixArr[i+windowSize] = std::max(prefixArr[i+windowSize-1], inputWaveform[i]);
    }
  }

  for (size_t i=N+paddingSize; i!=0; --i) {
    if (i > N) {
      // Compensate for divisibility padding (must be -inf)
      continue;
    }
    else if (i % k == 0) {
      suffixArr[i+windowSize-1] = inputWaveform[i-1];
    } 
    else {
      suffixArr[i+windowSize-1] = std::max(suffixArr[i+windowSize], inputWaveform[i-1]);
    }
  }

  int prefixIndex = 0;
  int suffixIndex = 0;

  for (size_t i=0; i<N; ++i) {
    prefixIndex = i + 2 * windowSize;
    suffixIndex = i;
    dilationVec[i] = std::max(prefixArr[prefixIndex],
      suffixArr[suffixIndex]);
  }
  return;
}


void sigproc_tools::Morph1DFast::getErosion(
  const Waveform<short>& waveform,
  const unsigned int structuringElement,
  Waveform<short>& erosionVec) const
{
  getErosion<short>(waveform, structuringElement, erosionVec);
  return;
}

void sigproc_tools::Morph1DFast::getErosion(
  const Waveform<float>& waveform,
  const unsigned int structuringElement,
  Waveform<float>& erosionVec) const
{
  getErosion<float>(waveform, structuringElement, erosionVec);
  return;
}

void sigproc_tools::Morph1DFast::getErosion(
  const Waveform<double>& waveform,
  const unsigned int structuringElement,
  Waveform<double>& erosionVec) const
{
  getErosion<double>(waveform, structuringElement, erosionVec);
  return;
}

template <typename T>
void sigproc_tools::Morph1DFast::getErosion(
  const Waveform<T>& inputWaveform,
  const unsigned int structuringElement,
  Waveform<T>& erosionVec) const
{
  size_t N = inputWaveform.size();
  size_t k = (size_t) structuringElement;
  if (N <= k) {
    std::cout << "Input array size " << N << " must be greater than structuring element size " << k << std::endl;
    return;
  }
  size_t bufferSize = N + 2 * (k/2) + (k - (N % k));
  size_t windowSize = k/2;
  size_t paddingSize = k - (N % k);
  std::vector<T> suffixArr(bufferSize);
  std::vector<T> prefixArr(bufferSize);
  erosionVec.resize(N);

  // Padding Operations on Buffers
  for (size_t i=0; i<windowSize; ++i) {
    suffixArr[i] = std::numeric_limits<T>::max();
    prefixArr[i] = std::numeric_limits<T>::max();
  }

  for (size_t i=N+windowSize; i<bufferSize; ++i) {
    suffixArr[i] = std::numeric_limits<T>::max();
    prefixArr[i] = std::numeric_limits<T>::max();
  }

  // Compute Prefix and Suffix Buffers
  for (size_t i=0; i<N+paddingSize; ++i) {
    if (i % k == 0) {
      prefixArr[i+windowSize] = inputWaveform[i];
    } else {
      prefixArr[i+windowSize] = std::min(prefixArr[i+windowSize-1], inputWaveform[i]);
    }
  }

  for (size_t i=N+paddingSize; i!=0; --i) {
    if (i > N) {
      // Compensate for divisibility padding (must be -inf)
      continue;
    }
    else if (i % k == 0) {
      suffixArr[i+windowSize-1] = inputWaveform[i-1];
    } 
    else {
      suffixArr[i+windowSize-1] = std::min(suffixArr[i+windowSize], inputWaveform[i-1]);
    }
  }

  int prefixIndex = 0;
  int suffixIndex = 0;

  for (size_t i=0; i<N; ++i) {
    prefixIndex = i + 2 * windowSize;
    suffixIndex = i;
    erosionVec[i] = std::min(prefixArr[prefixIndex],
      suffixArr[suffixIndex]);
  }
  return;
}


void sigproc_tools::Morph1DFast::getGradient(
  const Waveform<short>& waveform,
  const unsigned int structuringElement,
  Waveform<short>& gradientVec) const
{
  getGradient<short>(waveform, structuringElement, gradientVec);
  return;
}

void sigproc_tools::Morph1DFast::getGradient(
  const Waveform<float>& waveform,
  const unsigned int structuringElement,
  Waveform<float>& gradientVec) const
{
  getGradient<float>(waveform, structuringElement, gradientVec);
  return;
}

void sigproc_tools::Morph1DFast::getGradient(
  const Waveform<double>& waveform,
  const unsigned int structuringElement,
  Waveform<double>& gradientVec) const
{
  getGradient<double>(waveform, structuringElement, gradientVec);
  return;
}

template <typename T>
void sigproc_tools::Morph1DFast::getGradient(
  const Waveform<T>& inputWaveform,
  const unsigned int structuringElement,
  Waveform<T>& gradientVec) const
{
  std::vector<T> erosionVec;
  std::vector<T> dilationVec;
  getDilation<T>(inputWaveform, structuringElement, dilationVec);
  getErosion<T>(inputWaveform, structuringElement, erosionVec);

  size_t N = inputWaveform.size();
  gradientVec.resize(N);
  for (size_t i=0; i<N; ++i) {
    gradientVec[i] = dilationVec[i] - erosionVec[i];
  }
  return;
}


void sigproc_tools::Morph1DFast::getOpening(
  const Waveform<short>& waveform,
  const unsigned int structuringElement,
  Waveform<short>& openingVec) const
{
  getOpening<short>(waveform, structuringElement, openingVec);
  return;
}

void sigproc_tools::Morph1DFast::getOpening(
  const Waveform<float>& waveform,
  const unsigned int structuringElement,
  Waveform<float>& openingVec) const
{
  getOpening<float>(waveform, structuringElement, openingVec);
  return;
}

void sigproc_tools::Morph1DFast::getOpening(
  const Waveform<double>& waveform,
  const unsigned int structuringElement,
  Waveform<double>& openingVec) const
{
  getOpening<double>(waveform, structuringElement, openingVec);
  return;
}

template <typename T>
void sigproc_tools::Morph1DFast::getOpening(
  const Waveform<T>& inputWaveform,
  const unsigned int structuringElement,
  Waveform<T>& openingVec) const
{
  std::vector<T> tempVec;
  getErosion<T>(inputWaveform, structuringElement, tempVec);
  getDilation<T>(tempVec, structuringElement, openingVec);
  return;
}


void sigproc_tools::Morph1DFast::getClosing(
  const Waveform<short>& waveform,
  const unsigned int structuringElement,
  Waveform<short>& closingVec) const
{
  getClosing<short>(waveform, structuringElement, closingVec);
  return;
}

void sigproc_tools::Morph1DFast::getClosing(
  const Waveform<float>& waveform,
  const unsigned int structuringElement,
  Waveform<float>& closingVec) const
{
  getClosing<float>(waveform, structuringElement, closingVec);
  return;
}

void sigproc_tools::Morph1DFast::getClosing(
  const Waveform<double>& waveform,
  const unsigned int structuringElement,
  Waveform<double>& closingVec) const
{
  getClosing<double>(waveform, structuringElement, closingVec);
  return;
}

template <typename T>
void sigproc_tools::Morph1DFast::getClosing(
  const Waveform<T>& inputWaveform,
  const unsigned int structuringElement,
  Waveform<T>& closingVec) const
{
  std::vector<T> tempVec;
  getDilation<T>(inputWaveform, structuringElement, tempVec);
  getErosion<T>(tempVec, structuringElement, closingVec);
  return;
}

#endif

#ifndef __SIGPROC_TOOLS_DECONVOLUTION_CXX__
#define __SIGPROC_TOOLS_DECONVOLUTION_CXX__

#include "Deconvolution.h"

// 1D Inverse Filtering. Probably we should not use it...

void sigproc_tools::Deconvolution::Inverse1D(
  std::vector<std::vector<float>>& outputWaveform,
  const std::vector<std::vector<float>>& inputWaveform,
  const std::vector<float>& responseFunction)
{
  Inverse1D<float>(outputWaveform, inputWaveform, responseFunction);
}

void sigproc_tools::Deconvolution::Inverse1D(
  std::vector<std::vector<double>>& outputWaveform,
  const std::vector<std::vector<double>>& inputWaveform,
  const std::vector<double>& responseFunction)
{
  Inverse1D<double>(outputWaveform, inputWaveform, responseFunction);
}

template <typename T>
void sigproc_tools::Deconvolution::Inverse1D(
  std::vector<std::vector<T>>& outputWaveform,
  const std::vector<std::vector<T>>& inputWaveform,
  const std::vector<T>& responseFunction)
{
  size_t numChannels = inputWaveform.size();
  size_t nTicks = inputWaveform.at(0).size();
  outputWaveform.resize(numChannels);
  for (size_t i=0; i<numChannels; ++i) {
    outputWaveform[i].resize(nTicks);
  }
  for (size_t i=0; i<numChannels; ++i) {
    Eigen::FFT<T> fft;
    fft.SetFlag(fft.HalfSpectrum);
    std::vector<std::complex<T>> freqVec;
    std::vector<std::complex<T>> responseFFT;
    fft.fwd(freqVec, inputWaveform[i]);
    fft.fwd(responseFFT, responseFunction);
    for (size_t j=0; j<nTicks; ++j) {
      freqVec[j] = freqVec[j] / responseFFT[j];
    }
    fft.inv(outputWaveform[i], freqVec);
    fft.ClearFlag(fft.HalfSpectrum);
  }
  return;
}


// 1D Wiener Deconvolution.

void sigproc_tools::Deconvolution::Wiener1D(
  std::vector<std::vector<float>>& outputWaveform,
  const std::vector<std::vector<float>>& inputWaveform,
  const std::vector<float>& responseFunction,
  const float noiseVar)
{
  Wiener1D<float>(outputWaveform, inputWaveform, responseFunction, noiseVar);
}

void sigproc_tools::Deconvolution::Wiener1D(
  std::vector<std::vector<double>>& outputWaveform,
  const std::vector<std::vector<double>>& inputWaveform,
  const std::vector<double>& responseFunction,
  const float noiseVar)
{
  Wiener1D<double>(outputWaveform, inputWaveform, responseFunction, noiseVar);
}

template <typename T>
void sigproc_tools::Deconvolution::Wiener1D(
  std::vector<std::vector<T>>& outputWaveform,
  const std::vector<std::vector<T>>& inputWaveform,
  const std::vector<T>& responseFunction,
  const float noiseVar)
{
  size_t numChannels = inputWaveform.size();
  size_t nTicks = inputWaveform.at(0).size();
  std::complex<T> noisePower(noiseVar, 0);

  outputWaveform.resize(numChannels);
  for (size_t i=0; i<numChannels; ++i) {
    outputWaveform[i].resize(nTicks);
  }

  for (size_t i=0; i<numChannels; ++i) {
    Eigen::FFT<T> fft;
    fft.SetFlag(fft.HalfSpectrum);
    std::complex<T> wiener;
    std::vector<std::complex<T>> freqVec;
    std::vector<std::complex<T>> responseFFT;
    fft.fwd(freqVec, inputWaveform[i]);
    fft.fwd(responseFFT, responseFunction);
    for (size_t j=0; j<responseFFT.size(); ++j) {
      wiener = std::conj(responseFFT[j]) / 
      ((std::complex<T>) std::pow(std::abs(responseFFT[j]), 2.0) + 
        (std::complex<T>) noiseVar / (std::complex<T>) std::pow(std::abs(freqVec[j]), 2.0));
      freqVec[j] = freqVec[j] * wiener;
    }
    fft.inv(outputWaveform[i], freqVec);
    fft.ClearFlag(fft.HalfSpectrum);
  }
  return;
}

// 1D Adaptive Wiener Deconvolution

void sigproc_tools::Deconvolution::AdaptiveWiener1D(
  std::vector<std::vector<float>>& outputWaveform,
  const std::vector<std::vector<float>>& inputWaveform,
  const std::vector<float>& responseFunction,
  const std::vector<std::vector<bool>>& selectVals)
{
  AdaptiveWiener1D<float>(outputWaveform, inputWaveform, responseFunction, selectVals);
}

void sigproc_tools::Deconvolution::AdaptiveWiener1D(
  std::vector<std::vector<double>>& outputWaveform,
  const std::vector<std::vector<double>>& inputWaveform,
  const std::vector<double>& responseFunction,
  const std::vector<std::vector<bool>>& selectVals)
{
  AdaptiveWiener1D<double>(outputWaveform, inputWaveform, responseFunction, selectVals);
}

template <typename T>
void sigproc_tools::Deconvolution::AdaptiveWiener1D(
  std::vector<std::vector<T>>& outputWaveform,
  const std::vector<std::vector<T>>& inputWaveform,
  const std::vector<T>& responseFunction,
  const std::vector<std::vector<bool>>& selectVals)
{
  size_t numChannels = inputWaveform.size();
  size_t nTicks = inputWaveform.at(0).size();

  sigproc_tools::MiscUtils utils;

  outputWaveform.resize(numChannels);
  for (size_t i=0; i<numChannels; ++i) {
    outputWaveform[i].resize(nTicks);
  }

  std::vector<T> additiveNoise;
  additiveNoise.reserve(numChannels * nTicks);
  std::vector<T> additiveNoiseSq;
  additiveNoiseSq.reserve(numChannels * nTicks);
  // Estimate noise power
  for (size_t i=0; i<numChannels; ++i) {
    for (size_t j=0; j<nTicks; ++j) {
      // For each center pixel, compute 2D local estimate of noise.
      if (!selectVals[i][j]) {
        additiveNoise.push_back(inputWaveform[i][j]);
        additiveNoiseSq.push_back(std::pow(inputWaveform[i][j], 2.0));
      }
    }
  }

  T noiseMean = std::accumulate(additiveNoise.begin(), 
    additiveNoise.end(), 0.0) / additiveNoise.size();
  T noiseMeanSq = std::accumulate(additiveNoiseSq.begin(), 
    additiveNoiseSq.end(), 0.0) / additiveNoiseSq.size();
  std::complex<T> noiseVar = noiseMeanSq - std::pow(noiseMean, 2.0);

  // Wiener Deconvolution for each channel 1-dimensionally
  for (size_t i=0; i<numChannels; ++i) {
    Eigen::FFT<T> fft;
    fft.SetFlag(fft.HalfSpectrum);
    std::complex<T> wiener;
    std::vector<std::complex<T>> freqVec;
    std::vector<std::complex<T>> responseFFT;
    fft.fwd(freqVec, inputWaveform[i]);
    fft.fwd(responseFFT, responseFunction);
    for (size_t j=0; j<responseFFT.size(); ++j) {
      wiener = std::conj(responseFFT[j]) / 
      ((std::complex<T>) std::pow(std::abs(responseFFT[j]), 2.0) + 
        noiseVar / (std::complex<T>) std::pow(std::abs(freqVec[j]), 2.0));
      freqVec[j] = freqVec[j] * wiener;
    }
    fft.inv(outputWaveform[i], freqVec);
    fft.ClearFlag(fft.HalfSpectrum);
  }
  return;
}

// Fourier Shrinkage

void sigproc_tools::Deconvolution::FourierShrinkage1D(
  std::vector<std::vector<float>>& outputWaveform,
  const std::vector<std::vector<float>>& inputWaveform,
  const std::vector<float>& responseFunction,
  const float regParam)
{
  size_t numChannels = inputWaveform.size();
  size_t nTicks = inputWaveform.at(0).size();
  outputWaveform.resize(numChannels);
  for (size_t i=0; i<numChannels; ++i) {
    outputWaveform[i].resize(nTicks);
  }
  for (size_t i=0; i<numChannels; ++i) {
    Eigen::FFT<float> fft;
    fft.SetFlag(fft.HalfSpectrum);
    std::complex<float> wienerReg;
    std::vector<std::complex<float>> freqVec;
    std::vector<std::complex<float>> responseFFT;
    fft.fwd(freqVec, inputWaveform[i]);
    fft.fwd(responseFFT, responseFunction);
    for (size_t j=0; j<responseFFT.size(); ++j) {
      wienerReg = std::conj(responseFFT[j]) / 
            ((std::complex<float>) std::pow(std::abs(responseFFT[j]), 2.0) + 
              (std::complex<float>) regParam);
      freqVec[j] = freqVec[j] * wienerReg;
    }
    fft.inv(outputWaveform[i], freqVec);
    fft.ClearFlag(fft.HalfSpectrum);
  }
  return;
}


#endif

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

// 2D PseudoWiener Filtering


#endif

#ifndef __SIGPROC_TOOLS_MORPH2D_CXX__
#define __SIGPROC_TOOLS_MORPH2D_CXX__

#include "Morph2D.h"


void icarus_signal_processing::Morph2D::getFilter2D(
  const std::vector<std::vector<short> >& waveform2D,
  const unsigned int structuringElementx,
  const unsigned int structuringElementy,
  std::vector<std::vector<short> >& dilation2D,
  std::vector<std::vector<short> >& erosion2D,
  std::vector<std::vector<short> >& average2D,
  std::vector<std::vector<short> >& gradient2D) const
{
  getFilter2D<short>(waveform2D, structuringElementx, structuringElementy,
    dilation2D, erosion2D, average2D, gradient2D);
  return;
}

void icarus_signal_processing::Morph2D::getFilter2D(
  const std::vector<std::vector<float> >& waveform2D,
  const unsigned int structuringElementx,
  const unsigned int structuringElementy,
  std::vector<std::vector<float> >& dilation2D,
  std::vector<std::vector<float> >& erosion2D,
  std::vector<std::vector<float> >& average2D,
  std::vector<std::vector<float> >& gradient2D) const
{
  getFilter2D<float>(waveform2D, structuringElementx, structuringElementy,
    dilation2D, erosion2D, average2D, gradient2D);
  return;
}

void icarus_signal_processing::Morph2D::getFilter2D(
  const std::vector<std::vector<double> >& waveform2D,
  const unsigned int structuringElementx,
  const unsigned int structuringElementy,
  std::vector<std::vector<double> >& dilation2D,
  std::vector<std::vector<double> >& erosion2D,
  std::vector<std::vector<double> >& average2D,
  std::vector<std::vector<double> >& gradient2D) const
{
  getFilter2D<double>(waveform2D, structuringElementx, structuringElementy,
    dilation2D, erosion2D, average2D, gradient2D);
  return;
}

template <typename T>
void icarus_signal_processing::Morph2D::getFilter2D(
  const std::vector<std::vector<T> >& waveform2D,
  const unsigned int structuringElementx,
  const unsigned int structuringElementy,
  std::vector<std::vector<T> >& dilation2D,
  std::vector<std::vector<T> >& erosion2D,
  std::vector<std::vector<T> >& average2D,
  std::vector<std::vector<T> >& gradient2D) const
{
  auto numChannels = waveform2D.size();
  auto nTicks = waveform2D.at(0).size();
  int xHalfWindowSize(structuringElementx / 2);
  int yHalfWindowSize(structuringElementy / 2);

  dilation2D.resize(waveform2D.size());
  erosion2D.resize(waveform2D.size());
  average2D.resize(waveform2D.size());
  gradient2D.resize(waveform2D.size());

  for (size_t i=0; i<waveform2D.size(); ++i) {
    dilation2D[i].resize(waveform2D.at(0).size());
    erosion2D[i].resize(waveform2D.at(0).size());
    average2D[i].resize(waveform2D.at(0).size());
    gradient2D[i].resize(waveform2D.at(0).size());
  }
  float dilation;
  float erosion;
  float gradient;
  float average;
  for (size_t i=0; i<numChannels; ++i) {
    for (size_t j=0; j<nTicks; ++j) {
      // For each center pixel, do 2D morphological filtering.
      int lbx = i - (int) xHalfWindowSize;
      int ubx = i + (int) xHalfWindowSize;
      int lby = j - (int) yHalfWindowSize;
      int uby = j + (int) yHalfWindowSize;
      size_t lowerBoundx = std::max(lbx, 0);
      size_t upperBoundx = std::min(ubx, (int) numChannels);
      size_t lowerBoundy = std::max(lby, 0);
      size_t upperBoundy = std::min(uby, (int) nTicks);
      std::vector<T> v;
      v.reserve(structuringElementx * structuringElementy);
      for (size_t ix=lowerBoundx; ix<upperBoundx; ++ix) {
        for (size_t iy=lowerBoundy; iy<upperBoundy; ++iy) {
          v.push_back(waveform2D[ix][iy]);
        }
      }
      dilation = *std::max_element(v.begin(), v.end());
      erosion = *std::min_element(v.begin(), v.end());
      average = 0.5 * (dilation + erosion);
      gradient = dilation - erosion;
      dilation2D[i][j] = dilation;
      erosion2D[i][j] = erosion;
      average2D[i][j] = average;
      gradient2D[i][j] = gradient;
    }
  }
  return;
}


void icarus_signal_processing::Morph2D::getDilation(
  const std::vector<std::vector<short> >& waveform2D,
  const unsigned int structuringElementx,
  const unsigned int structuringElementy,
  std::vector<std::vector<short> >& dilation2D) const
{
  getDilation<short>(waveform2D, structuringElementx, 
    structuringElementy, dilation2D);
  return;
}

void icarus_signal_processing::Morph2D::getDilation(
  const std::vector<std::vector<float> >& waveform2D,
  const unsigned int structuringElementx,
  const unsigned int structuringElementy,
  std::vector<std::vector<float> >& dilation2D) const
{
  getDilation<float>(waveform2D, structuringElementx, 
    structuringElementy, dilation2D);
  return;
}

void icarus_signal_processing::Morph2D::getDilation(
  const std::vector<std::vector<double> >& waveform2D,
  const unsigned int structuringElementx,
  const unsigned int structuringElementy,
  std::vector<std::vector<double> >& dilation2D) const
{
  getDilation<double>(waveform2D, structuringElementx, 
    structuringElementy, dilation2D);
  return;
}

template <typename T>
void icarus_signal_processing::Morph2D::getDilation(
  const std::vector<std::vector<T> >& waveform2D,
  const unsigned int structuringElementx,
  const unsigned int structuringElementy,
  std::vector<std::vector<T> >& dilation2D) const
{
  auto numChannels = waveform2D.size();
  auto nTicks = waveform2D.at(0).size();
  int xHalfWindowSize(structuringElementx / 2);
  int yHalfWindowSize(structuringElementy / 2);

  dilation2D.resize(waveform2D.size());

  for (size_t i=0; i<waveform2D.size(); ++i) {
    dilation2D[i].resize(waveform2D.at(0).size());
  }
  float dilation;
  for (size_t i=0; i<numChannels; ++i) {
    for (size_t j=0; j<nTicks; ++j) {
      // For each center pixel, do 2D morphological filtering.
      int lbx = i - (int) xHalfWindowSize;
      int ubx = i + (int) xHalfWindowSize;
      int lby = j - (int) yHalfWindowSize;
      int uby = j + (int) yHalfWindowSize;
      size_t lowerBoundx = std::max(lbx, 0);
      size_t upperBoundx = std::min(ubx, (int) numChannels);
      size_t lowerBoundy = std::max(lby, 0);
      size_t upperBoundy = std::min(uby, (int) nTicks);
      std::vector<T> v;
      v.reserve(structuringElementx * structuringElementy);
      for (size_t ix=lowerBoundx; ix<upperBoundx; ++ix) {
        for (size_t iy=lowerBoundy; iy<upperBoundy; ++iy) {
          v.push_back(waveform2D[ix][iy]);
        }
      }
      dilation = *std::max_element(v.begin(), v.end());
      dilation2D[i][j] = dilation;
    }
  }
  return;
}


void icarus_signal_processing::Morph2D::getErosion(
  const std::vector<std::vector<short> >& waveform2D,
  const unsigned int structuringElementx,
  const unsigned int structuringElementy,
  std::vector<std::vector<short> >& erosion2D) const
{
  getErosion<short>(waveform2D, structuringElementx, 
    structuringElementy, erosion2D);
  return;
}

void icarus_signal_processing::Morph2D::getErosion(
  const std::vector<std::vector<float> >& waveform2D,
  const unsigned int structuringElementx,
  const unsigned int structuringElementy,
  std::vector<std::vector<float> >& erosion2D) const
{
  getErosion<float>(waveform2D, structuringElementx, 
    structuringElementy, erosion2D);
  return;
}

void icarus_signal_processing::Morph2D::getErosion(
  const std::vector<std::vector<double> >& waveform2D,
  const unsigned int structuringElementx,
  const unsigned int structuringElementy,
  std::vector<std::vector<double> >& erosion2D) const
{
  getErosion<double>(waveform2D, structuringElementx, 
    structuringElementy, erosion2D);
  return;
}

template <typename T>
void icarus_signal_processing::Morph2D::getErosion(
  const std::vector<std::vector<T> >& waveform2D,
  const unsigned int structuringElementx,
  const unsigned int structuringElementy,
  std::vector<std::vector<T> >& erosion2D) const
{
  auto numChannels = waveform2D.size();
  auto nTicks = waveform2D.at(0).size();
  int xHalfWindowSize(structuringElementx / 2);
  int yHalfWindowSize(structuringElementy / 2);

  erosion2D.resize(waveform2D.size());

  for (size_t i=0; i<waveform2D.size(); ++i) {
    erosion2D[i].resize(waveform2D.at(0).size());
  }
  float erosion;
  for (size_t i=0; i<numChannels; ++i) {
    for (size_t j=0; j<nTicks; ++j) {
      // For each center pixel, do 2D morphological filtering.
      int lbx = i - (int) xHalfWindowSize;
      int ubx = i + (int) xHalfWindowSize;
      int lby = j - (int) yHalfWindowSize;
      int uby = j + (int) yHalfWindowSize;
      size_t lowerBoundx = std::max(lbx, 0);
      size_t upperBoundx = std::min(ubx, (int) numChannels);
      size_t lowerBoundy = std::max(lby, 0);
      size_t upperBoundy = std::min(uby, (int) nTicks);
      std::vector<T> v;
      v.reserve(structuringElementx * structuringElementy);
      for (size_t ix=lowerBoundx; ix<upperBoundx; ++ix) {
        for (size_t iy=lowerBoundy; iy<upperBoundy; ++iy) {
          v.push_back(waveform2D[ix][iy]);
        }
      }
      erosion = *std::min_element(v.begin(), v.end());
      erosion2D[i][j] = erosion;
    }
  }
  return;
}


void icarus_signal_processing::Morph2D::getGradient(
  const std::vector<std::vector<short> >& waveform2D,
  const unsigned int structuringElementx,
  const unsigned int structuringElementy,
  std::vector<std::vector<short> >& gradient2D) const
{
  getGradient<short>(waveform2D, structuringElementx, 
    structuringElementy, gradient2D);
  return;
}

void icarus_signal_processing::Morph2D::getGradient(
  const std::vector<std::vector<float> >& waveform2D,
  const unsigned int structuringElementx,
  const unsigned int structuringElementy,
  std::vector<std::vector<float> >& gradient2D) const
{
  getGradient<float>(waveform2D, structuringElementx, 
    structuringElementy, gradient2D);
  return;
}

void icarus_signal_processing::Morph2D::getGradient(
  const std::vector<std::vector<double> >& waveform2D,
  const unsigned int structuringElementx,
  const unsigned int structuringElementy,
  std::vector<std::vector<double> >& gradient2D) const
{
  getGradient<double>(waveform2D, structuringElementx, 
    structuringElementy, gradient2D);
  return;
}

template <typename T>
void icarus_signal_processing::Morph2D::getGradient(
  const std::vector<std::vector<T> >& waveform2D,
  const unsigned int structuringElementx,
  const unsigned int structuringElementy,
  std::vector<std::vector<T> >& gradient2D) const
{
  auto numChannels = waveform2D.size();
  auto nTicks = waveform2D.at(0).size();
  int xHalfWindowSize(structuringElementx / 2);
  int yHalfWindowSize(structuringElementy / 2);

  gradient2D.resize(waveform2D.size());

  for (size_t i=0; i<waveform2D.size(); ++i) {
    gradient2D[i].resize(waveform2D.at(0).size());
  }
  float dilation;
  float erosion;
  float gradient;
  for (size_t i=0; i<numChannels; ++i) {
    for (size_t j=0; j<nTicks; ++j) {
      // For each center pixel, do 2D morphological filtering.
      int lbx = i - (int) xHalfWindowSize;
      int ubx = i + (int) xHalfWindowSize;
      int lby = j - (int) yHalfWindowSize;
      int uby = j + (int) yHalfWindowSize;
      size_t lowerBoundx = std::max(lbx, 0);
      size_t upperBoundx = std::min(ubx, (int) numChannels);
      size_t lowerBoundy = std::max(lby, 0);
      size_t upperBoundy = std::min(uby, (int) nTicks);
      std::vector<T> v;
      v.reserve(structuringElementx * structuringElementy);
      for (size_t ix=lowerBoundx; ix<upperBoundx; ++ix) {
        for (size_t iy=lowerBoundy; iy<upperBoundy; ++iy) {
          v.push_back(waveform2D[ix][iy]);
        }
      }
      dilation = *std::max_element(v.begin(), v.end());
      erosion = *std::min_element(v.begin(), v.end());
      gradient = dilation - erosion;
      gradient2D[i][j] = gradient;
    }
  }
  return;
}


void icarus_signal_processing::Morph2D::getMedian(
  const std::vector<std::vector<short> >& waveform2D,
  const unsigned int structuringElementx,
  const unsigned int structuringElementy,
  std::vector<std::vector<short> >& median2D) const
{
  getMedian<short>(waveform2D, structuringElementx, 
    structuringElementy, median2D);
  return;
}

void icarus_signal_processing::Morph2D::getMedian(
  const std::vector<std::vector<float> >& waveform2D,
  const unsigned int structuringElementx,
  const unsigned int structuringElementy,
  std::vector<std::vector<float> >& median2D) const
{
  getMedian<float>(waveform2D, structuringElementx, 
    structuringElementy, median2D);
  return;
}

void icarus_signal_processing::Morph2D::getMedian(
  const std::vector<std::vector<double> >& waveform2D,
  const unsigned int structuringElementx,
  const unsigned int structuringElementy,
  std::vector<std::vector<double> >& median2D) const
{
  getMedian<double>(waveform2D, structuringElementx, 
    structuringElementy, median2D);
  return;
}

template <typename T>
void icarus_signal_processing::Morph2D::getMedian(
  const std::vector<std::vector<T> >& waveform2D,
  const unsigned int structuringElementx,
  const unsigned int structuringElementy,
  std::vector<std::vector<T> >& median2D) const
{
  auto numChannels = waveform2D.size();
  auto nTicks = waveform2D.at(0).size();
  int xHalfWindowSize(structuringElementx / 2);
  int yHalfWindowSize(structuringElementy / 2);

  median2D.resize(waveform2D.size());

  for (size_t i=0; i<waveform2D.size(); ++i) {
    median2D[i].resize(waveform2D.at(0).size());
  }
  float median;
  for (size_t i=0; i<numChannels; ++i) {
    for (size_t j=0; j<nTicks; ++j) {
      // For each center pixel, do 2D morphological filtering.
      int lbx = i - (int) xHalfWindowSize;
      int ubx = i + (int) xHalfWindowSize;
      int lby = j - (int) yHalfWindowSize;
      int uby = j + (int) yHalfWindowSize;
      size_t lowerBoundx = std::max(lbx, 0);
      size_t upperBoundx = std::min(ubx, (int) numChannels);
      size_t lowerBoundy = std::max(lby, 0);
      size_t upperBoundy = std::min(uby, (int) nTicks);
      std::vector<T> v;
      v.reserve(structuringElementx * structuringElementy);
      for (size_t ix=lowerBoundx; ix<upperBoundx; ++ix) {
        for (size_t iy=lowerBoundy; iy<upperBoundy; ++iy) {
          v.push_back(waveform2D[ix][iy]);
        }
      }
      if (v.size() % 2 == 0) {
        const auto m1 = v.begin() + v.size() / 2 - 1;
        const auto m2 = v.begin() + v.size() / 2;
        std::nth_element(v.begin(), m1, v.end());
        std::nth_element(v.begin(), m2, v.end());
        median = (*m1 + *m2) / 2.0;
      } else {
        median = v[v.size() / 2];
      }
      median2D[i][j] = median;
    }
  }
  return;
}


void icarus_signal_processing::Morph2D::getOpeningAndClosing(
  const std::vector<std::vector<short> >& waveform2D,
  const unsigned int structuringElementx,
  const unsigned int structuringElementy,
  std::vector<std::vector<short> >& opening2D,
  std::vector<std::vector<short> >& closing2D) const
{
  getOpeningAndClosing<short>(
    waveform2D, 
    structuringElementx, 
    structuringElementy, 
    opening2D,
    closing2D);
  return;
}

void icarus_signal_processing::Morph2D::getOpeningAndClosing(
  const std::vector<std::vector<float> >& waveform2D,
  const unsigned int structuringElementx,
  const unsigned int structuringElementy,
  std::vector<std::vector<float> >& opening2D,
  std::vector<std::vector<float> >& closing2D) const
{
  getOpeningAndClosing<float>(
    waveform2D, 
    structuringElementx, 
    structuringElementy, 
    opening2D,
    closing2D);
  return;
}

void icarus_signal_processing::Morph2D::getOpeningAndClosing(
  const std::vector<std::vector<double> >& waveform2D,
  const unsigned int structuringElementx,
  const unsigned int structuringElementy,
  std::vector<std::vector<double> >& opening2D,
  std::vector<std::vector<double> >& closing2D) const
{
  getOpeningAndClosing<double>(
    waveform2D, 
    structuringElementx, 
    structuringElementy, 
    opening2D,
    closing2D);
  return;
}


template <typename T>
void icarus_signal_processing::Morph2D::getOpeningAndClosing(
  const std::vector<std::vector<T> >& waveform2D,
  const unsigned int structuringElementx,
  const unsigned int structuringElementy,
  std::vector<std::vector<T> >& opening2D,
  std::vector<std::vector<T> >& closing2D) const
{
  auto numChannels = waveform2D.size();
  auto nTicks = waveform2D.at(0).size();
  int xHalfWindowSize(structuringElementx / 2);
  int yHalfWindowSize(structuringElementy / 2);

  opening2D.resize(waveform2D.size());
  closing2D.resize(waveform2D.size());

  std::vector<std::vector<T> > dilation2D;
  std::vector<std::vector<T> > erosion2D;
  std::vector<std::vector<T> > average2D;
  std::vector<std::vector<T> > gradient2D;

  getFilter2D(waveform2D, structuringElementx, structuringElementy,
              dilation2D, erosion2D, average2D, gradient2D);

  for (size_t i=0; i<waveform2D.size(); ++i) {
    opening2D[i].resize(waveform2D.at(0).size());
    closing2D[i].resize(waveform2D.at(0).size());
  }
  float opening;
  float closing;
  for (size_t i=0; i<numChannels; ++i) {
    for (size_t j=0; j<nTicks; ++j) {
      // For each center pixel, do 2D morphological filtering.
      int lbx = i - (int) xHalfWindowSize;
      int ubx = i + (int) xHalfWindowSize;
      int lby = j - (int) yHalfWindowSize;
      int uby = j + (int) yHalfWindowSize;
      size_t lowerBoundx = std::max(lbx, 0);
      size_t upperBoundx = std::min(ubx, (int) numChannels);
      size_t lowerBoundy = std::max(lby, 0);
      size_t upperBoundy = std::min(uby, (int) nTicks);
      std::vector<T> v1;
      std::vector<T> v2;
      v1.reserve(structuringElementx * structuringElementy);
      v2.reserve(structuringElementx * structuringElementy);
      for (size_t ix=lowerBoundx; ix<upperBoundx; ++ix) {
        for (size_t iy=lowerBoundy; iy<upperBoundy; ++iy) {
          v1.push_back(dilation2D[ix][iy]);
          v2.push_back(erosion2D[ix][iy]);
        }
      }
      opening = *std::max_element(v2.begin(), v2.end());
      closing = *std::min_element(v1.begin(), v1.end());
      opening2D[i][j] = opening;
      closing2D[i][j] = closing;
    }
  }
  return;
}

#endif

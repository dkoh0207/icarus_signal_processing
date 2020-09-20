#ifndef __SIGPROC_TOOLS_MISCUTILS_CXX__
#define __SIGPROC_TOOLS_MISCUTILS_CXX__

#include "MiscUtils.h"

short sigproc_tools::MiscUtils::computeMedian(const std::vector<short>& vec) {
  short median = computeMedian<short>(vec);
  return median;
}

float sigproc_tools::MiscUtils::computeMedian(const std::vector<float>& vec) {
  float median = computeMedian<float>(vec);
  return median;
}

double sigproc_tools::MiscUtils::computeMedian(const std::vector<double>& vec) {
  double median = computeMedian<double>(vec);
  return median;
}

template <typename T>
T sigproc_tools::MiscUtils::computeMedian(const std::vector<T>& vec)
{
  T median;
  std::vector<T> localVec = vec;
  if (localVec.size() % 2 == 0) {
    const auto m1 = localVec.begin() + localVec.size() / 2 - 1;
    const auto m2 = localVec.begin() + localVec.size() / 2;
    std::nth_element(localVec.begin(), m1, localVec.end());
    const auto e1 = *m1;
    std::nth_element(localVec.begin(), m2, localVec.end());
    const auto e2 = *m2;
    median = (e1 + e2) / 2.0;
  } else {
    const auto m = localVec.begin() + localVec.size() / 2;
    std::nth_element(localVec.begin(), m, localVec.end());
    median = *m;
  }
  return median;
}

float sigproc_tools::MiscUtils::estimateNoiseVariance(
  const std::vector<std::vector<float>>& waveLessCoherent,
  const std::vector<std::vector<bool>>& selectVals)
{
  size_t numChannels = waveLessCoherent.size();
  size_t nTicks = waveLessCoherent.at(0).size();

  float var = 0.0;
  float mean = 0.0;
  int count = 0;

  for (size_t i=0; i<numChannels; ++i) {
    for (size_t j=0; j<nTicks; ++j) {
      if (!selectVals[i][j]) {
        var += std::pow(waveLessCoherent[i][j], 2.0);
        mean += waveLessCoherent[i][j];
        count += 1;
      }
    }
  }
  mean = mean / ((float) count);
  var = var / ((float) count) - std::pow(mean, 2.0);
  return var;
}

float sigproc_tools::MiscUtils::estimateMAD(
  const std::vector<float>& wf)
{
  size_t numChannels = wf.size();

  std::vector<float> wksp;
  wksp.reserve(numChannels);

  for (size_t i=0; i<numChannels; ++i) {
    wksp.push_back(std::abs(wf[i]));
  }

  float median;
  median = computeMedian(wksp);
  return median / 0.6745;
}



#endif

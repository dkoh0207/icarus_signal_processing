#ifndef __icarus_signal_processing_MISCUTILS_CXX__
#define __icarus_signal_processing_MISCUTILS_CXX__

#include "MiscUtils.h"

short icarus_signal_processing::MiscUtils::computeMedian(const std::vector<short>& vec) {
  short median = computeMedian<short>(vec);
  return median;
}

float icarus_signal_processing::MiscUtils::computeMedian(const std::vector<float>& vec) {
  float median = computeMedian<float>(vec);
  return median;
}

double icarus_signal_processing::MiscUtils::computeMedian(const std::vector<double>& vec) {
  double median = computeMedian<double>(vec);
  return median;
}

template <typename T>
T icarus_signal_processing::MiscUtils::computeMedian(const std::vector<T>& vec)
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

short icarus_signal_processing::MiscUtils::computeMaximum(const std::vector<std::vector<short>>& input2D) {
    short res = computeMaximum<short>(input2D);
    return res;
}

float icarus_signal_processing::MiscUtils::computeMaximum(const std::vector<std::vector<float>>& input2D) {
    float res = computeMaximum<float>(input2D);
    return res;
}

double icarus_signal_processing::MiscUtils::computeMaximum(const std::vector<std::vector<double>>& input2D) {
    double res = computeMaximum<double>(input2D);
    return res;
}

template <typename T>
T icarus_signal_processing::MiscUtils::computeMaximum(const std::vector<std::vector<T>>& input2D) 
{

    T res = std::numeric_limits<T>::min();

    for (size_t i=0; i<input2D.size(); ++i) {
        for (size_t j=0; j<input2D.at(0).size(); ++j) {
            if (input2D[i][j] > res) res = input2D[i][j];
        }
    }

    return res;
}

short icarus_signal_processing::MiscUtils::computeMinimum(const std::vector<std::vector<short>>& input2D) {
    short res = computeMinimum<short>(input2D);
    return res;
}

float icarus_signal_processing::MiscUtils::computeMinimum(const std::vector<std::vector<float>>& input2D) {
    float res = computeMinimum<float>(input2D);
    return res;
}

double icarus_signal_processing::MiscUtils::computeMinimum(const std::vector<std::vector<double>>& input2D) {
    double res = computeMinimum<double>(input2D);
    return res;
}

template <typename T>
T icarus_signal_processing::MiscUtils::computeMinimum(const std::vector<std::vector<T>>& input2D) 
{
    T res = std::numeric_limits<T>::max();

    for (size_t i=0; i<input2D.size(); ++i) {
        for (size_t j=0; j<input2D.at(0).size(); ++j) {
            if (input2D[i][j] < res) res = input2D[i][j];
        }
    }

  return res;
}

float icarus_signal_processing::MiscUtils::estimateNoiseVariance(
  const std::vector<std::vector<float>>& waveLessCoherent,
  const std::vector<std::vector<bool>>& selectVals)
{
    size_t numChannels = waveLessCoherent.size();
    size_t nTicks      = waveLessCoherent.at(0).size();

    float var   = 0.0;
    float mean  = 0.0;
    int   count = 0;

    for (size_t i=0; i<numChannels; ++i) {
        for (size_t j=0; j<nTicks; ++j) {
            if (!selectVals[i][j]) {
                var   += std::pow(waveLessCoherent[i][j], 2.0);
                mean  += waveLessCoherent[i][j];
                count += 1;
            }
        }
    }

    mean = mean / ((float) count);
    var  = var  / ((float) count) - std::pow(mean, 2.0);
    
    return var;
}

float icarus_signal_processing::MiscUtils::estimateMAD(
  const std::vector<float>& wf)
{
    size_t numChannels = wf.size();

    std::vector<float> wksp;
    wksp.reserve(numChannels);

    for (size_t i=0; i<numChannels; ++i) wksp.push_back(std::abs(wf[i]));

    float median = computeMedian(wksp);
  
    return median / 0.6745;
}


unsigned icarus_signal_processing::MiscUtils::nChoosek( unsigned n, unsigned k ) const
{
    if (k > n) return 0;
    if (k * 2 > n) k = n-k;
    if (k == 0) return 1;

    unsigned result = n;
    for( unsigned i = 2; i <= k; ++i ) {
        result *= (n-i+1);
        result /= i;
    }
    return result;
}


void icarus_signal_processing::MiscUtils::drawIndextoImage(
  const std::vector<int>& rows,
  const std::vector<int>& cols,
  std::vector<std::vector<bool>>& output2D) const
{
    drawIndextoImage<int>(rows, cols, output2D);
    return;
}

void icarus_signal_processing::MiscUtils::drawIndextoImage(
  const std::vector<size_t>& rows,
  const std::vector<size_t>& cols,
  std::vector<std::vector<bool>>& output2D) const
{
    drawIndextoImage<size_t>(rows, cols, output2D);
    return;
}

void icarus_signal_processing::MiscUtils::drawIndextoImage(
  const std::vector<unsigned int>& rows,
  const std::vector<unsigned int>& cols,
  std::vector<std::vector<bool>>& output2D) const
{
    drawIndextoImage<unsigned int>(rows, cols, output2D);
    return;
}

void icarus_signal_processing::MiscUtils::drawIndextoImage(
  const std::vector<short>& rows,
  const std::vector<short>& cols,
  std::vector<std::vector<bool>>& output2D) const
{
    drawIndextoImage<short>(rows, cols, output2D);
    return;
}

template <typename T>
void icarus_signal_processing::MiscUtils::drawIndextoImage(
  const std::vector<T>& rows,
  const std::vector<T>& cols,
  std::vector<std::vector<bool>>& output2D) const
{
    if (rows.size() != cols.size())
      throw std::runtime_error("drawIndexToImage: rows and columns don't match");

    size_t N = rows.size();

    for (size_t i=0; i<N; ++i) {
        T ix = rows[i];
        T iy = cols[i];
        output2D[ix][iy] = true;
    }

    return;

}


#endif

#ifndef __SIGPROC_TOOLS_DENOISING_CXX__
#define __SIGPROC_TOOLS_DENOISING_CXX__

#include "Denoising.h"

#include <chrono>


void icarus_signal_processing::Denoising::getSelectVals(
  ArrayShort::const_iterator  waveforms,
  ArrayShort::const_iterator  morphedWaveforms,
  ArrayBool::iterator         selectVals,
  ArrayBool::iterator         roi,
  VectorFloat::const_iterator thresholdVec,
  const unsigned int          numChannels,
  const unsigned int          window)
{
  getSelectVals<short>(waveforms, morphedWaveforms, selectVals,
    roi, thresholdVec, numChannels, window);
}

void icarus_signal_processing::Denoising::getSelectVals(
  ArrayFloat::const_iterator waveforms,
  ArrayFloat::const_iterator morphedWaveforms,
  ArrayBool::iterator         selectVals,
  ArrayBool::iterator         roi,
  VectorFloat::const_iterator thresholdVec,
  const unsigned int          numChannels,
  const unsigned int          window)
{
  getSelectVals<float>(waveforms, morphedWaveforms, selectVals,
    roi, thresholdVec, numChannels, window);
}

void icarus_signal_processing::Denoising::getSelectVals(
  ArrayDouble::const_iterator  waveforms,
  ArrayDouble::const_iterator  morphedWaveforms,
  ArrayBool::iterator          selectVals,
  ArrayBool::iterator          roi,
  VectorFloat::const_iterator  thresholdVec,
  const unsigned int           numChannels,
  const unsigned int           window)
{
  getSelectVals<double>(waveforms, morphedWaveforms, selectVals,
    roi, thresholdVec, numChannels, window);
}

template <typename T>
void icarus_signal_processing::Denoising::getSelectVals(
  typename std::vector<std::vector<T>>::const_iterator waveformsItr,
  typename std::vector<std::vector<T>>::const_iterator morphedWaveformsItr,
  ArrayBool::iterator                                  selectValsItr,
  ArrayBool::iterator                                  roiItr,
  VectorFloat::const_iterator                          thresholdItr,
  const unsigned int                                   numChannels,
  const unsigned int                                   window)
{
    auto nTicks = waveformsItr[0].size();

    // Set a protection width
    int halfWidth = std::max(int(window/8),1);

    for (size_t i=0; i<numChannels; ++i) 
    {
        std::vector<T> localVec = morphedWaveformsItr[i];

        T median = getMedian(localVec, localVec.size());

        std::vector<T> baseVec;
        baseVec.resize(localVec.size());

        for (size_t j=0; j<baseVec.size(); ++j) baseVec[j] = morphedWaveformsItr[i][j] - median;

        float rms;
        rms = std::sqrt(std::inner_product(baseVec.begin(), baseVec.end(), baseVec.begin(), 0.) / float(baseVec.size()));
        float threshold;
        threshold = (*(thresholdItr + i)) * rms;

        for (size_t j=0; j<nTicks; ++j) {
            if (std::fabs(morphedWaveformsItr[i][j]) > threshold) {
                // Check Bounds
                selectValsItr[i][j] = true;
                //int lb = j - (int) window/4;
                //int ub = j + (int) window/4 + 1;
                int lb = j - halfWidth;
                int ub = j + halfWidth + 1;
                size_t lowerBound = std::max(lb, 0);
                size_t upperBound = std::min(ub, (int) nTicks);
                for (size_t k=lowerBound; k<upperBound; ++k) {
                    roiItr[i][k] = true;
                }
            } else {
                selectValsItr[i][j] = false;
            }
        }
    }

    return;
}


void icarus_signal_processing::Denoising::removeCoherentNoise1D(
  ArrayShort::iterator              waveLessCoherent,
  ArrayShort::const_iterator        filteredWaveforms,
  ArrayShort::iterator              morphedWaveforms,
  ArrayShort::iterator              intrinsicRMS,
  ArrayBool::iterator               selectVals,
  ArrayBool::iterator               roi,
  ArrayShort::iterator              correctedMedians,
  FilterFunctionVec::const_iterator filterFunctions,
  VectorFloat::const_iterator       thresholdVec,
  const unsigned int                numChannels,
  const unsigned int                grouping,
  const unsigned int                window)
{
  removeCoherentNoise1D<short>(
    waveLessCoherent, filteredWaveforms, morphedWaveforms, 
    intrinsicRMS, selectVals, roi, correctedMedians, filterFunctions, thresholdVec,
    numChannels, grouping, window);
  return;
}

void icarus_signal_processing::Denoising::removeCoherentNoise1D(
  ArrayFloat::iterator              waveLessCoherent,
  ArrayFloat::const_iterator        filteredWaveforms,
  ArrayFloat::iterator              morphedWaveforms,
  ArrayFloat::iterator              intrinsicRMS,
  ArrayBool::iterator               selectVals,
  ArrayBool::iterator               roi,
  ArrayFloat::iterator              correctedMedians,
  FilterFunctionVec::const_iterator filterFunctions,
  VectorFloat::const_iterator       thresholdVec,
  const unsigned int                numChannels,
  const unsigned int                grouping,
  const unsigned int                window)
{
  removeCoherentNoise1D<float>(
    waveLessCoherent, filteredWaveforms, morphedWaveforms, 
    intrinsicRMS, selectVals, roi, correctedMedians, filterFunctions, thresholdVec,
    numChannels, grouping, window);
  return;
}

void icarus_signal_processing::Denoising::removeCoherentNoise1D(
  ArrayDouble::iterator             waveLessCoherent,
  ArrayDouble::const_iterator       filteredWaveforms,
  ArrayDouble::iterator             morphedWaveforms,
  ArrayDouble::iterator             intrinsicRMS,
  ArrayBool::iterator               selectVals,
  ArrayBool::iterator               roi,
  ArrayDouble::iterator             correctedMedians,
  FilterFunctionVec::const_iterator filterFunctions,
  VectorFloat::const_iterator       thresholdVec,
  const unsigned int                numChannels,
  const unsigned int                grouping,
  const unsigned int                window)
{
  removeCoherentNoise1D<double>(
    waveLessCoherent, filteredWaveforms, morphedWaveforms, 
    intrinsicRMS, selectVals, roi, correctedMedians, filterFunctions, thresholdVec, 
    numChannels, grouping, window);
  return;
}

template <typename T>
void icarus_signal_processing::Denoising::removeCoherentNoise1D(typename std::vector<std::vector<T>>::iterator       waveLessCoherentItr,
                                                                typename std::vector<std::vector<T>>::const_iterator filteredWaveformsItr,
                                                                typename std::vector<std::vector<T>>::iterator       morphedWaveformsItr,
                                                                typename std::vector<std::vector<T>>::iterator       intrinsicRMSItr,
                                                                ArrayBool::iterator                                  selectValsItr,
                                                                ArrayBool::iterator                                  roiItr,
                                                                typename std::vector<std::vector<T>>::iterator       correctedMediansItr,
                                                                FilterFunctionVec::const_iterator                    filterFunctions,
                                                                VectorFloat::const_iterator                          thresholdItr,
                                                                const unsigned int                                   numChannels,
                                                                const unsigned int                                   grouping,
                                                                const unsigned int                                   window)
{
    auto nTicks  = filteredWaveformsItr->size();
    auto nGroups = numChannels / grouping;

    std::chrono::high_resolution_clock::time_point funcStartTime = std::chrono::high_resolution_clock::now();

    std::chrono::high_resolution_clock::time_point morphStart = funcStartTime;

    for(size_t funcIdx = 0; funcIdx < numChannels; funcIdx++)
    {
        const icarus_signal_processing::IMorphologicalFunctions1D* func = (*(filterFunctions + funcIdx)).get();

        if (!func) 
        {
            std::cout << "Found null function for funcIdx " << funcIdx << " of " << numChannels << std::endl;
            continue;
        }

        (*func)(*(filteredWaveformsItr + funcIdx), *(morphedWaveformsItr + funcIdx));
    }

    std::chrono::high_resolution_clock::time_point morphStop  = std::chrono::high_resolution_clock::now();
    std::chrono::high_resolution_clock::time_point selStart = morphStop;

    getSelectVals(filteredWaveformsItr, morphedWaveformsItr, selectValsItr, roiItr, thresholdItr, numChannels, window);

    std::chrono::high_resolution_clock::time_point selStop  = std::chrono::high_resolution_clock::now();
    std::chrono::high_resolution_clock::time_point noiseStart = selStop;

    std::vector<T> v(grouping);

    for (size_t i=0; i<nTicks; ++i) 
    {
        for (size_t j=0; j<nGroups; ++j) 
        {
            size_t group_start = j * grouping;
            size_t group_end = (j+1) * grouping;
            // Compute median.
            size_t idxV(0);

            for (size_t c=group_start; c<group_end; ++c) 
            {
                if (!selectValsItr[c][i]) v[idxV++] = filteredWaveformsItr[c][i];
            }

            T median = getMedian(v,idxV);

            correctedMediansItr[j][i] = median;
            for (auto k=group_start; k<group_end; ++k) waveLessCoherentItr[k][i] = filteredWaveformsItr[k][i] - median;
        }
    }

    std::chrono::high_resolution_clock::time_point noiseStop = std::chrono::high_resolution_clock::now();

    T rms = (T) 0;
    for (size_t i=0; i<nGroups; ++i) 
    {
        for (size_t j=0; j<nTicks; ++j) 
        {
            size_t idxV(0);
            for (size_t k=i*grouping; k<(i+1)*grouping; ++k) v[idxV++] = waveLessCoherentItr[k][j];
            rms = std::sqrt(std::inner_product(v.begin(), v.begin()+idxV, v.begin(), 0.) / T(v.size()));
            intrinsicRMSItr[i][j] = (T) rms;
        }
    }

    std::chrono::high_resolution_clock::time_point funcStopTime = std::chrono::high_resolution_clock::now();
  
    std::chrono::duration<double> funcTime   = std::chrono::duration_cast<std::chrono::duration<double>>(funcStopTime - funcStartTime);
    std::chrono::duration<double> morphTime  = std::chrono::duration_cast<std::chrono::duration<double>>(morphStop - morphStart);
    std::chrono::duration<double> selTime    = std::chrono::duration_cast<std::chrono::duration<double>>(selStop - selStart);
    std::chrono::duration<double> noiseTime  = std::chrono::duration_cast<std::chrono::duration<double>>(noiseStop - noiseStart);
  
    std::cout << "*** Denoising ***  - # channels: " << numChannels << ", ticks: " << nTicks << ", groups: " << nGroups << std::endl;
    std::cout << "                   - morph: " << morphTime.count() << ", sel: " << selTime.count() << ", noise: " << noiseTime.count() << ", total: " << funcTime.count() << std::endl;
  
    return;
}

template <typename T> T icarus_signal_processing::Denoising::getMedian(typename std::vector<T>& vals, const unsigned int nVals) const
{
    T median = T(0);

    if (nVals > 2) 
    {
        if (nVals % 2 == 0) 
        {
            const auto m1 = vals.begin() + nVals / 2 - 1;
            const auto m2 = vals.begin() + nVals / 2;
            std::nth_element(vals.begin(), m1, vals.begin() + nVals);
            const auto e1 = *m1;
            std::nth_element(vals.begin(), m2, vals.begin() + nVals);
            const auto e2 = *m2;
            median = (e1 + e2) / 2.0;
        } 
        else 
        {
            const auto m = vals.begin() + nVals / 2;
            std::nth_element(vals.begin(), m, vals.begin() + nVals);
            median = *m;
        }
    }

    return median;
}

template <typename T> T icarus_signal_processing::Denoising::getMostProbable(typename std::vector<T>& vals, const unsigned int nVals)
{
    T mostProbable = T(0);

    // Do a simple average if only a few bins
    if (nVals < 5)
    {
        mostProbable = std::accumulate(vals.begin(),vals.begin()+nVals,0.) / T(nVals);
    }
    // Otherwise try to form value around MP
    else
    {
        auto minVItr = std::min_element(vals.begin(),vals.begin()+nVals);
        auto maxVItr = std::max_element(vals.begin(),vals.begin()+nVals);

        T minValue = *minVItr;
        T maxValue = *maxVItr;

        size_t numBins = size_t(maxValue - minValue + 1);

        std::fill(fMPVec.begin(),fMPVec.begin() + numBins,0);
  
        for(typename std::vector<T>::iterator vItr = vals.begin(); vItr != vals.begin() + nVals; vItr++)
        {
            int binIdx = int(std::round((*vItr - minValue)));
  
            fMPVec[binIdx]++;
        }
  
        std::vector<int>::iterator maxItr = std::max_element(fMPVec.begin(),fMPVec.begin()+numBins);
  
        int count = *maxItr;
        T   mpVal = T(count) * T(std::distance(fMPVec.begin(),maxItr));

        auto cntItr = maxItr;

        while(cntItr != fMPVec.begin())
        {
            if (*(cntItr - 1) < 0.5 * *maxItr) break;

            count += *(cntItr - 1);
            mpVal += T(*(cntItr - 1)) * T(std::distance(fMPVec.begin(),cntItr - 1));

            cntItr--;
        }

        cntItr = maxItr;

        while(cntItr + 1 != fMPVec.end())
        {
            if (*(cntItr + 1) < 0.5 * *maxItr) break;

            count += *(cntItr + 1);
            mpVal += T(*(cntItr + 1)) * T(std::distance(fMPVec.begin(),cntItr + 1));

            cntItr++;
        }
  
        mostProbable = mpVal / T(count) + minValue;
    }

    return mostProbable;
}


void icarus_signal_processing::Denoising::removeCoherentNoise2D(
  ArrayShort& waveLessCoherent,
  const ArrayShort& filteredWaveforms,
  ArrayShort& morphedWaveforms,
  ArrayShort& intrinsicRMS,
  ArrayBool& selectVals,
  ArrayBool& roi,
  ArrayShort& correctedMedians,
  VectorFloat& thresholdVec,
  const char filterName,
  const unsigned int grouping,
  const unsigned int structuringElementx,
  const unsigned int structuringElementy,
  const unsigned int window)
{
  removeCoherentNoise2D<short>(
    waveLessCoherent, filteredWaveforms, morphedWaveforms, 
    intrinsicRMS, selectVals, roi, correctedMedians,
    thresholdVec, filterName, grouping, structuringElementx, structuringElementy, 
    window);
  return;
}

void icarus_signal_processing::Denoising::removeCoherentNoise2D(
  ArrayFloat& waveLessCoherent,
  const ArrayFloat& filteredWaveforms,
  ArrayFloat& morphedWaveforms,
  ArrayFloat& intrinsicRMS,
  ArrayBool& selectVals,
  ArrayBool& roi,
  ArrayFloat& correctedMedians,
  VectorFloat& thresholdVec,
  const char filterName,
  const unsigned int grouping,
  const unsigned int structuringElementx,
  const unsigned int structuringElementy,
  const unsigned int window)
{
  removeCoherentNoise2D<float>(
    waveLessCoherent, filteredWaveforms, morphedWaveforms, 
    intrinsicRMS, selectVals, roi, correctedMedians,
    thresholdVec, filterName, grouping, structuringElementx, structuringElementy, 
    window);
  return;
}

void icarus_signal_processing::Denoising::removeCoherentNoise2D(
  ArrayDouble& waveLessCoherent,
  const ArrayDouble& filteredWaveforms,
  ArrayDouble& morphedWaveforms,
  ArrayDouble& intrinsicRMS,
  ArrayBool& selectVals,
  ArrayBool& roi,
  ArrayDouble& correctedMedians,
  VectorFloat& thresholdVec,
  const char filterName,
  const unsigned int grouping,
  const unsigned int structuringElementx,
  const unsigned int structuringElementy,
  const unsigned int window)
{
  removeCoherentNoise2D<double>(
    waveLessCoherent, filteredWaveforms, morphedWaveforms, 
    intrinsicRMS, selectVals, roi, correctedMedians,
    thresholdVec, filterName, grouping, structuringElementx, structuringElementy, 
    window);
  return;
}

template <typename T>
void icarus_signal_processing::Denoising::removeCoherentNoise2D(
  std::vector<std::vector<T>>& waveLessCoherent,
  const std::vector<std::vector<T>>& filteredWaveforms,
  std::vector<std::vector<T>>& morphedWaveforms,
  std::vector<std::vector<T>>& intrinsicRMS,
  ArrayBool& selectVals,
  ArrayBool& roi,
  std::vector<std::vector<T>>& correctedMedians,
  VectorFloat& thresholdVec,
  const char filterName,
  const unsigned int grouping,
  const unsigned int structuringElementx,
  const unsigned int structuringElementy,
  const unsigned int window)
{
  auto numChannels = filteredWaveforms.size();
  auto nTicks = filteredWaveforms.at(0).size();
  auto nGroups = numChannels / grouping;

  // Coherent noise subtracted denoised waveforms
  waveLessCoherent.resize(filteredWaveforms.size());
  for (auto& v : waveLessCoherent) {
    v.resize(filteredWaveforms.at(0).size());
  }

  // Regions to protect waveform from coherent noise subtraction.
  selectVals.resize(filteredWaveforms.size());
  for (auto& v : selectVals) {
    v.resize(filteredWaveforms.at(0).size());
  }

  roi.resize(filteredWaveforms.size());
  for (auto& v : roi) {
    v.resize(filteredWaveforms.at(0).size());
  }

  correctedMedians.resize(nGroups);
  for (auto& v : correctedMedians) {
    v.resize(nTicks);
  }

  intrinsicRMS.resize(nGroups);
  for (auto& v : intrinsicRMS) {
    v.resize(nTicks);
  }

  icarus_signal_processing::Morph2D denoiser;

  std::vector<std::vector<T>> dilation;
  std::vector<std::vector<T>> erosion;
  std::vector<std::vector<T>> average;
  std::vector<std::vector<T>> gradient;

  denoiser.getFilter2D(filteredWaveforms, structuringElementx,
    structuringElementy, dilation, erosion, average, gradient);

  switch (filterName) {
    case 'd':
      getSelectVals(filteredWaveforms.begin(), dilation.begin(), 
        selectVals.begin(), roi.begin(), thresholdVec.begin(), filteredWaveforms.size(), window);
      morphedWaveforms = dilation;
      break;
    case 'e':
      getSelectVals(filteredWaveforms.begin(), erosion.begin(), 
        selectVals.begin(), roi.begin(), thresholdVec.begin(), filteredWaveforms.size(), window);
      morphedWaveforms = erosion;
      break;
    case 'a':
      getSelectVals(filteredWaveforms.begin(), average.begin(), 
        selectVals.begin(), roi.begin(), thresholdVec.begin(), filteredWaveforms.size(), window);
      morphedWaveforms = average;
      break;
    case 'g':
      getSelectVals(filteredWaveforms.begin(), gradient.begin(), 
        selectVals.begin(), roi.begin(), thresholdVec.begin(), filteredWaveforms.size(), window);
      morphedWaveforms = gradient;
      break;
    default:
      getSelectVals(filteredWaveforms.begin(), gradient.begin(), 
        selectVals.begin(), roi.begin(), thresholdVec.begin(), filteredWaveforms.size(), window);
      morphedWaveforms = gradient;
      break;
  }

  for (size_t i=0; i<nTicks; ++i) {
    for (size_t j=0; j<nGroups; ++j) {
      size_t group_start = j * grouping;
      size_t group_end = (j+1) * grouping;
      // Compute median.
      std::vector<T> v;
      for (size_t c=group_start; c<group_end; ++c) {
        if (!selectVals[c][i]) {
          v.push_back(filteredWaveforms[c][i]);
        }
      }
      T median = (T) 0.0;
      if (v.size() > 0) {
        if (v.size() % 2 == 0) {
          const auto m1 = v.begin() + v.size() / 2 - 1;
          const auto m2 = v.begin() + v.size() / 2;
          std::nth_element(v.begin(), m1, v.end());
          const auto e1 = *m1;
          std::nth_element(v.begin(), m2, v.end());
          const auto e2 = *m2;
          median = (e1 + e2) / 2.0;
        } else {
          const auto m = v.begin() + v.size() / 2;
          std::nth_element(v.begin(), m, v.end());
          median = *m;
        }
      }
      correctedMedians[j][i] = median;
      for (size_t k=group_start; k<group_end; ++k) {
        if (!selectVals[k][i]) {
          waveLessCoherent[k][i] = filteredWaveforms[k][i] - median;
        } else {
          waveLessCoherent[k][i] = filteredWaveforms[k][i];
        }
      }
    }
  }

  float rms = 0.0;
  for (size_t i=0; i<nGroups; ++i) {
    for (size_t j=0; j<nTicks; ++j) {
      std::vector<T> v;
      for (size_t k=i*grouping; k<(i+1)*grouping; ++k) {
        v.push_back(waveLessCoherent[k][j]);
      }
      rms = std::sqrt(
        std::inner_product(
          v.begin(), v.end(), v.begin(), 0.) / T(v.size()));
      intrinsicRMS[i][j] = (T) rms;
    }
  }
  return;
}


#endif

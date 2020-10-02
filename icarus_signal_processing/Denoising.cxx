#ifndef __SIGPROC_TOOLS_DENOISING_CXX__
#define __SIGPROC_TOOLS_DENOISING_CXX__

#include "Denoising.h"

#include <chrono>


void icarus_signal_processing::Denoising::getSelectVals(ArrayFloat::const_iterator  waveformsItr,
                                                        ArrayFloat::const_iterator  morphedWaveformsItr,
                                                        ArrayBool::iterator         selectValsItr,
                                                        ArrayBool::iterator         roiItr,
                                                        VectorFloat::const_iterator thresholdItr,
                                                        const unsigned int          numChannels,
                                                        const unsigned int          window)
{
    auto nTicks = waveformsItr[0].size();

    // Set a protection width
    int halfWidth = std::max(int(window/8),1);

    for (size_t i=0; i<numChannels; ++i) 
    {
        std::vector<float> localVec = morphedWaveformsItr[i];

        float median = getMedian(localVec, localVec.size());

        std::vector<float> baseVec;
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

void icarus_signal_processing::Denoising::removeCoherentNoise1D(ArrayFloat::iterator              waveLessCoherentItr,
                                                                ArrayFloat::const_iterator        filteredWaveformsItr,
                                                                ArrayFloat::iterator              morphedWaveformsItr,
                                                                ArrayFloat::iterator              intrinsicRMSItr,
                                                                ArrayBool::iterator               selectValsItr,
                                                                ArrayBool::iterator               roiItr,
                                                                ArrayFloat::iterator              correctedMediansItr,
                                                                FilterFunctionVec::const_iterator filterFunctionsItr,
                                                                VectorFloat::const_iterator       thresholdItr,
                                                                const unsigned int                numChannels,
                                                                const unsigned int                grouping,
                                                                const unsigned int                window)
{
    auto nTicks  = filteredWaveformsItr->size();
    auto nGroups = numChannels / grouping;

    std::chrono::high_resolution_clock::time_point funcStartTime = std::chrono::high_resolution_clock::now();

    std::chrono::high_resolution_clock::time_point morphStart = funcStartTime;

    for(size_t funcIdx = 0; funcIdx < numChannels; funcIdx++)
    {
        const icarus_signal_processing::IMorphologicalFunctions1D* func = (*(filterFunctionsItr + funcIdx)).get();

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

    std::vector<float> v(grouping);

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

            float median = getMedian(v,idxV);

            correctedMediansItr[j][i] = median;
            for (auto k=group_start; k<group_end; ++k) waveLessCoherentItr[k][i] = filteredWaveformsItr[k][i] - median;
        }
    }

    std::chrono::high_resolution_clock::time_point noiseStop = std::chrono::high_resolution_clock::now();

    float rms(0.);
    for (size_t i=0; i<nGroups; ++i) 
    {
        for (size_t j=0; j<nTicks; ++j) 
        {
            size_t idxV(0);
            for (size_t k=i*grouping; k<(i+1)*grouping; ++k) v[idxV++] = waveLessCoherentItr[k][j];
            rms = std::sqrt(std::inner_product(v.begin(), v.begin()+idxV, v.begin(), 0.) / float(v.size()));
            intrinsicRMSItr[i][j] = rms;
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

float icarus_signal_processing::Denoising::getMedian(std::vector<float>& vals, const unsigned int nVals) const
{
    float median(0.);

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

float icarus_signal_processing::Denoising::getMostProbable(std::vector<float>& vals, const unsigned int nVals)
{
    float mostProbable(0.);

    // Do a simple average if only a few bins
    if (nVals < 5)
    {
        mostProbable = std::accumulate(vals.begin(),vals.begin()+nVals,0.) / float(nVals);
    }
    // Otherwise try to form value around MP
    else
    {
        auto minVItr = std::min_element(vals.begin(),vals.begin()+nVals);
        auto maxVItr = std::max_element(vals.begin(),vals.begin()+nVals);

        float minValue = *minVItr;
        float maxValue = *maxVItr;

        size_t numBins = size_t(maxValue - minValue + 1);

        std::fill(fMPVec.begin(),fMPVec.begin() + numBins,0);
  
        for(typename std::vector<float>::iterator vItr = vals.begin(); vItr != vals.begin() + nVals; vItr++)
        {
            int binIdx = int(std::round((*vItr - minValue)));
  
            fMPVec[binIdx]++;
        }
  
        std::vector<int>::iterator maxItr = std::max_element(fMPVec.begin(),fMPVec.begin()+numBins);
  
        int count = *maxItr;
        float   mpVal = float(count) * float(std::distance(fMPVec.begin(),maxItr));

        auto cntItr = maxItr;

        while(cntItr != fMPVec.begin())
        {
            if (*(cntItr - 1) < 0.5 * *maxItr) break;

            count += *(cntItr - 1);
            mpVal += float(*(cntItr - 1)) * float(std::distance(fMPVec.begin(),cntItr - 1));

            cntItr--;
        }

        cntItr = maxItr;

        while(cntItr + 1 != fMPVec.end())
        {
            if (*(cntItr + 1) < 0.5 * *maxItr) break;

            count += *(cntItr + 1);
            mpVal += float(*(cntItr + 1)) * float(std::distance(fMPVec.begin(),cntItr + 1));

            cntItr++;
        }
  
        mostProbable = mpVal / float(count) + minValue;
    }

    return mostProbable;
}


void icarus_signal_processing::Denoising::removeCoherentNoise2D(ArrayFloat::iterator                                       waveLessCoherentItr,
                                                                ArrayFloat::const_iterator                                 filteredWaveformsItr,
                                                                ArrayFloat::iterator                                       morphedWaveformsItr,
                                                                ArrayFloat::iterator                                       intrinsicRMSItr,
                                                                ArrayBool::iterator                                        selectValsItr,
                                                                ArrayBool::iterator                                        roiItr,
                                                                ArrayFloat::iterator                                       correctedMediansItr,
                                                                const icarus_signal_processing::IMorphologicalFunctions2D* filterFunction,
                                                                VectorFloat::const_iterator                                thresholdItr,
                                                                const unsigned int                                         numChannels,
                                                                const unsigned int                                         grouping,
                                                                const unsigned int                                         window)
{
    auto nTicks  = filteredWaveformsItr->size();
    auto nGroups = numChannels / grouping;

    std::chrono::high_resolution_clock::time_point funcStartTime = std::chrono::high_resolution_clock::now();

    std::chrono::high_resolution_clock::time_point morphStart = funcStartTime;

    (*filterFunction)(filteredWaveformsItr, numChannels, morphedWaveformsItr);

    std::chrono::high_resolution_clock::time_point morphStop  = std::chrono::high_resolution_clock::now();
    std::chrono::high_resolution_clock::time_point selStart = morphStop;

    getSelectVals(filteredWaveformsItr, morphedWaveformsItr, selectValsItr, roiItr, thresholdItr, numChannels, window);

    std::chrono::high_resolution_clock::time_point selStop  = std::chrono::high_resolution_clock::now();
    std::chrono::high_resolution_clock::time_point noiseStart = selStop;

    std::vector<float> v(grouping);

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

            float median = getMedian(v,idxV);

            correctedMediansItr[j][i] = median;
            for (auto k=group_start; k<group_end; ++k) waveLessCoherentItr[k][i] = filteredWaveformsItr[k][i] - median;
        }
    }

    std::chrono::high_resolution_clock::time_point noiseStop = std::chrono::high_resolution_clock::now();

    float rms(0.);
    for (size_t i=0; i<nGroups; ++i) 
    {
        for (size_t j=0; j<nTicks; ++j) 
        {
            size_t idxV(0);
            for (size_t k=i*grouping; k<(i+1)*grouping; ++k) v[idxV++] = waveLessCoherentItr[k][j];
            rms = std::sqrt(std::inner_product(v.begin(), v.begin()+idxV, v.begin(), 0.) / float(v.size()));
            intrinsicRMSItr[i][j] = rms;
        }
    }

    std::chrono::high_resolution_clock::time_point funcStopTime = std::chrono::high_resolution_clock::now();
  
    std::chrono::duration<double> funcTime   = std::chrono::duration_cast<std::chrono::duration<double>>(funcStopTime - funcStartTime);
    std::chrono::duration<double> morphTime  = std::chrono::duration_cast<std::chrono::duration<double>>(morphStop - morphStart);
    std::chrono::duration<double> selTime    = std::chrono::duration_cast<std::chrono::duration<double>>(selStop - selStart);
    std::chrono::duration<double> noiseTime  = std::chrono::duration_cast<std::chrono::duration<double>>(noiseStop - noiseStart);
  
    std::cout << "*** Denoising 2D ***  - # channels: " << numChannels << ", ticks: " << nTicks << ", groups: " << nGroups << std::endl;
    std::cout << "                      - morph: " << morphTime.count() << ", sel: " << selTime.count() << ", noise: " << noiseTime.count() << ", total: " << funcTime.count() << std::endl;

    return;
}


#endif

#ifndef __icarus_signal_processing_DECONVOLVE_CXX__
#define __icarus_signal_processing_DECONVOLVE_CXX__

#include "AdaptiveWiener.h"

#endif

void icarus_signal_processing::AdaptiveWiener::filterLee(
    std::vector<std::vector<float>> &deconvolvedWaveform,
    const std::vector<std::vector<float>> &waveLessCoherent,
    const float noiseVar,
    const unsigned int sx,
    const unsigned int sy)
{
    size_t numChannels = waveLessCoherent.size();
    size_t nTicks = waveLessCoherent.at(0).size();
    int xHalfWindowSize(sx / 2);
    int yHalfWindowSize(sy / 2);

    deconvolvedWaveform.resize(numChannels);
    for (size_t i = 0; i < numChannels; ++i)
    {
        deconvolvedWaveform[i].resize(nTicks);
    }

    for (size_t i = 0; i < numChannels; ++i)
    {
        for (size_t j = 0; j < nTicks; ++j)
        {
            // For each center pixel, apply a adaptive local wiener filter.
            int lbx = i - (int)xHalfWindowSize;
            int ubx = i + (int)xHalfWindowSize;
            int lby = j - (int)yHalfWindowSize;
            int uby = j + (int)yHalfWindowSize;
            size_t lowerBoundx = std::max(lbx, 0);
            size_t upperBoundx = std::min(ubx, (int)numChannels);
            size_t lowerBoundy = std::max(lby, 0);
            size_t upperBoundy = std::min(uby, (int)nTicks);
            std::vector<float> x;
            x.reserve(sx * sy);
            std::vector<float> xsq;
            xsq.reserve(sx * sy);
            for (size_t ix = lowerBoundx; ix < upperBoundx; ++ix)
            {
                for (size_t iy = lowerBoundy; iy < upperBoundy; ++iy)
                {
                    x.push_back(waveLessCoherent[ix][iy]);
                    xsq.push_back(waveLessCoherent[ix][iy] * waveLessCoherent[ix][iy]);
                }
            }
            float localMean   = std::accumulate(x.begin(), x.end(), 0.0) / x.size();
            float localSquare = std::accumulate(xsq.begin(), xsq.end(), 0.0) / x.size();
            float localVar    = localSquare - localMean * localMean;

            if (noiseVar > localVar)
            {
                deconvolvedWaveform[i][j] = localMean;
            }
            else
            {
                deconvolvedWaveform[i][j] = localMean + (1 - noiseVar / localVar) *
                                                            (waveLessCoherent[i][j] - localMean);
            }
        }
    }
    return;
}

void icarus_signal_processing::AdaptiveWiener::MMWF(
    std::vector<std::vector<float>> &deconvolvedWaveform,
    const std::vector<std::vector<float>> &waveLessCoherent,
    const float noiseVar,
    const unsigned int sx,
    const unsigned int sy)
{
    size_t numChannels = waveLessCoherent.size();
    size_t nTicks = waveLessCoherent.at(0).size();
    int xHalfWindowSize(sx / 2);
    int yHalfWindowSize(sy / 2);

    deconvolvedWaveform.resize(numChannels);
    for (size_t i = 0; i < numChannels; ++i)
    {
        deconvolvedWaveform[i].resize(nTicks);
    }

    icarus_signal_processing::MiscUtils utils;

    for (size_t i = 0; i < numChannels; ++i)
    {
        for (size_t j = 0; j < nTicks; ++j)
        {
            // For each center pixel, apply a adaptive local wiener filter.
            int lbx = i - (int)xHalfWindowSize;
            int ubx = i + (int)xHalfWindowSize;
            int lby = j - (int)yHalfWindowSize;
            int uby = j + (int)yHalfWindowSize;
            size_t lowerBoundx = std::max(lbx, 0);
            size_t upperBoundx = std::min(ubx, (int)numChannels);
            size_t lowerBoundy = std::max(lby, 0);
            size_t upperBoundy = std::min(uby, (int)nTicks);
            std::vector<float> x;
            x.reserve(sx * sy);
            std::vector<float> xsq;
            xsq.reserve(sx * sy);
            for (size_t ix = lowerBoundx; ix < upperBoundx; ++ix)
            {
                for (size_t iy = lowerBoundy; iy < upperBoundy; ++iy)
                {
                    x.push_back(waveLessCoherent[ix][iy]);
                    xsq.push_back(waveLessCoherent[ix][iy] * waveLessCoherent[ix][iy]);
                }
            }
            float localMean   = std::accumulate(x.begin(), x.end(), 0.0) / x.size();
            float localSquare = std::accumulate(xsq.begin(), xsq.end(), 0.0) / x.size();
            float localVar    = localSquare - localMean * localMean;
            float localMedian = utils.computeMedian(x);

            if (noiseVar > localVar)
            {
                deconvolvedWaveform[i][j] = localMedian;
            }
            else
            {
                deconvolvedWaveform[i][j] = localMedian + (1 - noiseVar / localVar) *
                                                              (waveLessCoherent[i][j] - localMedian);
            }
        }
    }
    return;
}

void icarus_signal_processing::AdaptiveWiener::MMWFStar(
    std::vector<std::vector<float>> &deconvolvedWaveform,
    const std::vector<std::vector<float>> &waveLessCoherent,
    const unsigned int sx,
    const unsigned int sy)
{
    size_t numChannels = waveLessCoherent.size();
    size_t nTicks = waveLessCoherent.at(0).size();
    int xHalfWindowSize(sx / 2);
    int yHalfWindowSize(sy / 2);

    deconvolvedWaveform.resize(numChannels);
    for (size_t i = 0; i < numChannels; ++i)
    {
        deconvolvedWaveform[i].resize(nTicks);
    }

    icarus_signal_processing::MiscUtils utils;

    std::vector<std::vector<float>> localMedians;
    std::vector<std::vector<float>> localVars;

    localMedians.resize(numChannels);
    for (size_t i = 0; i < numChannels; ++i)
    {
        localMedians[i].resize(nTicks);
    }

    localVars.resize(numChannels);
    for (size_t i = 0; i < numChannels; ++i)
    {
        localVars[i].resize(nTicks);
    }

    std::vector<float> localVarTemp;
    localVarTemp.reserve(numChannels * nTicks);

    for (size_t i = 0; i < numChannels; ++i)
    {
        for (size_t j = 0; j < nTicks; ++j)
        {
            // For each center pixel, apply a adaptive local wiener filter.
            int lbx = i - (int)xHalfWindowSize;
            int ubx = i + (int)xHalfWindowSize;
            int lby = j - (int)yHalfWindowSize;
            int uby = j + (int)yHalfWindowSize;
            size_t lowerBoundx = std::max(lbx, 0);
            size_t upperBoundx = std::min(ubx, (int)numChannels);
            size_t lowerBoundy = std::max(lby, 0);
            size_t upperBoundy = std::min(uby, (int)nTicks);
            std::vector<float> x;
            x.reserve(sx * sy);
            std::vector<float> xsq;
            xsq.reserve(sx * sy);
            for (size_t ix = lowerBoundx; ix < upperBoundx; ++ix)
            {
                for (size_t iy = lowerBoundy; iy < upperBoundy; ++iy)
                {
                    x.push_back(waveLessCoherent[ix][iy]);
                    xsq.push_back(waveLessCoherent[ix][iy] * waveLessCoherent[ix][iy]);
                }
            }
            float localMean   = std::accumulate(x.begin(), x.end(), 0.0) / x.size();
            float localSquare = std::accumulate(xsq.begin(), xsq.end(), 0.0) / x.size();
            float localMedian = utils.computeMedian(x);
            float localVar    = localSquare - 2.0 * localMean * localMedian + std::pow(localMean, 2.0);

            localMedians[i][j] = localMedian;
            localVarTemp.push_back(localVar);
            localVars[i][j] = localVar;
        }
    }

    float noiseMedian = utils.computeMedian(localVarTemp);
    std::cout << noiseMedian << std::endl;

    for (size_t i = 0; i < numChannels; ++i)
    {
        for (size_t j = 0; j < nTicks; ++j)
        {
            if (noiseMedian > localVars[i][j])
            {
                deconvolvedWaveform[i][j] = localMedians[i][j];
            }
            else
            {
                deconvolvedWaveform[i][j] = localMedians[i][j] +
                                            (1.0 - noiseMedian / localVars[i][j]) * (waveLessCoherent[i][j] - localMedians[i][j]);
            }
        }
    }

    return;
}

void icarus_signal_processing::AdaptiveWiener::filterLeeEnhanced(
    std::vector<std::vector<float>> &deconvolvedWaveform,
    const std::vector<std::vector<float>> &waveLessCoherent,
    const float noiseVar,
    const unsigned int sx,
    const unsigned int sy,
    const float a,
    const float epsilon)
{
    size_t numChannels = waveLessCoherent.size();
    size_t nTicks = waveLessCoherent.at(0).size();
    int xHalfWindowSize(sx / 2);
    int yHalfWindowSize(sy / 2);

    deconvolvedWaveform.resize(numChannels);
    for (size_t i = 0; i < numChannels; ++i)
    {
        deconvolvedWaveform[i].resize(nTicks);
    }

    for (size_t i = 0; i < numChannels; ++i)
    {
        for (size_t j = 0; j < nTicks; ++j)
        {
            // For each center pixel, apply a adaptive local wiener filter.
            int lbx = i - (int)xHalfWindowSize;
            int ubx = i + (int)xHalfWindowSize;
            int lby = j - (int)yHalfWindowSize;
            int uby = j + (int)yHalfWindowSize;
            size_t lowerBoundx = std::max(lbx, 0);
            size_t upperBoundx = std::min(ubx, (int)numChannels);
            size_t lowerBoundy = std::max(lby, 0);
            size_t upperBoundy = std::min(uby, (int)nTicks);
            std::vector<float> x;
            x.reserve(sx * sy);
            std::vector<float> xsq;
            xsq.reserve(sx * sy);
            std::vector<float> weight;
            weight.reserve(sx * sy);
            for (size_t ix = lowerBoundx; ix < upperBoundx; ++ix)
            {
                for (size_t iy = lowerBoundy; iy < upperBoundy; ++iy)
                {
                    float eps = std::pow(epsilon * std::sqrt(noiseVar), 2);
                    float fsq = std::pow(waveLessCoherent[i][j] - waveLessCoherent[ix][iy], 2.0);
                    float w = 1.0 / (1.0 + a * std::max(eps, fsq));
                    x.push_back(waveLessCoherent[ix][iy]);
                    xsq.push_back(waveLessCoherent[ix][iy] * waveLessCoherent[ix][iy]);
                    weight.push_back(w);
                }
            }
            float normWeight = std::accumulate(weight.begin(), weight.end(), 0.0);
            for (auto &w : weight)
            {
                w = w / normWeight;
            }
            float localMean   = std::inner_product(x.begin(), x.end(), weight.begin(), 0.0) / x.size();
            float localSquare = std::inner_product(xsq.begin(), xsq.end(), weight.begin(), 0.0) / x.size();
            float localVar    = localSquare - localMean * localMean;

            if (noiseVar > localVar)
            {
                deconvolvedWaveform[i][j] = localMean;
            }
            else
            {
                deconvolvedWaveform[i][j] = localMean + (1 - noiseVar / localVar) *
                                                            (waveLessCoherent[i][j] - localMean);
            }
        }
    }
    return;
}

void icarus_signal_processing::AdaptiveWiener::adaptiveROIWiener(
    std::vector<std::vector<float>> &deconvolvedWaveform,
    const std::vector<std::vector<float>> &waveLessCoherent,
    const std::vector<std::vector<bool>> &selectVals,
    const float noiseVar,
    const unsigned int sx,
    const unsigned int sy,
    const float a,
    const float epsilon)
{
    size_t numChannels = waveLessCoherent.size();
    size_t nTicks = waveLessCoherent.at(0).size();
    int xHalfWindowSize(sx / 2);
    int yHalfWindowSize(sy / 2);

    deconvolvedWaveform.resize(numChannels);
    for (size_t i = 0; i < numChannels; ++i)
    {
        deconvolvedWaveform[i].resize(nTicks);
    }

    for (size_t i = 0; i < numChannels; ++i)
    {
        for (size_t j = 0; j < nTicks; ++j)
        {
            // For each center pixel, apply a adaptive local wiener filter.
            int lbx = i - (int)xHalfWindowSize;
            int ubx = i + (int)xHalfWindowSize;
            int lby = j - (int)yHalfWindowSize;
            int uby = j + (int)yHalfWindowSize;
            size_t lowerBoundx = std::max(lbx, 0);
            size_t upperBoundx = std::min(ubx, (int)numChannels);
            size_t lowerBoundy = std::max(lby, 0);
            size_t upperBoundy = std::min(uby, (int)nTicks);
            std::vector<float> x;
            x.reserve(sx * sy);
            std::vector<float> xsq;
            xsq.reserve(sx * sy);
            std::vector<float> weight;
            weight.reserve(sx * sy);
            for (size_t ix = lowerBoundx; ix < upperBoundx; ++ix)
            {
                for (size_t iy = lowerBoundy; iy < upperBoundy; ++iy)
                {
                    float eps = std::pow(epsilon * std::sqrt(noiseVar), 2);
                    float fsq = std::pow(waveLessCoherent[i][j] - waveLessCoherent[ix][iy], 2.0);
                    float w = 1.0 / (1.0 + a * std::max(eps, fsq));
                    x.push_back(waveLessCoherent[ix][iy]);
                    xsq.push_back(waveLessCoherent[ix][iy] * waveLessCoherent[ix][iy]);
                    weight.push_back(w);
                }
            }
            float normWeight = std::accumulate(weight.begin(), weight.end(), 0.0);
            for (auto &w : weight)
            {
                w = w / normWeight;
            }
            float localMean   = std::inner_product(x.begin(), x.end(), weight.begin(), 0.0) / x.size();
            float localSquare = std::inner_product(xsq.begin(), xsq.end(), weight.begin(), 0.0) /x.size();
            float localVar    = localSquare - localMean * localMean;

            if (noiseVar > localVar)
            {
                deconvolvedWaveform[i][j] = localMean;
            }
            else if (selectVals[i][j])
            {
                deconvolvedWaveform[i][j] = waveLessCoherent[i][j];
            }
            else
            {
                deconvolvedWaveform[i][j] = localMean + (1 - noiseVar / localVar) *
                                                            (waveLessCoherent[i][j] - localMean);
            }
        }
    }
    return;
}

void icarus_signal_processing::AdaptiveWiener::sigmaFilter(
    std::vector<std::vector<float>> &deconvolvedWaveform,
    const std::vector<std::vector<float>> &waveLessCoherent,
    const float noiseVar,
    const unsigned int sx,
    const unsigned int sy,
    const unsigned int K,
    const float sigmaFactor)
{
    size_t numChannels = waveLessCoherent.size();
    size_t nTicks = waveLessCoherent.at(0).size();
    int xHalfWindowSize(sx / 2);
    int yHalfWindowSize(sy / 2);

    deconvolvedWaveform.resize(numChannels);
    for (size_t i = 0; i < numChannels; ++i)
    {
        deconvolvedWaveform[i].resize(nTicks);
    }

    for (size_t i = 0; i < numChannels; ++i)
    {
        for (size_t j = 0; j < nTicks; ++j)
        {
            // For each center pixel, apply a adaptive local wiener filter.
            int lbx = i - (int)xHalfWindowSize;
            int ubx = i + (int)xHalfWindowSize;
            int lby = j - (int)yHalfWindowSize;
            int uby = j + (int)yHalfWindowSize;
            size_t lowerBoundx = std::max(lbx, 0);
            size_t upperBoundx = std::min(ubx, (int)numChannels);
            size_t lowerBoundy = std::max(lby, 0);
            size_t upperBoundy = std::min(uby, (int)nTicks);
            std::vector<float> x;
            x.reserve(sx * sy);
            for (size_t ix = lowerBoundx; ix < upperBoundx; ++ix)
            {
                for (size_t iy = lowerBoundy; iy < upperBoundy; ++iy)
                {
                    if (std::abs(waveLessCoherent[ix][iy]) < sigmaFactor * noiseVar)
                    {
                        x.push_back(waveLessCoherent[ix][iy]);
                    }
                }
            }

            float localMean = std::accumulate(x.begin(), x.end(), 0.0) / x.size();

            if (x.size() > K)
            {
                deconvolvedWaveform[i][j] = localMean;
            }
            else
            {
                deconvolvedWaveform[i][j] = waveLessCoherent[i][j];
            }
        }
    }
    return;
}

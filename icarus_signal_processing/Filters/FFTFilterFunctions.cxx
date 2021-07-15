#ifndef __SIGPROC_TOOLS_FFTFILTERFUNCTIONS_CXX__
#define __SIGPROC_TOOLS_FFTFILTERFUNCTIONS_CXX__

#define ICARUSFFT_FLOAT

#include "FFTFilterFunctions.h"
#include "icarus_signal_processing/Filters/ICARUSFFT.h"

#include <iostream>
#include <algorithm>
#include <iterator>
#include <numeric>
#include <cmath>
#include <functional>


namespace icarus_signal_processing 
{

HighPassFFTFilter::HighPassFFTFilter(const double& sigma, const double& offset)
{
    fFFT = std::make_unique<icarus_signal_processing::ICARUSFFT<float>>();

        fHighPassFilterKernel.resize(4096);

        std::fill(fHighPassFilterKernel.begin(),fHighPassFilterKernel.end(),std::complex<float>(1.,0.));

        for(int binIdx = 0; binIdx < int(offset); binIdx++)
        {
            float expVal = -pow((float(binIdx) - offset)/sigma,2.);
            float binVal = exp(expVal);

            fHighPassFilterKernel[binIdx] = std::complex<float>(binVal,0.);
        }

    return;
}

HighPassFFTFilter::~HighPassFFTFilter()
{
    return;
}

void HighPassFFTFilter::operator()(Waveform<float>& waveformVec) const
{
    highPassFilter(waveformVec);

    return;
}

void HighPassFFTFilter::operator()(Waveform<double>& waveformVec) const
{
    Waveform<float> locWaveform(waveformVec.size());

    std::copy(waveformVec.begin(), waveformVec.end(), locWaveform.begin());

    highPassFilter(locWaveform);

    std::copy(locWaveform.begin(),locWaveform.end(),waveformVec.begin());

    return;
}

void icarus_signal_processing::HighPassFFTFilter::highPassFilter(Waveform<float>& inputWaveform) const
{
    fFFT->convolute(inputWaveform, fHighPassFilterKernel, 0);

    return;
}

void LowPassFFTFilter::operator()(Waveform<float>& waveformVec) const 
{
    lowPassFilter(waveformVec);

    return;
}

void LowPassFFTFilter::operator()(Waveform<double>& waveformVec) const
{
    Waveform<float> locWaveform(waveformVec.size());

    std::copy(waveformVec.begin(), waveformVec.end(), locWaveform.begin());

    lowPassFilter(locWaveform);

    std::copy(locWaveform.begin(), locWaveform.end(), waveformVec.begin());

    return;
}


void icarus_signal_processing::LowPassFFTFilter::lowPassFilter(Waveform<float>& inputWaveform) const
{
    return;
}

WindowFFTFilter::WindowFFTFilter(const std::pair<double,double>& sigmaPair, 
                                 const std::pair<double,double>& offsetPair)
{
    fFFT = std::make_unique<icarus_signal_processing::ICARUSFFT<float>>();

    fWindowFilterKernel.resize(4096,std::complex<float>(0.,0.));

    int lowOffset = std::max(int(offsetPair.first),     1);
    int hiOffset  = std::min(int(offsetPair.second), 4096);

    // We skip the zero bin to make sure it is set to zero
    for(int binIdx = 1; binIdx < lowOffset; binIdx++)
    {
        float expVal = -pow((float(binIdx) - offsetPair.first)/sigmaPair.first,2.);
        float binVal = exp(expVal);

        fWindowFilterKernel[binIdx] = std::complex<double>(binVal,0.);
    }

    if (lowOffset < hiOffset) 
        std::fill(fWindowFilterKernel.begin() + lowOffset, fWindowFilterKernel.begin() + hiOffset, std::complex<double>(1.,0.));

    int endOffset = std::min(5 * int(sigmaPair.second), 4096);

    for(int binIdx = hiOffset; binIdx < endOffset; binIdx++)
    {
        float expVal = -pow((float(binIdx) - offsetPair.second)/sigmaPair.second,2.);
        float binVal = exp(expVal);

        fWindowFilterKernel[binIdx] = std::complex<double>(binVal,0.);
    }

    return;
}

WindowFFTFilter::~WindowFFTFilter()
{
    return;
}

void WindowFFTFilter::operator()(Waveform<float>& waveformVec) const
{
    windowFilter(waveformVec);

    return;
}

void WindowFFTFilter::operator()(Waveform<double>& waveformVec) const
{
    Waveform<float> locWaveform(waveformVec.size());

    std::copy(waveformVec.begin(), waveformVec.end(), locWaveform.begin());

    windowFilter(locWaveform);

    std::copy(locWaveform.begin(),locWaveform.end(),waveformVec.begin());

    return;
}

void WindowFFTFilter::windowFilter(Waveform<float>& inputWaveform) const
{
    fFFT->convolute(inputWaveform, fWindowFilterKernel, 0);

    return;
}

LowPassButterworthFilter::LowPassButterworthFilter(unsigned int threshold,
                                                   unsigned int order,
                                                   unsigned int size) :
                fThreshold(threshold),
                fOrder(order),
                fFrequencySize(size)
{

    if (fThreshold < 1 or fThreshold > size / 2) 
            throw std::runtime_error("FrequencyFilters: Thresholding frequency bin must be at least one and less than half of total bin count.");

    fFFT = std::make_unique<icarus_signal_processing::ICARUSFFT<float>>();

    fFilterVec.resize(fFrequencySize);

    float d0 = (float) fThreshold;

    for (size_t i = 0; i < fFrequencySize / 2; i++) 
    {
        float d           = (float) i;
        float filterValue = 1.0 / (1.0 + std::pow( d / d0, 2 * fOrder));

        fFilterVec[i]                  = std::complex<float>(filterValue,0.);
        fFilterVec[fFrequencySize-i-1] = std::complex<float>(filterValue,0.);
    }

    // Ad Hoc compensation for odd length frequency arrays.
    fFilterVec[fFrequencySize/2] = fFilterVec[fFrequencySize/2-1];

    return;
}

LowPassButterworthFilter::~LowPassButterworthFilter()
{
    return;
}

void LowPassButterworthFilter::operator()(Waveform<float>& waveformVec) const
{
    lowPassButterworthFilter(waveformVec);

    return;
}

void LowPassButterworthFilter::operator()(Waveform<double>& waveformVec) const
{
    Waveform<float> locWaveform(waveformVec.size());

    std::copy(waveformVec.begin(), waveformVec.end(), locWaveform.begin());

    lowPassButterworthFilter(locWaveform);

    std::copy(locWaveform.begin(),locWaveform.end(),waveformVec.begin());

    return;
}

void LowPassButterworthFilter::lowPassButterworthFilter(Waveform<float>& inputWaveform) const
{
    fFFT->convolute(inputWaveform, fFilterVec, 0);

    return;
}

HighPassButterworthFilter::HighPassButterworthFilter(unsigned int threshold,
                                                     unsigned int order,
                                                     unsigned int size) :
                fThreshold(threshold),
                fOrder(order),
                fFrequencySize(size)
{

    if (fThreshold < 1 or fThreshold > size / 2) 
            throw std::runtime_error("FrequencyFilters: Thresholding frequency bin must be at least one and less than half of total bin count.");
            
    fFFT = std::make_unique<icarus_signal_processing::ICARUSFFT<float>>();

    fFilterVec.resize(fFrequencySize);

    float d0 = (float) fThreshold;

    for (size_t i = 0; i < fFrequencySize / 2; i++) 
    {
        float d           = (float) i;
        float filterValue = 1.0 - 1.0 / (1.0 + std::pow( d / d0, 2 * fOrder));

        fFilterVec[i]                  = std::complex<float>(filterValue,0.);
        fFilterVec[fFrequencySize-i-1] = std::complex<float>(filterValue,0.);
    }

    // Ad Hoc compensation for odd length frequency arrays.
    fFilterVec[fFrequencySize/2] = fFilterVec[fFrequencySize/2-1];

    return;
}

HighPassButterworthFilter::~HighPassButterworthFilter()
{
    return;
}

void HighPassButterworthFilter::operator()(Waveform<float>& waveformVec) const
{
    highPassButterworthFilter(waveformVec);

    return;
}

void HighPassButterworthFilter::operator()(Waveform<double>& waveformVec) const
{
    Waveform<float> locWaveform(waveformVec.size());

    std::copy(waveformVec.begin(), waveformVec.end(), locWaveform.begin());

    highPassButterworthFilter(locWaveform);

    std::copy(locWaveform.begin(),locWaveform.end(),waveformVec.begin());

    return;
}

void HighPassButterworthFilter::highPassButterworthFilter(Waveform<float>& inputWaveform) const
{
    fFFT->convolute(inputWaveform, fFilterVec, 0);

    return;
}

}


#endif

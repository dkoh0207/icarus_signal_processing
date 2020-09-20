#ifndef __SIGPROC_TOOLS_FFTFILTERFUNCTIONS_CXX__
#define __SIGPROC_TOOLS_FFTFILTERFUNCTIONS_CXX__

#include "FFTFilterFunctions.h"
#include "icarus_signal_processing/ICARUSFFT.h"

#include <iostream>
#include <algorithm>
#include <iterator>
#include <numeric>
#include <cmath>
#include <functional>


namespace icarus_signal_processing 
{

HighPassFFTFilter::HighPassFFTFilter(const std::vector<double>& sigmaVec, const std::vector<double>& offsetVec)
{
    fFFT = std::make_unique<icarus_signal_processing::ICARUSFFT<double>>();

    fHighPassFilterKernels.resize(3);

    for(size_t planeIdx = 0; planeIdx < sigmaVec.size(); planeIdx++)
    {
        KernelVec& highPassFilterKernel = fHighPassFilterKernels[planeIdx];

        highPassFilterKernel.resize(4096);

        std::fill(highPassFilterKernel.begin(),highPassFilterKernel.end(),std::complex<double>(1.,0.));

        for(int binIdx = 0; binIdx < int(offsetVec[planeIdx]); binIdx++)
        {
            float expVal = -pow((float(binIdx) - offsetVec[planeIdx])/sigmaVec[planeIdx],2.);
            float binVal = exp(expVal);

            highPassFilterKernel[binIdx] = std::complex<double>(binVal,0.);

            std::cout << "--Bin: " << binIdx << ", expVal: " << expVal << ", binVal: " << binVal << std::endl;
        }
    }

    return;
}

HighPassFFTFilter::~HighPassFFTFilter()
{
    return;
}

void HighPassFFTFilter::operator()(Waveform<float>& waveformVec, int plane) const
{
    Waveform<double> locWaveform(waveformVec.size());

    std::copy(waveformVec.begin(), waveformVec.end(), locWaveform.begin());

    highPassFilter(locWaveform, plane);

    std::copy(locWaveform.begin(),locWaveform.end(),waveformVec.begin());

    return;
}

void HighPassFFTFilter::operator()(Waveform<double>& waveformVec, int plane) const
{
    highPassFilter(waveformVec, plane);
}

void icarus_signal_processing::HighPassFFTFilter::highPassFilter(Waveform<double>& inputWaveform, int plane) const
{
    fFFT->convolute(inputWaveform, fHighPassFilterKernels[plane], 0);

    return;
}

void LowPassFFTFilter::operator()(Waveform<float>& waveformVec, int plane) const 
{
    Waveform<double> locWaveform(waveformVec.size());

    std::copy(waveformVec.begin(), waveformVec.end(), locWaveform.begin());

    lowPassFilter(locWaveform, plane);

    std::copy(locWaveform.begin(), locWaveform.end(), waveformVec.begin());

    return;
}

void LowPassFFTFilter::operator()(Waveform<double>& waveformVec, int plane) const
{
    lowPassFilter(waveformVec, plane);
}


void icarus_signal_processing::LowPassFFTFilter::lowPassFilter(Waveform<double>& inputWaveform, int plane) const
{
    return;
}

WindowFFTFilter::WindowFFTFilter(const std::vector<std::pair<double,double>>& sigmaPairVec, 
                                 const std::vector<std::pair<double,double>>& offsetPairVec)
{
    fFFT = std::make_unique<icarus_signal_processing::ICARUSFFT<double>>();

    fWindowFilterKernels.resize(3);

    for(size_t planeIdx = 0; planeIdx < sigmaPairVec.size(); planeIdx++)
    {
        KernelVec& windowFilterKernel = fWindowFilterKernels[planeIdx];

        windowFilterKernel.resize(4096,std::complex<double>(0.,0.));

        int lowOffset = std::min(int(offsetPairVec[planeIdx].first),  4096);
        int hiOffset  = std::min(int(offsetPairVec[planeIdx].second), 4096);

        // We skip the zero bin to make sure it is set to zero
        for(int binIdx = 1; binIdx < lowOffset; binIdx++)
        {
            float expVal = -pow((float(binIdx) - offsetPairVec[planeIdx].first)/sigmaPairVec[planeIdx].first,2.);
            float binVal = exp(expVal);

            windowFilterKernel[binIdx] = std::complex<double>(binVal,0.);

            std::cout << "--Bin: " << binIdx << ", expVal: " << expVal << ", binVal: " << binVal << std::endl;
        }

        if (lowOffset < hiOffset) 
            std::fill(windowFilterKernel.begin() + lowOffset, windowFilterKernel.begin() + hiOffset, std::complex<double>(1.,0.));

        int endOffset = std::min(5 * int(sigmaPairVec[planeIdx].second), 4096);

        for(int binIdx = hiOffset; binIdx < endOffset; binIdx++)
        {
            float expVal = -pow((float(binIdx) - offsetPairVec[planeIdx].second)/sigmaPairVec[planeIdx].second,2.);
            float binVal = exp(expVal);

            windowFilterKernel[binIdx] = std::complex<double>(binVal,0.);

            std::cout << "--Bin: " << binIdx << ", expVal: " << expVal << ", binVal: " << binVal << std::endl;

        }
    }

    return;
}

WindowFFTFilter::~WindowFFTFilter()
{
    return;
}

void WindowFFTFilter::operator()(Waveform<float>& waveformVec, int plane) const
{
    Waveform<double> locWaveform(waveformVec.size());

    std::copy(waveformVec.begin(), waveformVec.end(), locWaveform.begin());

    windowFilter(locWaveform, plane);

    std::copy(locWaveform.begin(),locWaveform.end(),waveformVec.begin());

    return;
}

void WindowFFTFilter::operator()(Waveform<double>& waveformVec, int plane) const
{
    windowFilter(waveformVec, plane);
}

void WindowFFTFilter::windowFilter(Waveform<double>& inputWaveform, int plane) const
{
    fFFT->convolute(inputWaveform, fWindowFilterKernels[plane], 0);

    return;
}


}


#endif

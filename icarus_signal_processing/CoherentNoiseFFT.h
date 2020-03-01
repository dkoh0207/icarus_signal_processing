#ifndef CoherentNoiseFFT_H
#define CoherentNoiseFFT_H
////////////////////////////////////////////////////////////////////////
//
// Class:       CoherentNoiseFFT
// Module Type: producer
// File:        CoherentNoiseFFT.h
//
//              This module provides some basic Fast Fourier Transform
//              algorithms for operating on RawDigit waveforms
//
// Configuration parameters:
//
// FillHistograms        - Turn on histogram filling for diagnostics
// RunFFTInputWires      - FFT analyze the input RawDigits if true - diagnostics
// RunFFTCorrectedWires  - FFT analyze the output RawDigits if true - diagnostics
//
//
// Created by Tracy Usher (usher@slac.stanford.edu) on January 6, 2016
//
////////////////////////////////////////////////////////////////////////

#include "icarus_signal_processing/ICARUSFFT.h"
#include "icarus_signal_processing/WaveformTools.h"

namespace icarus_signal_processing
{
    
template <class T> class CoherentNoiseFFT
{
public:

    // Copnstructors, destructor.
    CoherentNoiseFFT(size_t = 4096);
    ~CoherentNoiseFFT();
    
    void getFFTCorrection(std::vector<T>&, double) const;
    
private:
    
    using ICARUSFFTDEF = icarus_signal_processing::ICARUSFFT<double>;
    using FFTPointer = std::unique_ptr<ICARUSFFTDEF>;

    size_t                            fNumberTimeSamples;     //< Number time samples for FFT

    // Try to optimize the filter FFT function with static memory...
    std::vector<float>                fFFTInputVec;
    std::vector<std::complex<float>>  fFFTOutputVec;
    std::vector<float>                fPowerVec;

    WaveformTools<double>             fWaveformTools;
    FFTPointer                        fFFT;                   //< Object to handle thread safe FFT
};
    
} // end caldata namespace

#endif

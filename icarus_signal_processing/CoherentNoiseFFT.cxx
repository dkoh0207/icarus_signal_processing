
#include "icarus_signal_processing/CoherentNoiseFFT.h"

#include <complex.h>
#include <cmath>
#include <algorithm>
#include <complex>

namespace icarus_signal_processing
{

//----------------------------------------------------------------------------
/// Constructor.
///
/// Arguments:
///
/// pset - Fcl parameters.
///
template <class T> CoherentNoiseFFT<T>::CoherentNoiseFFT(size_t numberTimeSamples) :
    fNumberTimeSamples(numberTimeSamples),
    fWaveformTools(numberTimeSamples)
{ 
    fFFT = std::make_unique<ICARUSFFTDEF>(numberTimeSamples);

    return;
}

//----------------------------------------------------------------------------
/// Destructor.
template <class T> CoherentNoiseFFT<T>::~CoherentNoiseFFT()
{}
    
template <class T> void CoherentNoiseFFT<T>::getFFTCorrection(std::vector<T>& corValVec, double minPowerThreshold) const
{
    // This version will take FFT of input waveform and then remove bins in the time domain with a power less
    // than the threshold input above.
    size_t const fftDataSize = corValVec.size();
    
    std::vector<std::complex<T>> fftOutputVec(corValVec.size());

    // Now we compute the convolution kernel which is a straigtforward operation
    fFFT->forwardFFT(corValVec, fftOutputVec);
    
    size_t halfFFTDataSize(fftDataSize/2 + 1);
    
    std::vector<T> powerVec(halfFFTDataSize);
    
    std::transform(fftOutputVec.begin(), fftOutputVec.begin() + halfFFTDataSize, powerVec.begin(), [](const auto& val){return std::abs(val);});
    
    // Want the first derivative
    std::vector<T> firstDerivVec(powerVec.size());
    
    fWaveformTools.firstDerivative(powerVec, firstDerivVec);
    
    // Find the peaks
    typename icarus_signal_processing::WaveformTools<T>::PeakTupleVec peakTupleVec;
    
    fWaveformTools.findPeaks(firstDerivVec.begin(),firstDerivVec.end(),peakTupleVec,minPowerThreshold,0);
    
    if (!peakTupleVec.empty())
    {
        for(const auto& peakTuple : peakTupleVec)
        {
            size_t startTick = std::get<0>(peakTuple);
            size_t stopTick  = std::get<2>(peakTuple);
            
            if (stopTick > startTick)
            {
                std::complex<T> slope = (fftOutputVec[stopTick] - fftOutputVec[startTick]) / T(stopTick - startTick);
                
                for(size_t tick = startTick; tick < stopTick; tick++)
                {
                    std::complex<T> interpVal = fftOutputVec[startTick] + T(tick - startTick) * slope;
                    
                    fftOutputVec[tick]                   = interpVal;
                    //fftOutputVec[fftDataSize - tick - 1] = interpVal;
                }
            }
        }
        
        std::vector<T> tmpVec(corValVec.size());
        
        fFFT->inverseFFT(fftOutputVec,tmpVec);
        
        std::transform(corValVec.begin(),corValVec.end(),tmpVec.begin(),corValVec.begin(),std::minus<T>());
    }
    
    return;
}

}

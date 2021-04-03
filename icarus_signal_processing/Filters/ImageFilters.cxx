#ifndef __SIGPROC_TOOLS_ImageFilters_CXX__
#define __SIGPROC_TOOLS_ImageFilters_CXX__

#include "ImageFilters.h"


namespace icarus_signal_processing 
{

void ImageFilters::operator()(Waveform2D<float>& waveformVec) const
{
    for(auto& waveform : waveformVec) (*fFilterFunction)(waveform);

    return;
}

void ImageFilters::operator()(Waveform2D<double>& waveformVec) const
{
    for(auto& waveform : waveformVec) (*fFilterFunction)(waveform);

    return;
}

}

#endif

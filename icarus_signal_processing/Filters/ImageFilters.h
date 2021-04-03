/**
 * \file ImageFilters.h
 *
 * \ingroup sigproc_tools
 * 
 * \brief Class def header for a class ImageFilters
 *
 * @author koh0207
 */

/** \addtogroup sigproc_tools

    @{*/
#ifndef __SIGPROC_TOOLS_ImageFilters_H__
#define __SIGPROC_TOOLS_ImageFilters_H__

#include "FFTFilterFunctions.h"

namespace icarus_signal_processing 
{

    /**
     \class ImageFilters
     User defined class ImageFilters ... these comments are used to generate
     doxygen documentation!
  */
    class ImageFilters
    {
    public:
        template <class T> using Waveform   = std::vector<T>;
        template <class T> using Waveform2D = std::vector<std::vector<T>>;

        /// Default constructor
        ImageFilters(const IFFTFilterFunction* function) : fFilterFunction(function) {}

        /// Default destructor
        ~ImageFilters() {}
 
        /**
        *  @brief Interface functions to perform the waveform filtering
        */
        void operator()(Waveform2D<float>&)  const;
        void operator()(Waveform2D<double>&) const;

    private:
        const IFFTFilterFunction* fFilterFunction;
    };
} // namespace sigproc_tools

#endif
/** @} */ // end of doxygen group

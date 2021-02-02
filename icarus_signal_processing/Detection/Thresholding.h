/**
 * \file LocalThresholding.h
 *
 * \ingroup icuars_signal_processing
 * 
 * \brief Class def header for a class LocalThresholding
 *
 * @author koh0207
 */

/** \addtogroup icuars_signal_processing

    @{*/
#ifndef __icuars_signal_processing_LOCALTHRESHOLDING_H__
#define __icuars_signal_processing_LOCALTHRESHOLDING_H__

#include <vector>
#include <cmath>
#include <algorithm>
#include <functional>
#include <iostream>

namespace icuars_signal_processing
{

    /**
     \class Thresholding
     User defined class LocalThresholding ... these comments are used to generate
     doxygen documentation!
  */

    template <class T>
    using Array2D = std::vector<std::vector<T>>;

    class Thresholding
    {

    public:
        /// Default constructor
        Thresholding() {}

        void Niblack(const Array2D<float> &waveform2D,
                     Array2D<bool> &outputROI,
                     const float k,
                     const int sx,
                     const int sy) const;

        void Sauvola(const Array2D<float> &waveform2D,
                     Array2D<bool> &outputROI,
                     const float k,
                     const float R,
                     const int sx,
                     const int sy) const;

        void Otsu(const Array2D<float> &waveform2D,
                  Array2D<bool> &outputBinary) const;

        void globalMean(const Array2D<float> &waveform2D,
                        Array2D<bool> &outputBinary,
                        float k = 2.0) const;

        void computeIntegralImage(const Array2D<float> &waveform2D,
                                  Array2D<float> &integral2D) const;

        /// Default destructor
        ~Thresholding() {}
    };
} // namespace icuars_signal_processing

#endif
/** @} */ // end of doxygen group

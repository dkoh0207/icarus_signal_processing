/**
 * \file FindROI2D.h
 *
 * \ingroup icarus_signal_processing
 * 
 * \brief Class def header for a class FindROI2D
 *
 * @author koh0207
 */

/** \addtogroup icarus_signal_processing

    @{*/
#ifndef __icarus_signal_processing_FINDROI2D_H__
#define __icarus_signal_processing_FINDROI2D_H__


#include <vector>
#include <algorithm>
#include <cmath>
#include <iostream>
#include <numeric>
#include <iterator>

namespace icarus_signal_processing {

  /**
     \class FindROI2D
     User defined class FindROI2D ... these comments are used to generate
     doxygen documentation!
  */
  template <class T> using Array2D = std::vector<std::vector<T>>;

  class FindROI2D{
    
  public:
    
    /// Default constructor
    FindROI2D(){}

    void applyChainFilter(
        const Array2D<float>& waveform2D,
        Array2D<float>& fullEvent,
        Array2D<bool>& outputROI,
        Array2D<float>& waveLessCoherent,
        Array2D<float>& morphedWaveform2D,
        Array2D<float>& finalErosion2D,
        // Default Parameters
        size_t FREQUENCY_THRESHOLD,
        size_t FREQUENCY_FILTER_SMOOTHNESS_ORDER,
        size_t FREQUENCY_FILTER_MODE,

        char MORPHOLOGICAL_FILTER_NAME,
        const unsigned int CHANNEL_GROUPING,
        const unsigned int STRUCTURING_ELEMENT_X,
        const unsigned int STRUCTURING_ELEMENT_Y,
        const unsigned int ROI_EXPAND_WINDOW_SIZE,
        const float MORPHOLOGICAL_THRESHOLD_FACTOR,

        const size_t THETASTEPS,
        const unsigned int HOUGH_THRESHOLD,
        const unsigned int NMS_WINDOW_SIZE,
        const unsigned int ANGLE_WINDOW,

        // float NOISE_VARIANCE = 20.0;
        const unsigned int ADFILTER_SX,
        const unsigned int ADFILTER_SY,

        const unsigned int BINARY_CLOSING_SX,
        const unsigned int BINARY_CLOSING_SY,

        const float GLOBAL_THRESHOLDING_FACTOR) const;

    void applyCannyFilter(
        const Array2D<float>& waveform2D,
        Array2D<float>& fullEvent,
        Array2D<bool>& outputROI,
        Array2D<float>& waveLessCoherent,
        Array2D<float>& morphedWaveform2D,
        // Default Parameters
        size_t FREQUENCY_THRESHOLD,
        size_t FREQUENCY_FILTER_SMOOTHNESS_ORDER,
        size_t FREQUENCY_FILTER_MODE,

        char MORPHOLOGICAL_FILTER_NAME,
        const unsigned int CHANNEL_GROUPING,
        const unsigned int STRUCTURING_ELEMENT_X,
        const unsigned int STRUCTURING_ELEMENT_Y,
        const unsigned int ROI_EXPAND_WINDOW_SIZE,
        const float MORPHOLOGICAL_THRESHOLD_FACTOR,

        const size_t THETASTEPS,
        const unsigned int HOUGH_THRESHOLD,
        const unsigned int NMS_WINDOW_SIZE,
        const unsigned int ANGLE_WINDOW,

        // float NOISE_VARIANCE = 20.0;
        const unsigned int ADFILTER_SX,
        const unsigned int ADFILTER_SY,
        const float sigma_x, 
        const float sigma_y, 
        const float sigma_r, 
        const float lowThreshold,
        const float highThreshold,

        const unsigned int BINARY_CLOSING_SX,
        const unsigned int BINARY_CLOSING_SY) const;
    
    /// Default destructor
    ~FindROI2D(){}
    
  };
}

#endif
/** @} */ // end of doxygen group 


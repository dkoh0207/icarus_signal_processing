/**
 * \file Deconvolve.h
 *
 * \ingroup icarus_signal_processing
 *
 * \brief Class def header for a class Deconvolve
 *
 * @author koh0207
 */

/** \addtogroup icarus_signal_processing

    @{*/
#ifndef __icarus_signal_processing_DECONVOLVE_H__
#define __icarus_signal_processing_DECONVOLVE_H__

#include <vector>
#include <iostream>
#include <string>
#include <algorithm>
#include <iterator>
#include <numeric>
#include <cmath>
#include <functional>
#include "MiscUtils.h"

namespace icarus_signal_processing {

  /**
     \class Deconvolve
     User defined class Deconvolve ... these comments are used to generate
     doxygen documentation!
  */
  class AdaptiveWiener{

    public:

      /// Default constructor
      AdaptiveWiener(){}

      void filterLee(
        std::vector<std::vector<float>>&,
        const std::vector<std::vector<float>>&,
        const float,
        const unsigned int,
        const unsigned int
      );

      void MMWF(
        std::vector<std::vector<float>>&,
        const std::vector<std::vector<float>>&,
        const float,
        const unsigned int,
        const unsigned int
      );

      void MMWFStar(
        std::vector<std::vector<float>>&,
        const std::vector<std::vector<float>>&,
        const unsigned int,
        const unsigned int
      );

      void filterLeeEnhanced(
        std::vector<std::vector<float>>&,
        const std::vector<std::vector<float>>&,
        const float,
        const unsigned int,
        const unsigned int,
        const float,
        const float
      );

      void adaptiveROIWiener(
        std::vector<std::vector<float>>&,
        const std::vector<std::vector<float>>&,
        const std::vector<std::vector<bool>>&,
        const float,
        const unsigned int,
        const unsigned int,
        const float,
        const float
      );

      void sigmaFilter(
        std::vector<std::vector<float>>&,
        const std::vector<std::vector<float>>&,
        const float,
        const unsigned int,
        const unsigned int,
        const unsigned int,
        const float);

      /// Default destructor
      ~AdaptiveWiener(){}

    private:
  };
}

#endif
/** @} */ // end of doxygen group

/**
 * \file Deconvolve.h
 *
 * \ingroup sigproc_tools
 *
 * \brief Class def header for a class Deconvolve
 *
 * @author koh0207
 */

/** \addtogroup sigproc_tools

    @{*/
#ifndef __SIGPROC_TOOLS_DECONVOLVE_H__
#define __SIGPROC_TOOLS_DECONVOLVE_H__

#include <vector>
#include <iostream>
#include <string>
#include <algorithm>
#include <iterator>
#include <numeric>
#include <cmath>
#include <functional>
#include "MiscUtils.h"

namespace sigproc_tools {

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
        std::vector<std::vector<short>>&,
        const std::vector<std::vector<short>>&,
        const float,
        const unsigned int,
        const unsigned int
      );

      void filterLee(
        std::vector<std::vector<float>>&,
        const std::vector<std::vector<float>>&,
        const float,
        const unsigned int,
        const unsigned int
      );

      void filterLee(
        std::vector<std::vector<double>>&,
        const std::vector<std::vector<double>>&,
        const float,
        const unsigned int,
        const unsigned int
      );


      void MMWF(
        std::vector<std::vector<short>>&,
        const std::vector<std::vector<short>>&,
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

      void MMWF(
        std::vector<std::vector<double>>&,
        const std::vector<std::vector<double>>&,
        const float,
        const unsigned int,
        const unsigned int
      );


      void MMWFStar(
        std::vector<std::vector<short>>&,
        const std::vector<std::vector<short>>&,
        const unsigned int,
        const unsigned int
      );

      void MMWFStar(
        std::vector<std::vector<float>>&,
        const std::vector<std::vector<float>>&,
        const unsigned int,
        const unsigned int
      );

      void MMWFStar(
        std::vector<std::vector<double>>&,
        const std::vector<std::vector<double>>&,
        const unsigned int,
        const unsigned int
      );


      void filterLeeEnhanced(
        std::vector<std::vector<short>>&,
        const std::vector<std::vector<short>>&,
        const float,
        const unsigned int,
        const unsigned int,
        const float,
        const float
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

      void filterLeeEnhanced(
        std::vector<std::vector<double>>&,
        const std::vector<std::vector<double>>&,
        const float,
        const unsigned int,
        const unsigned int,
        const float,
        const float
      );


      void adaptiveROIWiener(
        std::vector<std::vector<short>>&,
        const std::vector<std::vector<short>>&,
        const std::vector<std::vector<bool>>&,
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

      void adaptiveROIWiener(
        std::vector<std::vector<double>>&,
        const std::vector<std::vector<double>>&,
        const std::vector<std::vector<bool>>&,
        const float,
        const unsigned int,
        const unsigned int,
        const float,
        const float
      );


      void sigmaFilter(
        std::vector<std::vector<short>>&,
        const std::vector<std::vector<short>>&,
        const float,
        const unsigned int,
        const unsigned int,
        const unsigned int,
        const float);

      void sigmaFilter(
        std::vector<std::vector<float>>&,
        const std::vector<std::vector<float>>&,
        const float,
        const unsigned int,
        const unsigned int,
        const unsigned int,
        const float);

      void sigmaFilter(
        std::vector<std::vector<double>>&,
        const std::vector<std::vector<double>>&,
        const float,
        const unsigned int,
        const unsigned int,
        const unsigned int,
        const float);

      /// Default destructor
      ~AdaptiveWiener(){}

    private:

      template <typename T>
      void filterLee(
        std::vector<std::vector<T> >& deconvolvedWaveform,
        const std::vector<std::vector<T> >& waveLessCoherent,
        const float noiseVar,
        const unsigned int sx=7,
        const unsigned int sy=7);

      template <typename T>
      void MMWF(
        std::vector<std::vector<T> >& deconvolvedWaveform,
        const std::vector<std::vector<T> >& waveLessCoherent,
        const float noiseVar,
        const unsigned int sx=7,
        const unsigned int sy=7);

      template <typename T>
      void MMWFStar(
        std::vector<std::vector<T> >& deconvolvedWaveform,
        const std::vector<std::vector<T> >& waveLessCoherent,
        const unsigned int sx=7,
        const unsigned int sy=7);

      template <typename T>
      void filterLeeEnhanced(
        std::vector<std::vector<T> >& deconvolvedWaveform,
        const std::vector<std::vector<T> >& waveLessCoherent,
        const float noiseVar,
        const unsigned int sx=3,
        const unsigned int sy=3,
        const float a=1,
        const float epsilon=2.5);

      template <typename T>
      void adaptiveROIWiener(
        std::vector<std::vector<T>>& deconvolvedWaveform,
        const std::vector<std::vector<T>>& waveLessCoherent,
        const std::vector<std::vector<bool>>& selectVals,
        const float noiseVar,
        const unsigned int sx=3,
        const unsigned int sy=3,
        const float a=1,
        const float epsilon=2.5);

      template <typename T>
      void sigmaFilter(
        std::vector<std::vector<T>>& deconvolvedWaveform,
        const std::vector<std::vector<T> >& waveLessCoherent,
        const float noiseVar,
        const unsigned int sx=7,
        const unsigned int sy=7,
        const unsigned int K=5,
        const float sigmaFactor=2.0);
  };
}

#endif
/** @} */ // end of doxygen group

/**
 * \file MorphologicalCNC.h
 *
 * \ingroup sigproc_tools
 *
 * \brief Class def header for a class MorphologicalCNC
 *
 * @author koh0207
 */

/** \addtogroup sigproc_tools

    @{*/
#ifndef __SIGPROC_TOOLS_MORPHOLOGICALCNC_H__
#define __SIGPROC_TOOLS_MORPHOLOGICALCNC_H__

#include <vector>
#include <iostream>
#include <string>
#include <algorithm>
#include <iterator>
#include <numeric>
#include <cmath>
#include <functional>
#include "Morph1DFast.h"
#include "Morph2DFast.h"
#include "MiscUtils.h"

namespace sigproc_tools {

  /**
     \class MorphologicalCNC
     User defined class MorphologicalCNC ... these comments are used to generate
     doxygen documentation!
  */
  class MorphologicalCNC{

    public:

      /// Default constructor
      MorphologicalCNC(){}

      void getSelectVals(
        const std::vector<std::vector<short>>&,
        std::vector<std::vector<bool>>&,
        std::vector<std::vector<bool>>&,
        const unsigned int,
        const float);

      void getSelectVals(
        const std::vector<std::vector<float>>&,
        std::vector<std::vector<bool>>&,
        std::vector<std::vector<bool>>&,
        const unsigned int,
        const float);

      void getSelectVals(
        const std::vector<std::vector<double>>&,
        std::vector<std::vector<bool>>&,
        std::vector<std::vector<bool>>&,
        const unsigned int,
        const float);

      void denoiseCoherent1D(
        std::vector<std::vector<short>>&,
        const std::vector<std::vector<short>>&,
        std::vector<std::vector<bool>>&,
        std::vector<std::vector<bool>>&,
        const char,
        const unsigned int,
        const unsigned int,
        const unsigned int,
        const float);

      void denoiseCoherent1D(
        std::vector<std::vector<float>>&,
        const std::vector<std::vector<float>>&,
        std::vector<std::vector<bool>>&,
        std::vector<std::vector<bool>>&,
        const char,
        const unsigned int,
        const unsigned int,
        const unsigned int,
        const float);

      void denoiseCoherent1D(
        std::vector<std::vector<double>>&,
        const std::vector<std::vector<double>>&,
        std::vector<std::vector<bool>>&,
        std::vector<std::vector<bool>>&,
        const char,
        const unsigned int,
        const unsigned int,
        const unsigned int,
        const float);


      void denoiseCoherent2D(
        std::vector<std::vector<short>>&,
        const std::vector<std::vector<short>>&,
        std::vector<std::vector<bool>>&,
        std::vector<std::vector<bool>>&,
        const char,
        const unsigned int,
        const unsigned int,
        const unsigned int,
        const unsigned int,
        const float);

      void denoiseCoherent2D(
        std::vector<std::vector<float>>&,
        const std::vector<std::vector<float>>&,
        std::vector<std::vector<bool>>&,
        std::vector<std::vector<bool>>&,
        const char,
        const unsigned int,
        const unsigned int,
        const unsigned int,
        const unsigned int,
        const float);

      void denoiseCoherent2D(
        std::vector<std::vector<double>>&,
        const std::vector<std::vector<double>>&,
        std::vector<std::vector<bool>>&,
        std::vector<std::vector<bool>>&,
        const char,
        const unsigned int,
        const unsigned int,
        const unsigned int,
        const unsigned int,
        const float);

    private:

      template <typename T>
      void getSelectVals(
        const std::vector<std::vector<T> >& morphedWaveforms,
        std::vector<std::vector<bool> >& selectVals,
        std::vector<std::vector<bool> >& roi,
        const unsigned int window,
        const float thresholdFactor);

      template <typename T>
      void denoiseCoherent1D(
        std::vector<std::vector<T> >& waveLessCoherent,
        const std::vector<std::vector<T> >& fullEvent,
        std::vector<std::vector<bool> >& selectVals,
        std::vector<std::vector<bool> >& roi,
        const char filterName='g',
        const unsigned int grouping=64,
        const unsigned int structuringElement=5,
        const unsigned int window=0,
        const float thresholdFactor=2.5);

      template <typename T>
      void denoiseCoherent2D(
        std::vector<std::vector<T> >& waveLessCoherent,
        const std::vector<std::vector<T> >& fullEvent,
        std::vector<std::vector<bool> >& selectVals,
        std::vector<std::vector<bool> >& roi,
        const char filterName='g',
        const unsigned int grouping=64,
        const unsigned int structuringElementx=5,
        const unsigned int structuringElementy=20,
        const unsigned int window=0,
        const float thresholdFactor=2.5);
      /// Default destructor
      ~MorphologicalCNC(){}

  };
}

#endif
/** @} */ // end of doxygen group

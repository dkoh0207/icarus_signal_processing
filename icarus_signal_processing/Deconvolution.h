/**
 * \file Deconvolution.h
 *
 * \ingroup sigproc_tools
 * 
 * \brief Class def header for a class Deconvolution
 *
 * @author koh0207
 */

/** \addtogroup sigproc_tools

    @{*/
#ifndef __SIGPROC_TOOLS_DECONVOLUTION_H__

#include <vector>
#include <iostream>
#include <algorithm>
#include <iterator>
#include <numeric>
#include <cmath>
#include <functional>
#include "MiscUtils.h"

#include <Eigen/Core>
#include <unsupported/Eigen/FFT>


#define __SIGPROC_TOOLS_DECONVOLUTION_H__

namespace sigproc_tools {

  /**
     \class Deconvolution
     User defined class Deconvolution ... these comments are used to generate
     doxygen documentation!
  */
  class Deconvolution{
    
    public:
      
      /// Default constructor
      Deconvolution(){}

      void Inverse1D(
        std::vector<std::vector<float>>&,
        const std::vector<std::vector<float>>&,
        const std::vector<float>&
      );

      void Inverse1D(
        std::vector<std::vector<double>>&,
        const std::vector<std::vector<double>>&,
        const std::vector<double>&
      );


      void Wiener1D(
        std::vector<std::vector<float>>&,
        const std::vector<std::vector<float>>&,
        const std::vector<float>&,
        const float
      );

      void Wiener1D(
        std::vector<std::vector<double>>&,
        const std::vector<std::vector<double>>&,
        const std::vector<double>&,
        const float
      );

      
      /// Default destructor
      ~Deconvolution(){}

    private:

      template <typename T>
      void Inverse1D(
        std::vector<std::vector<T>>& outputWaveform,
        const std::vector<std::vector<T>>& inputWaveform,
        const std::vector<T>& responseFunction
      );

      template <typename T>
      void Wiener1D(
        std::vector<std::vector<T>>& outputWaveform,
        const std::vector<std::vector<T>>& inputWaveform,
        const std::vector<T>& responseFunction,
        const float noiseVar
      );
      
    };
}

#endif
/** @} */ // end of doxygen group 


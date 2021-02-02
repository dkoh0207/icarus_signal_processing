/**
 * \file BilateralFilters.h
 *
 * \ingroup sigproc_tools
 * 
 * \brief Class def header for a class BilateralFilters
 *
 * @author koh0207
 */

/** \addtogroup sigproc_tools

    @{*/
#ifndef __SIGPROC_TOOLS_BILATERALFILTERS_H__
#define __SIGPROC_TOOLS_BILATERALFILTERS_H__

#include <vector>
#include <iostream>
#include <string>
#include <algorithm>
#include <iterator>
#include <numeric>
#include <cmath>
#include <functional>
#include <complex>

namespace icarus_signal_processing {

  /**
     \class BilateralFilters
     User defined class BilateralFilters ... these comments are used to generate
     doxygen documentation!
  */

  template <class T> using Array2D = std::vector<std::vector<T>>;
  template <class T> using Array3D = std::vector<Array2D<T>>;
  // std::complex<float> I(0, 1);

  class BilateralFilters{
    
    public:
      
      /// Default constructor
      BilateralFilters(){}

      void bilateral(
        const Array2D<float>& input2D,
        Array2D<float>& output2D,
        const unsigned int sx,
        const unsigned int sy,
        float sigma_x,
        float sigma_y,
        float sigma_r
      ) const;

      void directional(
        const Array2D<float>& input2D,
        const Array2D<float>& direction,
        Array2D<float>& output2D,
        const unsigned int sx,
        const unsigned int sy,
        float sigma_x,
        float sigma_y,
        float sigma_theta,
        const int angleMesh
      ) const;

      void Gaussian(
        const Array2D<float>& input2D,
        const float sigma_x,
        const float sigma_y,
        Array2D<float>& output2D
      ) const;

      void Gaussian(
        const Array2D<std::complex<float>>& input2D,
        const float sigma_x,
        const float sigma_y,
        Array2D<std::complex<float>>& output2D
      ) const;

      void Gaussian1D(
        const std::vector<float>& input1D,
        const std::vector<float>& filterCoeffs,
        std::vector<float>& output1D) const;

      void Gaussian1D(
        const std::vector<std::complex<float>>& input1D,
        const std::vector<float>& filterCoeffs,
        std::vector<std::complex<float>>& output1D) const;

      float Discretize(
        const Array2D<float>& input2D,
        const unsigned int maxPixelValue,
        Array2D<short>& output2D
      ) const;

      float Discretize(
        const Array2D<float>& input2D,
        const unsigned int maxPixelValue,
        Array2D<float>& output2D
      ) const;

      void CSUBilateral(
        const Array2D<float>& input2D,
        Array2D<float>& output2D,
        const unsigned int T,
        const float sigma_s,
        const float sigma_r
      ) const;

      void GPABilateral(
        const Array2D<float>& input2D,
        Array2D<float>& output2D,
        const unsigned int T,
        const float sigma_x,
        const float sigma_y,
        const float sigma_r,
        const float epsilon
      ) const;

      unsigned int estimateGPAOrder(
        float sigma_r,
        float epsilon,
        float T
      ) const;

      std::vector<float> computeGausCoeffs(const float sigma) const;

      void guidedBilateral(
        const Array2D<float>& input2D,
        Array2D<float>& output2D,
        const unsigned int sx,
        const unsigned int sy
      ) const;
      
      /// Default destructor
      ~BilateralFilters(){}

    private:

      template <typename T>
      void Gaussian1D(
        const std::vector<T>& input1D,
        const std::vector<float>& filterCoeffs,
        std::vector<T>& output1D) const;

      template <typename T>
      void Gaussian(
        const Array2D<T>& input2D,
        const float sigma_x,
        const float sigma_y,
        Array2D<T>& output2D
      ) const;

      template <typename T>
      float Discretize(
        const Array2D<float>& input2D,
        const unsigned int maxPixelValue,
        Array2D<T>& output2D
      ) const;
      
      
    };
}

#endif
/** @} */ // end of doxygen group 


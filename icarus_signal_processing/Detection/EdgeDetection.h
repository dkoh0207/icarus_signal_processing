/**
 * \file EdgeDetection.h
 *
 * \ingroup icarus_signal_processing
 * 
 * \brief Class def header for a class EdgeDetection
 *
 * @author koh0207
 */

/** \addtogroup icarus_signal_processing

    @{*/
#ifndef __SIGPROC_TOOLS_EDGEDETECTION_H__
#define __SIGPROC_TOOLS_EDGEDETECTION_H__
#define _USE_MATH_DEFINES

#include <vector>
#include <iostream>
#include <string>
#include <algorithm>
#include <iterator>
#include <cmath>
#include <assert.h>
#include <queue>

namespace icarus_signal_processing
{

  /**
     \class EdgeDetection
     User defined class EdgeDetection ... these comments are used to generate
     doxygen documentation!
  */

  template <class T>
  using Array2D = std::vector<std::vector<T>>;

  class EdgeDetection
  {

  public:
    /// Default constructor
    EdgeDetection() {}

    void Convolve2D(
        const Array2D<float> &input2D,
        Array2D<float> &output2D,
        const Array2D<float> &kernel) const;

    void SobelX(
        const Array2D<float> &input2D,
        Array2D<float> &gradient) const;

    void SobelY(
        const Array2D<float> &input2D,
        Array2D<float> &gradient) const;

    void Sobel(
        const Array2D<float> &input2D,
        Array2D<float> &sobelX,
        Array2D<float> &sobelY,
        Array2D<float> &gradient,
        Array2D<float> &direction) const;

    void LSDGradX(
        const Array2D<float> &input2D,
        Array2D<float> &gradient) const;

    void LSDGradY(
        const Array2D<float> &input2D,
        Array2D<float> &gradient) const;

    void LSDGrad(
        const Array2D<float> &input2D,
        Array2D<float> &gradX,
        Array2D<float> &gradY,
        Array2D<float> &gradient,
        Array2D<float> &direction) const;

    void Gauss(
        const Array2D<float> &input2D,
        Array2D<float> &output2D,
        unsigned int radius) const;

    void EdgeNMS(
        const Array2D<float> &gradient2D,
        const Array2D<float> &radians2D,
        Array2D<float> &output2D) const;

    // TODO: NEEDS DEBUGGING
    void EdgeNMSInterpolation(
        const Array2D<float> &gradient2D,
        const Array2D<float> &gradX,
        const Array2D<float> &gradY,
        const Array2D<float> &degrees2D,
        Array2D<float> &output2D) const;

    void DoubleThresholding(
        const Array2D<float> &doneNMS2D,
        Array2D<bool> &binary2D,
        std::vector<int> &strongEdgeRows,
        std::vector<int> &strongEdgeCols,
        std::vector<int> &weakEdgeRows,
        std::vector<int> &weakEdgeCols,
        float lowThreshold,
        float highThreshold) const;

    void HysteresisThresholding(
        const Array2D<bool> &binary2D,
        const std::vector<int> &strongEdgeRows,
        const std::vector<int> &strongEdgeCols,
        const std::vector<int> &weakEdgeRows,
        const std::vector<int> &weakEdgeCols,
        Array2D<bool> &output2D) const;

    void Canny(
        const Array2D<float> &waveLessCoherent,
        Array2D<bool> &output2D,
        const unsigned int sx,
        const unsigned int sy,
        const float sigma_x,
        const float sigma_y,
        const float sigma_r,
        const float lowThreshold,
        const float highThreshold,
        const char mode) const;

    void gradientRegionGrow(
        const Array2D<float> &direction,
        const int anchorX,
        const int anchorY,
        const int regionID,
        const float tolerance,
        const unsigned int windowX,
        const unsigned int windowY,
        Array2D<int> &partitions) const;

    void regionGrow2D(
        const Array2D<float> &direction,
        const std::vector<int> &anchorsX,
        const std::vector<int> &anchorsY,
        const float tolerance,
        const unsigned int windowX,
        const unsigned int windowY,
        Array2D<int> &partitions) const;

    /// Default destructor
    ~EdgeDetection() {}

  private:
    void getDilation1D(const std::vector<bool>&,
                       const unsigned int,
                       std::vector<bool>&) const;

    void getDilation2D(const std::vector<std::vector<bool>>&,
                       const unsigned int,
                       const unsigned int,
                       std::vector<std::vector<bool>>&) const;
  };

} // namespace icarus_signal_processing

#endif
/** @} */ // end of doxygen group

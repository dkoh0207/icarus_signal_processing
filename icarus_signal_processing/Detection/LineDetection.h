/**
 * \file LineDetection.h
 *
 * \ingroup icarus_signal_processing
 * 
 * \brief Class def header for a class LineDetection
 *
 * @author koh0207
 */

/** \addtogroup icarus_signal_processing

    @{*/
#ifndef __icarus_signal_processing_LINEDETECTION_H__
#define __icarus_signal_processing_LINEDETECTION_H__

#include <vector>
#include <cmath>
#include <algorithm>
#include <iostream>
#include <functional>

namespace icarus_signal_processing {

  /**
     \class LineDetection
     User defined class LineDetection ... these comments are used to generate
     doxygen documentation!
  */

  template <class T> using Array2D = std::vector<std::vector<T>>;

  struct Point {
    int x;
    int y;
  };

  struct Line {
    int x0;
    int y0;
    float slope;
  };

  class LineDetection{
    
    public:
      
      /// Default constructor
      LineDetection(){}
      
      /// Default destructor
      ~LineDetection(){}
      
      void HoughTransform(const Array2D<bool>& binary2D,
                          Array2D<int>& accumulator2D,
                          const unsigned int thetaSteps = 360) const;

      
      void FindPeaksNMS(const Array2D<int>& accumulator2D,
                        std::vector<int>& rhoIndex,
                        std::vector<int>& thetaIndex,
                        const unsigned int threshold = 100,
                        const unsigned int angleWindow = 10,
                        const unsigned int maxLines = 20,
                        const unsigned int windowSize = 2) const;

      Line getLine(const Array2D<bool>& binary2D,
                   const float theta,
                   const float r,
                   const float eps = 0.001,
                   const unsigned int angleWindow = 10) const;

      void drawLine(Array2D<bool>& newSelectVals,
                    const Line line,
                    const unsigned int dilationX = 7,
                    const unsigned int dilationY = 20) const;

      void refineSelectVals(const Array2D<bool>& selectVals,
                            Array2D<bool>& refinedSelectVals,
                            const size_t thetaSteps = 360,
                            const unsigned int threshold = 100,
                            const unsigned int angleWindow = 10,
                            const unsigned int maxLines = 20,
                            const unsigned int windowSize = 10,
                            const unsigned int dilationX = 7,
                            const unsigned int dilationY = 20,
                            const float eps = 0.00001) const;

      void FindVerticalSegments(const Array2D<bool>& selectVals,
                                Array2D<bool>& newMask2D,
                                const size_t thetaSteps,
                                const std::vector<int>& rhoIndex,
                                const std::vector<int>& thetaIndex,
                                const unsigned int maxGap = 10,
                                const unsigned int minLineLength = 100,
                                const unsigned int margin = 3) const;

      
      // void WindowedHoughTransform(const Array2D<bool>& binary2D,
      //                             Array2D<int>& accumulator2D,
      //                             const size_t mSteps = 360,
      //                             const size_t thetaSteps = 360) const;

      
      // void FindSegmentInWindow(const Array2D<bool>& binary2D,
      //                          const size_t centerx, 
      //                          const size_t centery,
      //                          unsigned int& r,
      //                          float& theta,
      //                          const size_t windowSize = 3) const;
    
  };
}

#endif
/** @} */ // end of doxygen group 


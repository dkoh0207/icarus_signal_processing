/**
 * \file MiscUtils.h
 *
 * \ingroup sigproc_tools
 * 
 * \brief Class def header for a class MiscUtils
 *
 * @author koh0207
 */

/** \addtogroup sigproc_tools

    @{*/
#ifndef __SIGPROC_TOOLS_MISCUTILS_H__
#define __SIGPROC_TOOLS_MISCUTILS_H__

#include <vector>
#include <algorithm>
#include <iterator>
#include <numeric>
#include <cmath>
#include <functional>

namespace sigproc_tools {

  /**
     \class MiscUtils
     Miscellaneous utility functions. 
  */
  class MiscUtils{
    
    public:
      
      /// Default constructor
      MiscUtils(){}

      short computeMedian(const std::vector<short>& vec);
      float computeMedian(const std::vector<float>& vec);
      double computeMedian(const std::vector<double>& vec);

      float estimateNoiseVariance(
        const std::vector<std::vector<float>>& waveLessCoherent,
        const std::vector<std::vector<bool>>& selectVals);

      float estimateMAD(
        const std::vector<float>& wf);

      
      /// Default destructor
      ~MiscUtils(){}
    
    private:

      template <typename T>
      T computeMedian(const std::vector<T>& vec);


      // template <typename T> T computeMedian(
      //   const std::vector<T>& waveform
      // );
  };
}

#endif
/** @} */ // end of doxygen group 


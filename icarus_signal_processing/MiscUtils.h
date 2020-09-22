/**
 * \file MiscUtils.h
 *
 * \ingroup icarus_signal_processing
 * 
 * \brief Class def header for a class MiscUtils
 *
 * @author koh0207
 */

/** \addtogroup icarus_signal_processing

    @{*/
#ifndef __icarus_signal_processing_MISCUTILS_H__
#define __icarus_signal_processing_MISCUTILS_H__

#include <vector>
#include <algorithm>
#include <iterator>
#include <numeric>
#include <cmath>
#include <functional>

namespace icarus_signal_processing {

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


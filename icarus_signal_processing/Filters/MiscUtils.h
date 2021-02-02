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

      short  computeMedian(const std::vector<short>&);
      float  computeMedian(const std::vector<float>&);
      double computeMedian(const std::vector<double>&);

      short  computeMaximum(const std::vector<std::vector<short>>&);
      float  computeMaximum(const std::vector<std::vector<float>>&);
      double computeMaximum(const std::vector<std::vector<double>>&);

      short  computeMinimum(const std::vector<std::vector<short>>&);
      float  computeMinimum(const std::vector<std::vector<float>>&);
      double computeMinimum(const std::vector<std::vector<double>>&);

      float estimateNoiseVariance(
        const std::vector<std::vector<float>>& waveLessCoherent,
        const std::vector<std::vector<bool>>& selectVals);

      float estimateMAD(
        const std::vector<float>& wf);

      unsigned nChoosek(unsigned n, unsigned k) const;

      void drawIndextoImage(const std::vector<int>&,
                            const std::vector<int>&,
                            std::vector<std::vector<bool>>&) const;

      void drawIndextoImage(const std::vector<size_t>&,
                            const std::vector<size_t>&,
                            std::vector<std::vector<bool>>&) const;

      void drawIndextoImage(const std::vector<unsigned int>&,
                            const std::vector<unsigned int>&,
                            std::vector<std::vector<bool>>&) const;

      void drawIndextoImage(const std::vector<short>&,
                            const std::vector<short>&,
                            std::vector<std::vector<bool>>&) const;

      
      /// Default destructor
      ~MiscUtils(){}
    
    private:

      template <typename T> 
      T computeMedian(const std::vector<T>&);

      template <typename T>
      T computeMaximum(const std::vector<std::vector<T>>&);

      template <typename T>
      T computeMinimum(const  std::vector<std::vector<T>>&);

      template <typename T>
      void drawIndextoImage(const std::vector<T>&,
                            const std::vector<T>&,
                            std::vector<std::vector<bool>>&) const;

      // template <typename T> T computeMedian(
      //   const std::vector<T>& waveform
      // );
  };
}

#endif
/** @} */ // end of doxygen group 


/**
 * \file Morph2D.h
 *
 * \ingroup icarus_signal_processing
 * 
 * \brief Class def header for a class Morph2D
 *
 * @author koh0207
 */

/** \addtogroup icarus_signal_processing

    @{*/
#ifndef __icarus_signal_processing_MORPH2D_H__
#define __icarus_signal_processing_MORPH2D_H__

#include <iostream>
#include <vector>
#include <algorithm>
#include <iterator>
#include <numeric>
#include <cmath>
#include <functional>

namespace icarus_signal_processing {

  /**
     \class Morph2D
     2D Morphological Filters
  */
  class Morph2D{
    
  public:
    
    /// Default constructor
    Morph2D(){}

      void getFilter2D(const std::vector<std::vector<short> >&,
                      const unsigned int,
                      const unsigned int,
                      std::vector<std::vector<short> >&,
                      std::vector<std::vector<short> >&,
                      std::vector<std::vector<short> >&,
                      std::vector<std::vector<short> >&) const;

      void getFilter2D(const std::vector<std::vector<float> >&,
                      const unsigned int,
                      const unsigned int,
                      std::vector<std::vector<float> >&,
                      std::vector<std::vector<float> >&,
                      std::vector<std::vector<float> >&,
                      std::vector<std::vector<float> >&) const;

      void getFilter2D(const std::vector<std::vector<double> >&,
                      const unsigned int,
                      const unsigned int,
                      std::vector<std::vector<double> >&,
                      std::vector<std::vector<double> >&,
                      std::vector<std::vector<double> >&,
                      std::vector<std::vector<double> >&) const;


      void getDilation(const std::vector<std::vector<short> >&,
                      const unsigned int,
                      const unsigned int,
                      std::vector<std::vector<short> >&) const;

      void getDilation(const std::vector<std::vector<float> >&,
                      const unsigned int,
                      const unsigned int,
                      std::vector<std::vector<float> >&) const;

      void getDilation(const std::vector<std::vector<double> >&,
                      const unsigned int,
                      const unsigned int,
                      std::vector<std::vector<double> >&) const;


      void getErosion(const std::vector<std::vector<short> >&,
                      const unsigned int,
                      const unsigned int,
                      std::vector<std::vector<short> >&) const;

      void getErosion(const std::vector<std::vector<float> >&,
                      const unsigned int,
                      const unsigned int,
                      std::vector<std::vector<float> >&) const;

      void getErosion(const std::vector<std::vector<double> >&,
                      const unsigned int,
                      const unsigned int,
                      std::vector<std::vector<double> >&) const;


      void getGradient(const std::vector<std::vector<short> >&,
                      const unsigned int,
                      const unsigned int,
                      std::vector<std::vector<short> >&) const;

      void getGradient(const std::vector<std::vector<float> >&,
                      const unsigned int,
                      const unsigned int,
                      std::vector<std::vector<float> >&) const;

      void getGradient(const std::vector<std::vector<double> >&,
                      const unsigned int,
                      const unsigned int,
                      std::vector<std::vector<double> >&) const;


      void getMedian(const std::vector<std::vector<short> >&,
                      const unsigned int,
                      const unsigned int,
                      std::vector<std::vector<short> >&) const;

      void getMedian(const std::vector<std::vector<float> >&,
                      const unsigned int,
                      const unsigned int,
                      std::vector<std::vector<float> >&) const;

      void getMedian(const std::vector<std::vector<double> >&,
                      const unsigned int,
                      const unsigned int,
                      std::vector<std::vector<double> >&) const;


      void getOpeningAndClosing(const std::vector<std::vector<short> >&,
                      const unsigned int,
                      const unsigned int,
                      std::vector<std::vector<short> >&,
                      std::vector<std::vector<short> >&) const;

      void getOpeningAndClosing(const std::vector<std::vector<float> >&,
                      const unsigned int,
                      const unsigned int,
                      std::vector<std::vector<float> >&,
                      std::vector<std::vector<float> >&) const;

      void getOpeningAndClosing(const std::vector<std::vector<double> >&,
                      const unsigned int,
                      const unsigned int,
                      std::vector<std::vector<double> >&,
                      std::vector<std::vector<double> >&) const;

      /// Default destructor
      ~Morph2D(){}
      
    private:

      template <typename T> 
      void getFilter2D(
        const std::vector<std::vector<T> >& waveform2D,
        const unsigned int structuringElementx,
        const unsigned int structuringElementy,
        std::vector<std::vector<T> >& dilation2D,
        std::vector<std::vector<T> >& erosion2D,
        std::vector<std::vector<T> >& average2D,
        std::vector<std::vector<T> >& gradient2D) const;

      template <typename T> 
      void getGradient(
        const std::vector<std::vector<T> >& waveform2D,
        const unsigned int structuringElementx,
        const unsigned int structuringElementy,
        std::vector<std::vector<T> >& gradient2D) const;

      template <typename T> 
      void getDilation(
        const std::vector<std::vector<T> >& waveform2D,
        const unsigned int structuringElementx,
        const unsigned int structuringElementy,
        std::vector<std::vector<T> >& dilation2D) const;

      template <typename T> 
      void getErosion(
        const std::vector<std::vector<T> >& waveform2D,
        const unsigned int structuringElementx,
        const unsigned int structuringElementy,
        std::vector<std::vector<T> >& erosion2D) const;

      template <typename T> 
      void getMedian(
        const std::vector<std::vector<T> >& waveform2D,
        const unsigned int structuringElementx,
        const unsigned int structuringElementy,
        std::vector<std::vector<T> >& median2D) const;

      template <typename T> 
      void getOpeningAndClosing(
        const std::vector<std::vector<T> >& waveform2D,
        const unsigned int structuringElementx,
        const unsigned int structuringElementy,
        std::vector<std::vector<T> >& opening2D,
        std::vector<std::vector<T> >& closing2D) const;
    
  };
}

#endif
/** @} */ // end of doxygen group 


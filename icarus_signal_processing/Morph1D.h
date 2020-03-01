/**
 * \file Morph1D.h
 *
 * \ingroup icarus_signal_processing
 * 
 * \brief Class def header for a class Morph1D
 *
 * @author koh0207
 */

/** \addtogroup icarus_signal_processing

    @{*/
#ifndef __SIGPROC_TOOLS_MORPH1D_H__
#define __SIGPROC_TOOLS_MORPH1D_H__

#include <iostream>
#include <vector>
#include <algorithm>
#include <iterator>
#include <numeric>
#include <cmath>
#include <functional>

namespace icarus_signal_processing {

  /**
     \class Morph1D
     Parent class for 1D Morphological operations for signal processing.
  */

  template <class T> using Waveform = std::vector<T>;

  class Morph1D{
    
    public:
      
      /// Default constructor
      Morph1D(){}

      void getWaveformParams(const std::vector<short>&,
                             float&,
                             float&,
                             float&);

      void getWaveformParams(const std::vector<float>&,
                             float&,
                             float&,
                             float&);

      void getWaveformParams(const std::vector<double>&,
                             float&,
                             float&,
                             float&);


      void getDilation(const Waveform<short>&,
                      const unsigned int,
                      Waveform<short>&) const;

      void getDilation(const Waveform<float>&,
                      const unsigned int,
                      Waveform<float>&) const;

      void getDilation(const Waveform<double>&,
                      const unsigned int,
                      Waveform<double>&) const;


      void getErosion(const Waveform<short>&,
                      const unsigned int,
                      Waveform<short>&) const;

      void getErosion(const Waveform<float>&,
                      const unsigned int,
                      Waveform<float>&) const;

      void getErosion(const Waveform<double>&,
                      const unsigned int,
                      Waveform<double>&) const;


      void getGradient(const Waveform<short>&,
                      const unsigned int,
                      Waveform<short>&) const;

      void getGradient(const Waveform<float>&,
                      const unsigned int,
                      Waveform<float>&) const;

      void getGradient(const Waveform<double>&,
                      const unsigned int,
                      Waveform<double>&) const;


      void getAverage(const Waveform<short>&,
                      const unsigned int,
                      Waveform<short>&) const;

      void getAverage(const Waveform<float>&,
                      const unsigned int,
                      Waveform<float>&) const;

      void getAverage(const Waveform<double>&,
                      const unsigned int,
                      Waveform<double>&) const;


      void getMedian(const Waveform<short>&,
                      const unsigned int,
                      Waveform<short>&) const;

      void getMedian(const Waveform<float>&,
                      const unsigned int,
                      Waveform<float>&) const;

      void getMedian(const Waveform<double>&,
                      const unsigned int,
                      Waveform<double>&) const;


      void getOpeningAndClosing(const Waveform<short>&,
                      const unsigned int,
                      Waveform<short>&,
                      Waveform<short>&) const;

      void getOpeningAndClosing(const Waveform<float>&,
                      const unsigned int,
                      Waveform<float>&,
                      Waveform<float>&) const;

      void getOpeningAndClosing(const Waveform<double>&,
                      const unsigned int,
                      Waveform<double>&,
                      Waveform<double>&) const;

      /// Default destructor
      ~Morph1D(){}
      
    private:

      template <typename T> 
      void getWaveformParams(
        const std::vector<T>& waveform,
        float& mean,
        float& median,
        float& rms);

      template <typename T> 
      void getDilation(
        const Waveform<T>& inputWaveform,
        const unsigned int structuringElement,
        Waveform<T>& dilationVec) const;

      template <typename T> 
      void getErosion(
        const Waveform<T>& inputWaveform,
        const unsigned int structuringElement,
        Waveform<T>& erosionVec) const;

      template <typename T> 
      void getGradient(
        const Waveform<T>& inputWaveform,
        const unsigned int structuringElement,
        Waveform<T>& gradientVec) const;
    
      template <typename T> 
      void getAverage(
        const Waveform<T>& inputWaveform,
        const unsigned int structuringElement,
        Waveform<T>& averageVec) const;

      template <typename T> 
      void getMedian(
        const Waveform<T>& inputWaveform,
        const unsigned int structuringElement,
        Waveform<T>& medianVec) const;

      template <typename T> 
      void getOpeningAndClosing(
        const Waveform<T>& inputWaveform,
        const unsigned int structuringElement,
        Waveform<T>& openingVec,
        Waveform<T>& closingVec) const;
  };
}

#endif
/** @} */ // end of doxygen group 


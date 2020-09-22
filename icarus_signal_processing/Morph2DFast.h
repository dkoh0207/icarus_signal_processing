/**
 * \file Morph2DFast.h
 *
 * \ingroup icarus_signal_processing
 *
 * \brief Class def header for a class Morph2DFast
 *
 * @author koh0207
 */

/** \addtogroup icarus_signal_processing

    @{*/
#ifndef __icarus_signal_processing_MORPH2DFAST_H__
#define __icarus_signal_processing_MORPH2DFAST_H__

#include <iostream>
#include <vector>
#include <algorithm>
#include <iterator>
#include <numeric>
#include <cmath>
#include <functional>
#include "Morph1DFast.h"

namespace icarus_signal_processing {

  /**
     \class Morph2DFast
     User defined class Morph2DFast ... these comments are used to generate
     doxygen documentation!
  */
  class Morph2DFast{

    public:

      /// Default constructor
      Morph2DFast(){}

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


        void getOpening(const std::vector<std::vector<short> >&,
                        const unsigned int,
                        const unsigned int,
                        std::vector<std::vector<short> >&) const;

        void getOpening(const std::vector<std::vector<float> >&,
                        const unsigned int,
                        const unsigned int,
                        std::vector<std::vector<float> >&) const;

        void getOpening(const std::vector<std::vector<double> >&,
                        const unsigned int,
                        const unsigned int,
                        std::vector<std::vector<double> >&) const;

        void getClosing(const std::vector<std::vector<short> >&,
                        const unsigned int,
                        const unsigned int,
                        std::vector<std::vector<short> >&) const;

        void getClosing(const std::vector<std::vector<float> >&,
                        const unsigned int,
                        const unsigned int,
                        std::vector<std::vector<float> >&) const;

        void getClosing(const std::vector<std::vector<double> >&,
                        const unsigned int,
                        const unsigned int,
                        std::vector<std::vector<double> >&) const;
      /// Default destructor
      ~Morph2DFast(){}


    private:

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
      void getGradient(
        const std::vector<std::vector<T> >& waveform2D,
        const unsigned int structuringElementx,
        const unsigned int structuringElementy,
        std::vector<std::vector<T> >& gradient2D) const;

      template <typename T>
      void getMedian(
        const std::vector<std::vector<T> >& waveform2D,
        const unsigned int structuringElementx,
        const unsigned int structuringElementy,
        std::vector<std::vector<T> >& median2D) const;

      template <typename T>
      void getClosing(
        const std::vector<std::vector<T> >& waveform2D,
        const unsigned int structuringElementx,
        const unsigned int structuringElementy,
        std::vector<std::vector<T> >& closing2D) const;

      template <typename T>
      void getOpening(
        const std::vector<std::vector<T> >& waveform2D,
        const unsigned int structuringElementx,
        const unsigned int structuringElementy,
        std::vector<std::vector<T> >& opening2D) const;

  };
}

#endif
/** @} */ // end of doxygen group

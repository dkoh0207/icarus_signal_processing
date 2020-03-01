/**
 * \file Denoising.h
 *
 * \ingroup icarus_signal_processing
 * 
 * \brief Class def header for a class Denoising
 *
 * @author koh0207
 */

/** \addtogroup icarus_signal_processing

    @{*/
#ifndef __SIGPROC_TOOLS_DENOISING_H__
#define __SIGPROC_TOOLS_DENOISING_H__

#include <vector>
#include <iostream>
#include <string>
#include <algorithm>
#include <iterator>
#include <numeric>
#include <cmath>
#include <functional>
#include "Morph1D.h"
#include "Morph2D.h"
#include "ICARUSSigProcDefs.h"

namespace icarus_signal_processing {

/**
   \class Denoising
   User defined class Denoising ... these comments are used to generate
   doxygen documentation!
*/
class Denoising {
  
public:
  
  /// Default constructor
  Denoising(){}

  void getSelectVals(const ArrayShort&,
                     const ArrayShort&,
                     ArrayBool&,
                     ArrayBool&,
                     const unsigned int,
                     const float);

  void getSelectVals(
    const ArrayFloat&,
    const ArrayFloat&,
    ArrayBool&,
    ArrayBool&,
    const unsigned int,
    const float);

  void getSelectVals(
    const ArrayDouble&,
    const ArrayDouble&,
    ArrayBool&,
    ArrayBool&,
    const unsigned int,
    const float);

  void removeCoherentNoise1D(
    ArrayShort&,
    const ArrayShort&,
    ArrayShort&,
    ArrayShort&,
    ArrayBool&,
    ArrayBool&,
    ArrayShort&,
    const char,
    const unsigned int,
    const unsigned int,
    const unsigned int,
    const float 
  );

  void removeCoherentNoise1D(
    ArrayFloat&,
    const ArrayFloat&,
    ArrayFloat&,
    ArrayFloat&,
    ArrayBool&,
    ArrayBool&,
    ArrayFloat&,
    const char,
    const unsigned int,
    const unsigned int,
    const unsigned int,
    const float 
  );

  void removeCoherentNoise1D(
    ArrayDouble&, 
    const ArrayDouble&,
    ArrayDouble&,
    ArrayDouble&,
    ArrayBool&,
    ArrayBool&,
    ArrayDouble&,
    const char,
    const unsigned int,
    const unsigned int,
    const unsigned int,
    const float 
  );

  void removeCoherentNoise2D(
    ArrayShort&,
    const ArrayShort&, 
    ArrayShort&,
    ArrayShort&,
    ArrayBool&,
    ArrayBool&,
    ArrayShort&,
    const char, 
    const unsigned int,
    const unsigned int,
    const unsigned int,
    const unsigned int,
    const float);

  void removeCoherentNoise2D(
    ArrayFloat&,
    const ArrayFloat&, 
    ArrayFloat&,
    ArrayFloat&,
    ArrayBool&,
    ArrayBool&,
    ArrayFloat&,
    const char, 
    const unsigned int,
    const unsigned int,
    const unsigned int,
    const unsigned int,
    const float);

  void removeCoherentNoise2D(
    ArrayDouble&, 
    const ArrayDouble&, 
    ArrayDouble&,
    ArrayDouble&,
    ArrayBool&,
    ArrayBool&,
    ArrayDouble&,
    const char, 
    const unsigned int,
    const unsigned int,
    const unsigned int,
    const unsigned int,
    const float);

/// Default destructor
~Denoising(){}

private:
  template <typename T>
  void getSelectVals(
    const std::vector<std::vector<T> >& waveforms,
    const std::vector<std::vector<T> >& morphedWaveforms,
    ArrayBool& selectVals,
    ArrayBool& roi,
    const unsigned int window,
    const float thresholdFactor
  );

  template <typename T>
  void removeCoherentNoise1D(
    std::vector<std::vector<T> >& waveLessCoherent, 
    const std::vector<std::vector<T> >& filteredWaveforms, 
    std::vector<std::vector<T> >& morphedWaveforms, 
    std::vector<std::vector<T> >& intrinsicRMS,
    ArrayBool& selectVals,
    ArrayBool& roi,
    std::vector<std::vector<T> >& correctedMedians,
    const char filterName='d', 
    const unsigned int grouping=64,
    const unsigned int structuringElement=5,
    const unsigned int window=0,
    const float thresholdFactor=2.5);

  template <typename T>
  void removeCoherentNoise2D(
    std::vector<std::vector<T> >& waveLessCoherent, 
    const std::vector<std::vector<T> >& filteredWaveforms,
    std::vector<std::vector<T> >& morphedWaveforms, 
    std::vector<std::vector<T> >& intrinsicRMS,
    ArrayBool& selectVals,
    ArrayBool& roi,
    std::vector<std::vector<T> >& correctedMedians,
    const char filterName='g',
    const unsigned int grouping=64, 
    const unsigned int structuringElementx=5,
    const unsigned int structuringElementy=20,
    const unsigned int window=0,
    const float thresholdFactor=2.5);
  
};
}

#endif
/** @} */ // end of doxygen group 


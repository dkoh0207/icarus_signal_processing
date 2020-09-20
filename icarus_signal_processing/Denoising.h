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
#include "MorphologicalFunctions1D.h"
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
  Denoising() : fMPVec(4096) {}

  void getSelectVals(ArrayShort::const_iterator,
                     ArrayShort::const_iterator,
                     ArrayBool::iterator,
                     ArrayBool::iterator,
                     VectorFloat::const_iterator,
                     const unsigned int,
                     const unsigned int);

  void getSelectVals(ArrayFloat::const_iterator,
                     ArrayFloat::const_iterator,
                     ArrayBool::iterator,
                     ArrayBool::iterator,
                     VectorFloat::const_iterator,
                     const unsigned int,
                     const unsigned int);

  void getSelectVals(ArrayDouble::const_iterator,
                     ArrayDouble::const_iterator,
                     ArrayBool::iterator,
                     ArrayBool::iterator,
                     VectorFloat::const_iterator,
                     const unsigned int,
                     const unsigned int);

  void removeCoherentNoise1D(
    ArrayShort::iterator,
    ArrayShort::const_iterator,
    ArrayShort::iterator,
    ArrayShort::iterator,
    ArrayBool::iterator,
    ArrayBool::iterator,
    ArrayShort::iterator,
    FilterFunctionVec::const_iterator,
    VectorFloat::const_iterator,
    const unsigned int,
    const unsigned int,
    const unsigned int 
  );

  void removeCoherentNoise1D(
    ArrayFloat::iterator,
    ArrayFloat::const_iterator,
    ArrayFloat::iterator,
    ArrayFloat::iterator,
    ArrayBool::iterator,
    ArrayBool::iterator,
    ArrayFloat::iterator,
    FilterFunctionVec::const_iterator,
    VectorFloat::const_iterator,
    const unsigned int,
    const unsigned int,
    const unsigned int 
  );

  void removeCoherentNoise1D(
    ArrayDouble::iterator, 
    ArrayDouble::const_iterator,
    ArrayDouble::iterator,
    ArrayDouble::iterator,
    ArrayBool::iterator,
    ArrayBool::iterator,
    ArrayDouble::iterator,
    FilterFunctionVec::const_iterator,
    VectorFloat::const_iterator,
    const unsigned int,
    const unsigned int,
    const unsigned int
  );

  void removeCoherentNoise2D(
    ArrayShort&,
    const ArrayShort&, 
    ArrayShort&,
    ArrayShort&,
    ArrayBool&,
    ArrayBool&,
    ArrayShort&,
    VectorFloat&,
    const char, 
    const unsigned int,
    const unsigned int,
    const unsigned int,
    const unsigned int);

  void removeCoherentNoise2D(
    ArrayFloat&,
    const ArrayFloat&, 
    ArrayFloat&,
    ArrayFloat&,
    ArrayBool&,
    ArrayBool&,
    ArrayFloat&,
    VectorFloat&,
    const char, 
    const unsigned int,
    const unsigned int,
    const unsigned int,
    const unsigned int);

  void removeCoherentNoise2D(
    ArrayDouble&, 
    const ArrayDouble&, 
    ArrayDouble&,
    ArrayDouble&,
    ArrayBool&,
    ArrayBool&,
    ArrayDouble&,
    VectorFloat&,
    const char, 
    const unsigned int,
    const unsigned int,
    const unsigned int,
    const unsigned int);

/// Default destructor
~Denoising(){}

private:
  template <typename T>
  void getSelectVals(typename std::vector<std::vector<T>>::const_iterator waveforms,
                     typename std::vector<std::vector<T>>::const_iterator morphedWaveforms,
                     ArrayBool::iterator                                  selectVals,
                     ArrayBool::iterator                                  roi,
                     VectorFloat::const_iterator                          thresholdVec,
                     const unsigned int                                   numChannels,
                     const unsigned int                                   window
  );

  template <typename T> T getMedian(      typename std::vector<T>&, const unsigned int) const;
  template <typename T> T getMostProbable(typename std::vector<T>&, const unsigned int);

  template <typename T>
  void removeCoherentNoise1D(
    typename std::vector<std::vector<T>>::iterator       waveLessCoherent, 
    typename std::vector<std::vector<T>>::const_iterator filteredWaveforms, 
    typename std::vector<std::vector<T>>::iterator       morphedWaveforms, 
    typename std::vector<std::vector<T>>::iterator       intrinsicRMS,
    ArrayBool::iterator                                  selectVals,
    ArrayBool::iterator                                  roi,
    typename std::vector<std::vector<T> >::iterator      correctedMedians,
    FilterFunctionVec::const_iterator                    filterFunctions,
    VectorFloat::const_iterator                          thresholdFactorVec,
    const unsigned int                                   numChannels=64,
    const unsigned int                                   grouping=64,
    const unsigned int                                   window=0);

  template <typename T>
  void removeCoherentNoise2D(
    std::vector<std::vector<T> >& waveLessCoherent, 
    const std::vector<std::vector<T> >& filteredWaveforms,
    std::vector<std::vector<T> >& morphedWaveforms, 
    std::vector<std::vector<T> >& intrinsicRMS,
    ArrayBool& selectVals,
    ArrayBool& roi,
    std::vector<std::vector<T> >& correctedMedians,
    VectorFloat& thresholdVec,
    const char filterName='g',
    const unsigned int grouping=64, 
    const unsigned int structuringElementx=5,
    const unsigned int structuringElementy=20,
    const unsigned int window=0);

    // The code for the most probable calculation will need a std vector
    // We don't wnat to allocated/deallocate each call so have master copy here
    std::vector<int> fMPVec;
  
};
}

#endif
/** @} */ // end of doxygen group 


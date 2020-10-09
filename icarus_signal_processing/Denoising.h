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
#include "MorphologicalFunctions2D.h"
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

    void removeCoherentNoise1D(ArrayFloat::iterator,
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

    void removeCoherentNoise2D(ArrayFloat::iterator,
                               ArrayFloat::const_iterator,
                               ArrayFloat::iterator,
                               ArrayFloat::iterator,
                               ArrayBool::iterator,
                               ArrayBool::iterator,
                               ArrayFloat::iterator,
                               const IMorphologicalFunctions2D*,
                               VectorFloat::const_iterator,
                               const unsigned int,
                               const unsigned int,
                               const unsigned int 
    );

    void removeCoherentNoiseHough(ArrayFloat::iterator,
                                  ArrayFloat::const_iterator,
                                  ArrayFloat::iterator,
                                  ArrayFloat::iterator,
                                  ArrayBool::iterator,
                                  ArrayBool::iterator,
                                  ArrayFloat::iterator,
                                  const IMorphologicalFunctions2D*,
                                  VectorFloat::const_iterator,
                                  const unsigned int,
                                  const unsigned int,
                                  const unsigned int 
    );

    /// Default destructor
    ~Denoising(){}

private:
    void getSelectVals(ArrayFloat::const_iterator,
                       ArrayFloat::const_iterator,
                       ArrayBool::iterator,
                       ArrayBool::iterator,
                       VectorFloat::const_iterator,
                       const unsigned int,
                       const unsigned int);

    void removeCoherentNoise(ArrayFloat::iterator,
                             ArrayFloat::const_iterator,
                             ArrayFloat::iterator,
                             ArrayBool::iterator,
                             ArrayFloat::iterator,
                             const unsigned int,
                             const unsigned int 
    );
    
    float getMedian(      std::vector<float>&, const unsigned int) const;
    float getMostProbable(std::vector<float>&, const unsigned int);
    
    // The code for the most probable calculation will need a std vector
    // We don't wnat to allocated/deallocate each call so have master copy here
    std::vector<int> fMPVec;
  
};
}

#endif
/** @} */ // end of doxygen group 


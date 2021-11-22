#ifndef icarus_signal_processingDEFS_H
#define icarus_signal_processingDEFS_H
////////////////////////////////////////////////////////////////////////
//
// File:        icarus_signal_processingDefs.h
//
//              This file provides data structure definitions for filtering
//
// Created by Tracy Usher (usher@slac.stanford.edu) on January 6, 2016
//
////////////////////////////////////////////////////////////////////////

#include <vector>
#include "tbb/concurrent_vector.h"
#include <map>

namespace icarus_signal_processing
{
      using VectorShort   = std::vector<short>;
      using VectorInt     = std::vector<int>;
      using VectorFloat   = std::vector<float>;
      using VectorDouble  = std::vector<double>;
      using VectorBool    = std::vector<bool>;

      using ArrayShort    = std::vector<VectorShort>;
      using ArrayInt      = std::vector<VectorInt>;
      using ArrayFloat    = std::vector<VectorFloat>;
      using ArrayDouble   = std::vector<VectorDouble>;
      using ArrayBool     = std::vector<VectorBool>;

      using ConcurrentVectorShort = tbb::concurrent_vector<short>;
      using ConcurrentVectorInt = tbb::concurrent_vector<int>;
      using ConcurrentVectorFloat = tbb::concurrent_vector<float>;
      using ConcurrentVectorDouble = tbb::concurrent_vector<double>;
      using ConcurrentVectorBool = tbb::concurrent_vector<bool>;

      using ConcurrentArrayShort = tbb::concurrent_vector<ConcurrentVectorShort>;
      using ConcurrentArrayInt = tbb::concurrent_vector<ConcurrentVectorInt>;
      using ConcurrentArrayFloat = tbb::concurrent_vector<ConcurrentVectorFloat>;
      using ConcurrentArrayDouble = tbb::concurrent_vector<ConcurrentVectorDouble>;
      using ConcurrentArrayBool = tbb::concurrent_vector<ConcurrentVectorBool>;

} // end sigproc namespace
#endif

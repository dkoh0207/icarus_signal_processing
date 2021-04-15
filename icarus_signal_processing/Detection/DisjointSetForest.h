/**
 * \file DisjointSetForest.h
 *
 * \ingroup sigproc_tools
 * 
 * \brief Class def header for a class DisjointSetForest
 *
 * @author koh0207
 */

/** \addtogroup sigproc_tools

    @{*/
#ifndef __SIGPROC_TOOLS_DISJOINTSETFOREST_H__
#define __SIGPROC_TOOLS_DISJOINTSETFOREST_H__

#include <vector>
#include <iostream>
#include <string>

namespace icarus_signal_processing {

  /**
     \class DisjointSetForest
     User defined class DisjointSetForest ... these comments are used to generate
     doxygen documentation!
  */

  class DisjointSetForest{
    
    public:
      
      /// Default constructor
      DisjointSetForest(const unsigned int size) {
        parent.resize(size+1);
        rank.resize(size+1);
        this->size = size+1;
      }
      
      /// Default destructor
      ~DisjointSetForest(){}

      std::vector<int> parent;
      std::vector<int> rank;
      size_t size = 0;

      void MakeSet();

      void MakeSet(const std::vector<int>& strongEdges);

      void Union(const int x, const int y);

      int Find(const int x);
  };
}

#endif
/** @} */ // end of doxygen group 


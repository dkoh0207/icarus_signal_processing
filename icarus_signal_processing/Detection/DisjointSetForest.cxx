#ifndef __SIGPROC_TOOLS_DISJOINTSETFOREST_CXX__
#define __SIGPROC_TOOLS_DISJOINTSETFOREST_CXX__

#include "DisjointSetForest.h"


int icarus_signal_processing::DisjointSetForest::Find(const int x) 
{
    if (x == parent[x]) return x;
    else {
        int rep = Find(parent[x]);
        parent[x] = rep; // Path compression
        return rep;
    }
}

void icarus_signal_processing::DisjointSetForest::MakeSet()
{
    for (int i=0; i< (int) size; ++i) {
        parent[i] = i;
    }
    return;
}

void icarus_signal_processing::DisjointSetForest::MakeSet(const std::vector<int>& strongEdges)
{
    if (strongEdges.size() < 1) {
        std::string msg = "When constructing disjoint set parent with list of root "
                          "set members, list must contain at least one entry. Returning";
        std::cout << msg << std::endl;
        return;
    }
    int root = strongEdges[0];

    for (int i=0; i< (int) size; ++i) {
        parent[i] = i;
    }

    for (const int& x : strongEdges) {
        parent[x] = root;
    }
    return;
}

void icarus_signal_processing::DisjointSetForest::Union(const int x, const int y)
{
    int repX = Find(x);
    int repY = Find(y);

    if (repX == repY) return;

    else if (repX == (int) (size-1)) parent[repY] = repX;

    else if (repY == (int) (size-1)) parent[repX] = repY;

    else {
        int rankX = rank[repX];
        int rankY = rank[repY];

        if (rankX < rankY) {
            parent[repX] = repY;
            // std::cout << repX << " -> " << repY << std::endl;
        }
        else if (rankX > rankY) {
            parent[repY] = repX;
            // std::cout << repY << " -> " << repX << std::endl;
        }
        else {
            parent[repX] = repY;
            rank[repY] = rank[repY] + 1;
            // std::cout << repX << " -> " << repY << std::endl;
        }
    }
} 

#endif

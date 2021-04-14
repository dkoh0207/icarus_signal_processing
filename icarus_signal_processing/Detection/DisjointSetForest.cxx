#ifndef __SIGPROC_TOOLS_DISJOINTSETFOREST_CXX__
#define __SIGPROC_TOOLS_DISJOINTSETFOREST_CXX__

#include "DisjointSetForest.h"


int icarus_signal_processing::DisjointSetForest::Find(int x) 
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

void icarus_signal_processing::DisjointSetForest::MakeSet(std::vector<int>& strongEdges)
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

    for (int& x : strongEdges) {
        parent[x] = root;
    }
    return;
}

void icarus_signal_processing::DisjointSetForest::Union(int x, int y)
{
    int repX = Find(x);
    int repY = Find(y);

    if (repX == repY) return;

    else {
        int rankX = rank[repX];
        int rankY = rank[repY];

        if (rankX < rankY) parent[repX] = repY;
        else if (rankX > rankY) parent[repY] = repX;
        else {
            parent[repX] = repY;
            rank[repY] = rank[repY] + 1;
        }
    }
} 

#endif

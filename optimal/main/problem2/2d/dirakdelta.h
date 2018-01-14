#ifndef DIRAKDELTA_H
#define DIRAKDELTA_H

#include <grid/grid.h>
#include <vector>
#include <math.h>

class DirakDelta
{
public:

    void payla(const SpaceNodePDE &point, const Dimension &dimension, std::vector<SpaceNodePDE> &nodes);

    SpaceNodePDE point;
    Dimension dimension;
};

#endif // DIRAKDELTA_H

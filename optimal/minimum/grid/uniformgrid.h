#ifndef UNIFORMGRID_H
#define UNIFORMGRID_H

#include <global.h>

struct GridNode
{
    double x;
    double y;
    double z;
    double t;

    int xi;
    int yi;
    int zi;
    int ti;
};

class MINIMUMSHARED_EXPORT UniformGrid
{
public:
    UniformGrid(double step = 0.0, int min = 0, int max = 0);

    double step() const;
    int minN() const;
    int maxN() const;
    int sizeN() const;

    bool isGridSet() const;

private:
    double mstep;
    int mminN;
    int mmaxN;
};

#endif // UNIFORMGRID_H

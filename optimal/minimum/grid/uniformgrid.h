#ifndef UNIFORMGRID_H
#define UNIFORMGRID_H

#include <global.h>

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

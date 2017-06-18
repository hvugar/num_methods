#ifndef GRID_H
#define GRID_H

#include <global.h>
#include <vector>

struct MINIMUMSHARED_EXPORT SpaceNode
{
    unsigned int i;
    unsigned int j;
    unsigned int k;

    double x;
    double y;
    double z;
};

struct MINIMUMSHARED_EXPORT TimeNode
{
    unsigned int i;
    double t;
};

class MINIMUMSHARED_EXPORT Dimension
{
public:
    enum SpaceDimension
    {
        Dim1 = 0,
        Dim2 = 1,
        Dim3 = 2
    };

    Dimension(double step=0.01, unsigned int maxN=100, unsigned int minN = 0);

    virtual double step() const;
    virtual unsigned int minN() const;
    virtual unsigned int maxN() const;
    virtual unsigned int sizeN() const;

protected:
    double mstep;
    unsigned int mmaxN;
    unsigned int mminN;
};

#endif // GRID_H

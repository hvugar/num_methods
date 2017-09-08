#ifndef LOADEDHEATEQUATION_H
#define LOADEDHEATEQUATION_H

#include "iloadedheatequation.h"

#define SAMPLE_1

class LoadedHeatEquation : public ILoadedHeatEquation
{
public:
    static void Main(int argc, char** argv);
    double U(const SpaceNodePDE &sn,const TimeNodePDE &tn) const;

    virtual double initial(const SpaceNodePDE &sn) const;
    virtual double boundary(const SpaceNodePDE &sn, const TimeNodePDE &tn, BoundaryType boundary) const;
    virtual double f(const SpaceNodePDE &sn, const TimeNodePDE &tn) const;

    virtual double g(const TimeNodePDE &tn) const;
    virtual double h(const TimeNodePDE &tn) const;

    virtual void layerInfo(const DoubleVector &, unsigned int) const;
};

#endif // LOADEDHEATEQUATION_H

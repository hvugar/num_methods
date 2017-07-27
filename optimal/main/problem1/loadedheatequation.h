#ifndef LOADEDHEATEQUATION_H
#define LOADEDHEATEQUATION_H

#include "iloadedheatequation.h"

#define SAMPLE_1

class LoadedHeatEquation : public ILoadedHeatEquation
{
public:
    static void Main(int argc, char** argv);
    double U(const SpaceNode &sn,const TimeNode &tn) const;

    virtual double initial(const SpaceNode &sn) const;
    virtual double boundary(const SpaceNode &sn, const TimeNode &tn, BoundaryType boundary) const;
    virtual double f(const SpaceNode &sn, const TimeNode &tn) const;

    virtual double g(const TimeNode &tn) const;
    virtual double h(const TimeNode &tn) const;

    virtual void layerInfo(const DoubleVector &, unsigned int) const;
};

#endif // LOADEDHEATEQUATION_H

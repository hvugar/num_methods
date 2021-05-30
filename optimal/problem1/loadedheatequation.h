#ifndef LOADEDHEATEQUATION_H
#define LOADEDHEATEQUATION_H

#include "iloadedheatequation.h"

#define SAMPLE_1

class LoadedHeatEquation : public ILoadedHeatEquation
{
public:
    static void Main(int argc, char** argv);
    double U(const SpaceNodePDE &sn,const TimeNodePDE &tn) const;

protected:
    virtual double initial(const SpaceNodePDE &sn, InitialCondition condition) const;
    virtual double boundary(const SpaceNodePDE &sn, const TimeNodePDE &tn, BoundaryConditionPDE &condition) const;
    virtual double f(const SpaceNodePDE &sn, const TimeNodePDE &tn) const;

public:
    virtual void layerInfo(const DoubleVector &, const TimeNodePDE &) const;
    virtual void layerInfo(const DoubleMatrix &, const TimeNodePDE &) const {}

    virtual double g(const TimeNodePDE &tn) const;
    virtual double h(const TimeNodePDE &tn) const;

    virtual Dimension timeDimension() const { return Dimension(); }
    virtual Dimension spaceDimensionX() const { return Dimension(); }
    virtual Dimension spaceDimensionY() const { return Dimension(); }
    virtual Dimension spaceDimensionZ() const { return Dimension(); }
};

#endif // LOADEDHEATEQUATION_H

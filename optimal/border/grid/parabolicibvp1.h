#ifndef PARABOLICIBVP1_H
#define PARABOLICIBVP1_H

#include <grid/pibvp.h>
#include <time.h>

#define SAMPLE_0

class MINIMUMSHARED_EXPORT ParabolicIBVP1 : public ParabolicIBVP
{
public:
    void static Main(int argc, char* argv[]);

public:
    double U(unsigned int i, unsigned int j) const;
    double U(const SpaceNodePDE &sn, const TimeNodePDE& tn) const;

protected:
    virtual double initial(const SpaceNodePDE &sn) const;
    virtual double boundary(const SpaceNodePDE &sn, const TimeNodePDE& tn, BoundaryType boundary = Unused) const;
    virtual double f(const SpaceNodePDE &sn, const TimeNodePDE &tn) const;
    virtual double a(const SpaceNodePDE &sn, const TimeNodePDE &tn) const;

    virtual void layerInfo(const DoubleVector &, unsigned int) const;
};

#endif // PARABOLICIBVP1_H

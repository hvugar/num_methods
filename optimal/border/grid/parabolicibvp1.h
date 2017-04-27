#ifndef PARABOLICIBVP1_H
#define PARABOLICIBVP1_H

#include <grid/pibvp.h>
#include <time.h>

#define SAMPLE_6

class MINIMUMSHARED_EXPORT ParabolicIBVP1 : public ParabolicIBVP
{
public:
    ParabolicIBVP1();
    double U(unsigned int i, unsigned int j) const;

protected:
    virtual double initial(const SpaceNode &sn) const;
    virtual double boundary(const SpaceNode &sn, const TimeNode& tn, BoundaryType boundary = Unused) const;
    virtual double f(const SpaceNode &sn, const TimeNode &tn) const;
    virtual double a(const SpaceNode &sn, const TimeNode &tn) const;

    virtual void layerInfo(const DoubleVector &, unsigned int) const;

public:
    void static Main(int argc, char* argv[]);
};

#endif // PARABOLICIBVP1_H

#ifndef PARABOLICIBVP2_H
#define PARABOLICIBVP2_H

#include <grid/pibvp.h>
#include <time.h>

#define SAMPLE_3

class MINIMUMSHARED_EXPORT ParabolicIBVP2 : public ParabolicIBVP
{
public:
    ParabolicIBVP2();

protected:
    virtual double initial(const SpaceNodePDE &sn) const;
    virtual double boundary(const SpaceNodePDE &sn, const TimeNodePDE& tn, BoundaryType boundary = Unused) const;
    virtual double f(const SpaceNodePDE &sn, const TimeNodePDE &tn) const;
    virtual double a(const SpaceNodePDE &sn, const TimeNodePDE &tn) const;

public:
    void static Main(int argc, char* argv[]);
    double U(unsigned int i, unsigned int j, unsigned int k) const;
};

#endif // PARABOLICIBVP2_H

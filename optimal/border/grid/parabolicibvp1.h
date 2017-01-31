#ifndef PARABOLICIBVP1_H
#define PARABOLICIBVP1_H

#include <grid/pibvp.h>

#define SAMPLE_4

class MINIMUMSHARED_EXPORT ParabolicIBVP1 : public ParabolicIBVP
{
public:
    ParabolicIBVP1(const GridPDE &grid);

    double U(unsigned int i, unsigned int j) const;

protected:
    virtual double initial(unsigned int n) const;
    virtual double boundary(unsigned int m, BoundaryType boundary) const;
    virtual double f(unsigned int n, unsigned int m) const;
    virtual double a(unsigned int n, unsigned int m) const;

    virtual double initial(const SpaceNode &sn) const;
    virtual double boundary(const SpaceNode &sn, const TimeNode& tn) const;
    virtual double f(const SpaceNode &sn, const TimeNode &tn) const;
    virtual double a(const SpaceNode &sn, const TimeNode &tn) const;

public:
    void static Main(int agrc, char* argv[]);
};

#endif // PARABOLICIBVP1_H

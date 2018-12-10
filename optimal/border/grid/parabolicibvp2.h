#ifndef PARABOLICIBVP2_H
#define PARABOLICIBVP2_H

#include <grid/pibvp.h>

class MINIMUMSHARED_EXPORT ParabolicIBVP2 : public ParabolicIBVP
{
protected:
    virtual double initial(const SpaceNodePDE &sn) const;
    virtual double boundary(const SpaceNodePDE &sn, const TimeNodePDE& tn) const;
    virtual double f(const SpaceNodePDE &sn, const TimeNodePDE &tn) const;
    virtual double a(const SpaceNodePDE &sn, const TimeNodePDE &tn) const;

public:
    void static Main(int argc, char* argv[]);
    double U(const SpaceNodePDE &sn, const TimeNodePDE &tn) const;
};

class MINIMUMSHARED_EXPORT CCParabolicIBVP1 : public CCParabolicIBVP
{
protected:
    virtual double initial(const SpaceNodePDE &sn) const;
    virtual double boundary(const SpaceNodePDE &sn, const TimeNodePDE &tn) const;
    virtual double f(const SpaceNodePDE &sn, const TimeNodePDE &tn) const;

    virtual double U(const SpaceNodePDE &sn, const TimeNodePDE &tn) const;

public:
    static void Main(int argc, char* argv[]);
};

#endif // PARABOLICIBVP2_H

#ifndef HYPERBOLICIBVP1_H
#define HYPERBOLICIBVP1_H

#include <grid/hibvp.h>
#include <time.h>

class MINIMUMSHARED_EXPORT HyperbolicIBVP1 : public HyperbolicIBVP
{
public:
    HyperbolicIBVP1();
    double U(unsigned int i, unsigned int j) const;

protected:
    virtual double initial1(const SpaceNodePDE &sn) const;
    virtual double initial2(const SpaceNodePDE &sn) const;
    virtual double boundary(const SpaceNodePDE &sn, const TimeNodePDE &tn, BoundaryType boundary = Unused) const;
    virtual double f(const SpaceNodePDE &sn, const TimeNodePDE &tn) const;
    virtual double a(const SpaceNodePDE &sn, const TimeNodePDE &tn) const;
public:
    static void Main(int argc, char* argv[]);
};

class MINIMUMSHARED_EXPORT HyperbolicIBVP2 : public IHyperbolicIBVP
{
public:
    static void Main(int argc, char* argv[]);
protected:
    virtual double initial1(const SpaceNodePDE &sn) const;
    virtual double initial2(const SpaceNodePDE &sn) const;
    virtual double boundary(const SpaceNodePDE &sn, const TimeNodePDE &tn, BoundaryType boundary = Unused) const;
    virtual double f(const SpaceNodePDE &sn, const TimeNodePDE &tn) const;
private:
    double a = 1.0;
};

#endif // HYPERBOLICIBVP1_H

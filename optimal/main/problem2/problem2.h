#ifndef IPROBLEM2_H
#define IPROBLEM2_H

#include <grid/pibvp.h>


class IProblem2 : public IParabolicIBVP
{
public:
    IProblem2();

    void gridMethod(DoubleVector &u) const;


protected:
    virtual double initial(const SpaceNodePDE &sn) const = 0;
    virtual double boundary(const SpaceNodePDE &sn, const TimeNodePDE &tn, BoundaryType boundary = Unused) const = 0;
    virtual double f(const SpaceNodePDE &sn, const TimeNodePDE &tn) const = 0;

private:
    double a;
    double lambda0;
    double theta;

    unsigned int Lc;
    unsigned int Lo;
    DoubleMatrix k;
    DoubleMatrix z;
    DoubleVector xi;
    DoubleVector eta;

    double U(double x, double t) const;
};

#endif // IPROBLEM2_H

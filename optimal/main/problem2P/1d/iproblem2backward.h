#ifndef IPROBLEM2_BACKWARD_H
#define IPROBLEM2_BACKWARD_H

#include <grid/pibvp.h>
#include <printer.h>

class IProblem2Backward : public IParabolicIBVP
{
public:
    IProblem2Backward();

    void gridMethod(DoubleMatrix &p) const;

    DoubleVector uT;
    DoubleVector U;
    DoubleVector mu;

protected:
    virtual double initial(const SpaceNodePDE &sn) const;
    virtual double boundary(const SpaceNodePDE &sn, const TimeNodePDE &tn, BoundaryType boundary = Unused) const;
    virtual double f(const SpaceNodePDE &sn, const TimeNodePDE &tn) const;

    virtual double h(const SpaceNodePDE &sn) const;
    virtual double g0(const TimeNodePDE &tn) const;
    virtual double g1(const TimeNodePDE &tn) const;

    virtual double delta(const SpaceNodePDE &sn, unsigned int j) const;

    virtual double P(double x, double t) const;

public:
    double a;
    double lambda0;
    double lambda1;
    double lambda2;
    double theta;

    unsigned int Lc;
    unsigned int Lo;
    DoubleMatrix k;
    DoubleMatrix z;
    DoubleVector xi;
    DoubleVector eta;
};

#endif // IPROBLEM2_BACKWARD_H

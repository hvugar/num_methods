#ifndef IPROBLEM2_H
#define IPROBLEM2_H

#include <grid/pibvp.h>
#include <printer.h>


class IProblem2 : public IParabolicIBVP
{
public:
    static void Main(int argc, char* argv[]);

    IProblem2();

    void gridMethod(DoubleMatrix &u) const;


protected:
    virtual double initial(const SpaceNodePDE &sn) const;
    virtual double boundary(const SpaceNodePDE &sn, const TimeNodePDE &tn, BoundaryType boundary = Unused) const;
    virtual double f(const SpaceNodePDE &sn, const TimeNodePDE &tn) const;

    virtual double g0(const TimeNodePDE &tn) const;
    virtual double g1(const TimeNodePDE &tn) const;
    virtual double delta(const SpaceNodePDE &sn, unsigned int i) const;

private:
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

    double U(double x, double t) const;
};

#endif // IPROBLEM2_H

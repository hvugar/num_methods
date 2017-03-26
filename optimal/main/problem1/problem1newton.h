#ifndef PROBLEM1_NEWTON_H
#define PROBLEM1_NEWTON_H

#include <grid/pibvp.h>
#include <grid/bpibvp.h>

class Problem1L2;

class Problem1NewtonF : public IParabolicIBVP
{
public:
    Problem1NewtonF();
    virtual ~Problem1NewtonF();

    void calculate(DoubleMatrix &u, const DoubleVector &k, const DoubleVector &z, const DoubleVector &e);

    double fi;
    double tt;
    double a;
    double lambda0;
    double lambda1;
    double lambda2;

protected:
    virtual double initial(const SpaceNode &sn) const;
    virtual double boundary(const SpaceNode &sn, const TimeNode &tn, BoundaryType boundary = Unused) const;
    virtual double f(const SpaceNode &sn, const TimeNode &tn) const;

private:
    void qovmaFirstRowM(double *a, double *b, double *c, double *d, double *x, unsigned int n, double *e) const;
};

class Problem1NewtonB : public IBackwardParabolicIBVP
{
public:
    Problem1NewtonB();
    virtual ~Problem1NewtonB();

    void calculate(DoubleMatrix &p, const DoubleMatrix &u, const DoubleVector &k, const DoubleVector &z, const DoubleVector &e);

    double a;
    double lambda0;
    double lambda1;
    double lambda2;

    double alpha0;

protected:
    virtual double initial(const SpaceNode &sn) const;
    virtual double boundary(const SpaceNode &sn, const TimeNode &tn, BoundaryType boundary = Unused) const;
    virtual double f(const SpaceNode &sn, const TimeNode &tn) const;

private:
    void qovmaFirstColM(double *a, double *b, double *c, double *d, double *x, unsigned int n, double *e) const;
};

#endif // PROBLEM1NEWTONF_H

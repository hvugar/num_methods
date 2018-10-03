#ifndef PROBLEM2_H
#define PROBLEM2_H

#include "iproblem2forward.h"
#include "iproblem2backward.h"

#include <function.h>
#include <gradient.h>

class Problem2Forward : public IProblem2Forward
{
protected:
    virtual double f(const SpaceNodePDE &sn UNUSED_PARAM, const TimeNodePDE &tn UNUSED_PARAM) const { return 0.0; }
    virtual double g0(const TimeNodePDE &tn UNUSED_PARAM) const { return 0.0; }
    virtual double g1(const TimeNodePDE &tn UNUSED_PARAM) const { return 0.0; }
};

class Problem2Backward : public IProblem2Backward
{
protected:
    virtual double f(const SpaceNodePDE &sn UNUSED_PARAM, const TimeNodePDE &tn UNUSED_PARAM) const { return 0.0; }
    virtual double h(const SpaceNodePDE &sn UNUSED_PARAM) const { return 0.0; }
    virtual double g0(const TimeNodePDE &tn UNUSED_PARAM) const { return 0.0; }
    virtual double g1(const TimeNodePDE &tn UNUSED_PARAM) const { return 0.0; }
};

class Problem2 : public RnFunction, public IGradient
{
public:
    static void Main(int argc, char* argv[]);

    Problem2();

    DoubleVector U;

public:
    virtual double fx(const DoubleVector &x) const;
    virtual void gradient(const DoubleVector &x, DoubleVector &g) const;

    virtual double integral(const DoubleMatrix &u) const;

    virtual double mu(unsigned int n) const;

private:
    Dimension mTimeDimension;
    Dimension mSpaceDimension;

    Problem2Forward forward;
    Problem2Backward backward;

    double alpha0;

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

    void array2Parameters(const DoubleVector &prms, DoubleMatrix &k, DoubleMatrix &z, DoubleVector &xi, DoubleVector &eta) const;
    void paremeters2Array(const DoubleMatrix &k, const DoubleMatrix &z, const DoubleVector &xi, const DoubleVector &eta, DoubleVector &prms) const;
};

#endif // PROBLEM2_H

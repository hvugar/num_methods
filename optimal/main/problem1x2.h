#ifndef PROBLEM1X2_H
#define PROBLEM1X2_H

#include <cmethods.h>
#include <function.h>
#include <printer.h>
#include <gradient_cjt.h>
#include <vector>
#include <rungekutta.h>
#include <projection.h>
#include <limits.h>
#include <float.h>

class Problem1X2 : public RnFunction, public IGradient, public IPrinter, public ConjugateGradient, public IProjection
{
public:
    Problem1X2();
    virtual ~Problem1X2();

    virtual double fx(const DoubleVector &x);
    virtual void gradient(const DoubleVector &x, DoubleVector &g);
    virtual void print(unsigned int iteration, const DoubleVector &k, const DoubleVector &g, double alpha, RnFunction *fn) const;
    virtual void project(DoubleVector &x, int index);

    virtual double vm(unsigned int j) const;
    //virtual double vl(unsigned int j) const;
    virtual double vr(unsigned int j) const;

    virtual double initial(unsigned int i) const;
    virtual double mu(unsigned int i) const;
    void calculateU(DoubleMatrix &m, double ht, double hx, unsigned int M, unsigned int N, double alpha, double lambda0, double lambdal, double a);
    void calculateP(DoubleMatrix &m, double ht, double hx, unsigned int M, unsigned int N, double alpha, double lambda0, double lambdal, double a);
    void calculateV(const DoubleVector &k);
    //double vs(double j) const;

    static void Main(int argc, char* argv[]);

public:
    double t0;
    double t1;
    double x0;
    double x1;
    double hx;
    double ht;
    unsigned int M;
    unsigned int N;
    double a;

    double alpha;
    double lambda0;
    double lambdal;

    unsigned int L;
    DoubleVector k;
    DoubleVector z;
    double Te;
    double Ti;

    double alpha1;
    double alpha2;

    DoubleVector V;
    const DoubleVector *pxi;
    const DoubleMatrix *pu;
    const DoubleMatrix *pp;

    DoubleVector xis;
};

#endif // PROBLEM1X2_H

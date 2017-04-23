#ifndef IPROBLEM1_H
#define IPROBLEM1_H

#include <math.h>
#include <function.h>
#include <gradient.h>
#include <vector2d.h>
#include <matrix2d.h>
#include <cmethods.h>
#include <gradient/igradient.h>
#include <printer.h>

class IProblem1 : protected RnFunction, protected IGradient, public IPrinter
{
public:
    unsigned int L;
    unsigned int N;
    unsigned int M;
    double hx;
    double ht;
    double a;
    double lambda0;
    double lambda1;
    double lambda2;

    double alpha0 = 1.0;
    double alpha1;
    double alpha2;
    double alpha3;

    DoubleVector vfi;
    DoubleVector vtt;
    double fi;
    double tt;

    DoubleVector V;
    const DoubleVector *py;

    double R = 0.0;
    double vmin;
    double vmax;
    double d0;
    double d1;

    bool optimizeK;
    bool optimizeZ;
    bool optimizeE;

    DoubleVector K;
    DoubleVector Z;
    DoubleVector E;

    DoubleVector k0;
    DoubleVector z0;
    DoubleVector e0;
    unsigned int DD = 10;

    double hk;
    double hz;
    double he;

    bool withError = false;
    bool hello = false;
    double persent = 0.01;

protected:
    virtual double fx(const DoubleVector &y) const;
    virtual void gradient(const DoubleVector &y, DoubleVector &g);

    virtual double integral(const DoubleMatrix &u) const;
    virtual double norm(const DoubleVector &k, const DoubleVector &z, const DoubleVector &e) const;
    virtual double penalty(const DoubleVector &k, const DoubleVector &z, const DoubleVector &e, const DoubleMatrix &u) const;

    virtual void print(unsigned int i, const DoubleVector &y, const DoubleVector &g, double fx, GradientMethod::MethodResult result) const;

    virtual void calculateU(DoubleMatrix &u, const DoubleVector &k, const DoubleVector &z, const DoubleVector &e) const;
    virtual void calculateU1(DoubleMatrix &u, const DoubleVector &k, const DoubleVector &z, const DoubleVector &e) const;
    virtual void calculateU2(DoubleMatrix &u, const DoubleVector &k, const DoubleVector &z, const DoubleVector &e) const;
    virtual void calculateP(DoubleMatrix &p, const DoubleMatrix &u, const DoubleVector &k, const DoubleVector &z, const DoubleVector &e) const;

    virtual double initial(unsigned int n) const;
    virtual double mu(unsigned int n) const;

    virtual void qovmaFirstColM(double *a, double *b, double *c, double *d, double *x, unsigned int n, double *E) const;
    virtual void qovmaFirstRowM(double *a, double *b, double *c, double *d, double *x, unsigned int n, double *E) const;
    virtual void getParameters(DoubleVector &k, DoubleVector &z, DoubleVector &e, const DoubleVector &y) const;

    virtual double delta(unsigned int n, const DoubleVector &e, unsigned int s) const;
    virtual double vf(unsigned int m, const DoubleVector &k, const DoubleVector &z, const DoubleVector &e, const DoubleMatrix &u) const;
    virtual double sgn_min(unsigned int m, const DoubleVector &k, const DoubleVector &z, const DoubleVector &e, const DoubleMatrix &u) const;
    virtual double vd0(unsigned int m, const DoubleVector &k, const DoubleVector &z, const DoubleVector &e, const DoubleMatrix &u) const;
    virtual double gf(unsigned int m, const DoubleVector &k, const DoubleVector &z, const DoubleVector &e, const DoubleMatrix &u) const;
    virtual double sign(double x) const;

    virtual double u_xi(unsigned int m, double xi, const DoubleMatrix &u) const;
    virtual double u_xi_d(unsigned int m, double xi, const DoubleMatrix &u) const;

};

#endif // IPROBLEM1_H

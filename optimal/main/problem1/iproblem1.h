#ifndef IPROBLEM1_H
#define IPROBLEM1_H

#include <matrix2d.h>
#include <vector2d.h>
#include <math.h>
#include <function.h>

class IProblem1 : protected RnFunction, protected IGradient
{
public:
    unsigned int N;
    unsigned int M;
    double hx;
    double ht;
    double a;
    double lambda0;
    double lambda1;
    double lambda2;

    DoubleVector vfi;
    DoubleVector vtt;
    double fi;
    double tt;

    unsigned int L;

    double alpha0 = 1.0;

    double alpha1;
    double alpha2;
    double alpha3;

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


protected:
    virtual double fx(const DoubleVector &y) const;
    virtual void gradient(const DoubleVector &y, DoubleVector &g);

    void calculateU(DoubleMatrix &u, const DoubleVector &k, const DoubleVector &z, const DoubleVector &e) const;
    void calculateP(DoubleMatrix &p, const DoubleMatrix &u, const DoubleVector &k, const DoubleVector &z, const DoubleVector &e) const;

    virtual double initial(unsigned int n) const;
    virtual double mu(unsigned int n) const;

    void qovmaFirstColM(double *a, double *b, double *c, double *d, double *x, unsigned int n, double *E) const;
    void qovmaFirstRowM(double *a, double *b, double *c, double *d, double *x, unsigned int n, double *E) const;
    void getParameters(DoubleVector &k, DoubleVector &z, DoubleVector &e, const DoubleVector &y) const;

    double delta(unsigned int n, const DoubleVector &e, unsigned int s) const;
    double vf(unsigned int m, const DoubleVector &k, const DoubleVector &z, const DoubleVector &e, const DoubleMatrix &u) const;
    double sgn_min(unsigned int m, const DoubleVector &k, const DoubleVector &z, const DoubleVector &e, const DoubleMatrix &u) const;
    double vd0(unsigned int m, const DoubleVector &k, const DoubleVector &z, const DoubleVector &e, const DoubleMatrix &u) const;
    double gf(unsigned int m, const DoubleVector &k, const DoubleVector &z, const DoubleVector &e, const DoubleMatrix &u) const;
    double sign(double x) const;

    double u_xi(unsigned int m, double xi, const DoubleMatrix &u) const;
    double u_xi_d(unsigned int m, double xi, const DoubleMatrix &u) const;

};

#endif // IPROBLEM1_H

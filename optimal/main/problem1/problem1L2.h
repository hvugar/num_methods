#ifndef PROBLEM1L2_H
#define PROBLEM1L2_H

#include <function.h>
#include <gradient.h>
#include <vector2d.h>
#include <matrix2d.h>
#include <printer.h>
#include <projection.h>
#include <gradient_cjt.h>
#include <cmethods.h>

#define _OPTIMIZE_K_
#define _OPTIMIZE_Z_
#define _OPTIMIZE_E_

class Problem1L2 : protected RnFunction, protected IGradient, public IPrinter, public IProjection, public R1Function
{
public:
    Problem1L2();

protected:
    virtual double fx(double t) const;

    virtual double fx(const DoubleVector &x) const;
    double integral(const DoubleVector &prm, const DoubleMatrix &u) const;
    double norm(const DoubleVector &prm) const;
    double penalty(const DoubleVector &k, const DoubleVector &z, const DoubleVector &e, const DoubleMatrix &u) const;
    double vd0(unsigned int m, const DoubleVector &k, const DoubleVector &z, const DoubleVector &e, const DoubleMatrix &u) const;
    double gf(unsigned int m, const DoubleVector &k, const DoubleVector &z, const DoubleVector &e, const DoubleMatrix &u) const;

    virtual void gradient(const DoubleVector &prm, DoubleVector &g);

    virtual void print(unsigned int i, const DoubleVector &x, const DoubleVector &g, double fx, GradientMethod::MethodResult result) const;
    virtual void project(DoubleVector &x, int index);


    void calculateU(DoubleMatrix &u) const;
    void calculateP(DoubleMatrix &p, const DoubleMatrix &u);

    double initial(unsigned int n) const;
    double mu(unsigned int i UNUSED_PARAM) const { return 1.0; }

    void getComponents(DoubleVector &k, DoubleVector &z, DoubleVector &e, const DoubleVector &x) const;
    void qovmaFirstColM(double *a, double *b, double *c, double *d, double *x, unsigned int n, double *e) const;
    void qovmaFirstRowM(double *a, double *b, double *c, double *d, double *x, unsigned int n, double *e) const;

    double vf(unsigned int m, const DoubleVector &k, const DoubleVector &z, const DoubleVector &e, const DoubleMatrix &u) const;

private:
    unsigned int L;
    unsigned int N;
    unsigned int M;
    double hx;
    double ht;

    double hk;
    double hz;
    double he;

    double fi;
    double tt;

    DoubleVector vfi;
    DoubleVector vtt;

    double a;
    double lambda0;
    double lambda1;
    double lambda2;

    double alpha0;
    double alpha1;
    double alpha2;
    double alpha3;
    double R;

    double vmin;
    double vmax;
    double d0;
    double d1;

    double zmin;
    double zmax;

    DoubleVector V;
    const DoubleVector *px;

    DoubleVector K;
    DoubleVector z;
    DoubleVector e;

    bool optimizeK;
    bool optimizeZ;
    bool optimizeE;

    FILE* file;

public:
    static void Main(int argc, char* argv[]);
};

#endif // PROBLEM1L2_H

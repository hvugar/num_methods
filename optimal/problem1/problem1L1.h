#ifndef PROBLEM1L1_H
#define PROBLEM1L1_H

#include <gradient_cjt.h>
#include <gradient_sd.h>
#include <function.h>
#include <gradient.h>
#include <vector2d.h>
#include <matrix2d.h>
#include <printer.h>
#include <projection.h>
#include <cmethods.h>
#include <float.h>
#include <vector>

#define _OPTIMIZE_K_
#define _OPTIMIZE_Z_
#define _OPTIMIZE_E_

//#define __V_NORM__

void qovmaFirstCol(double *a, double *b, double *c, double *d, double *x, unsigned int n, double *e);
void qovmaFirstRow(double *a, double *b, double *c, double *d, double *x, unsigned int n, double *e);

class Problem1L1 : public RnFunction, public IGradient, public ConjugateGradient, public IPrinter, public IProjection
{
public:
    Problem1L1();
    virtual ~Problem1L1() {}

    void initialize();

    virtual double fx(const DoubleVector &x) const;
    virtual void gradient(const DoubleVector &x, DoubleVector &g) const;

    virtual void print(unsigned int i, const DoubleVector &x, const DoubleVector &g, double fx, double alpha, GradientMethod::MethodResult result) const;
    virtual void project(DoubleVector &x, unsigned int index);

    double initial(unsigned int i) const;
    double mu(unsigned int i UNUSED_PARAM) const { return 1.0; }

    // qovma 1
    void calculateU(DoubleMatrix &u, unsigned int N, unsigned int M, double hx, double ht) const;
    void calculateUN2L2R(DoubleMatrix &u, unsigned int N, unsigned int M, double hx, double ht);
    void calculateUN4L2R(DoubleMatrix &u, unsigned int N, unsigned int M, double hx, double ht);
    void calculateP(DoubleMatrix& p, const DoubleMatrix &u);

    void getComponents(DoubleVector &k, DoubleVector &z, DoubleVector &e, const DoubleVector &x) const;
    void printNAGradinets(const DoubleVector &x0);

    double v(unsigned int j, const DoubleMatrix &u) const;

private:
    double a = 1.0;
    unsigned int L = 2;
    unsigned int N = 100;
    unsigned int M = 100;
    double hx = 0.01;
    double ht = 0.01;
    double h  = 0.01;

    double Ti = 2.0;
    double Te = 3.0;
    double alpha = 1.0;
    double lambda0 = 1.0;
    double lambdal = 1.0;

    double alpha0 = 1.0;

#ifdef __V_NORM__
    double alpha4 = 1.0;
#else
#ifdef _OPTIMIZE_K_
    double alpha1 = 1.0;
#else
    double alpha1 = 0.0;
#endif

#ifdef _OPTIMIZE_Z_
    double alpha2 = 1.0;
#else
    double alpha2 = 0.0;
#endif

#ifdef _OPTIMIZE_E_
    double alpha3 = 1.0;
#else
    double alpha3 = 0.0;
#endif
#endif

    const DoubleVector *px;
    DoubleVector V;
    DoubleVector xs;

#ifndef _OPTIMIZE_K_
    DoubleVector k;
#endif

#ifndef _OPTIMIZE_Z_
    DoubleVector z;
#endif

#ifndef _OPTIMIZE_E_
    DoubleVector e;
#endif

public:
    static void Main(int argc, char* argv[]);
};

#endif // PROBLEM1L1_H

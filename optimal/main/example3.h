#ifndef EXAMPLE3_H
#define EXAMPLE3_H

#include <gradient_cjt.h>
#include <gradient_sd.h>
#include <function.h>
#include <gradient.h>
#include <vector2d.h>
#include <matrix2d.h>
#include <printer.h>
#include <projection.h>

#define _OPTIMIZE_K_
#define _OPTIMIZE_Z_
#define _OPTIMIZE_E_

void qovmaFirstCol(double *a, double *b, double *c, double *d, double *x, unsigned int n, double *e);
void qovmaFirstRow(double *a, double *b, double *c, double *d, double *x, unsigned int n, double *e);

class Example3 : public RnFunction, public IGradient, public ConjugateGradient, public IPrinter, public IProjection
{
public:
    Example3();
    virtual ~Example3() {}

    void initialize();

    virtual double fx(const DoubleVector &x);
    virtual void gradient(const DoubleVector &x, DoubleVector &g);
    virtual void print(unsigned int iteration, const DoubleVector &x, const DoubleVector &gradient, double alpha, RnFunction *fn) const;
    virtual void project(DoubleVector &x, int index);

    double initial(unsigned int i) const;
    double mu(unsigned int i UNUSED_PARAM) const { return 1.0; }
    // qovma 1
    void calculateU(DoubleMatrix &u);
    void calculateP(DoubleMatrix& p, const DoubleMatrix &u);

    void getComponents(DoubleVector &k, DoubleVector &z, DoubleVector &e, const DoubleVector &x);
    void printNAGradinets(const DoubleVector &x0);


    //    //GaussianElimination
    //    void calculateU1(DoubleMatrix &u, const DoubleVector &x);
    //    // qovma E
    //    void calculateU2(DoubleMatrix &u, const DoubleVector &x);
    //    // teze qovma
    //    void calculateU3(DoubleMatrix &u, const DoubleVector &x);
    //    void calculate1(unsigned int N, unsigned int K, const DoubleMatrix &a, DoubleMatrix &beta, DoubleVector &qamma, DoubleVector &x);

private:
    double a = 1.0;
    unsigned int L = 2;
    unsigned int N = 1000;
    unsigned int M = 1000;
    double hx = 0.001;
    double ht = 0.001;
    double h  = 0.001;

    double Ti = 2.0;
    double Te = 3.0;
    double alpha = 1.0;
    double lambda0 = 1.0;
    double lambdal = 1.0;

    double alpha0 = 1.0;

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

#endif // EXAMPLE3_H

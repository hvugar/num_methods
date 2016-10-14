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
    void calculateU1(DoubleMatrix &u);


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
    double alpha1 = 0.0;
    double alpha2 = 0.0;
    double alpha3 = 0.0;

    const DoubleVector *px;
    DoubleVector V;
    DoubleVector xs;

    //DoubleVector k;
    //DoubleVector z;
    //DoubleVector e;

public:
    static void Main(int argc, char* argv[]);
};

#endif // EXAMPLE3_H

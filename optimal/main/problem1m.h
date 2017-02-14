#ifndef PROBLEM1M_H
#define PROBLEM1M_H

#include <function.h>
#include <gradient.h>
#include <vector2d.h>
#include <matrix2d.h>
#include <printer.h>
#include <projection.h>
#include <gradient_cjt.h>
#include <cmethods.h>

void qovmaFirstColM(double *a, double *b, double *c, double *d, double *x, unsigned int n, double *e);
void qovmaFirstRowM(double *a, double *b, double *c, double *d, double *x, unsigned int n, double *e);

class Problem1M : protected RnFunction, protected IGradient, public IPrinter, public IProjection
{
public:
    Problem1M();

protected:
    virtual double fx(const DoubleVector &x);
    virtual void gradient(const DoubleVector &x, DoubleVector &g);
    virtual void print(unsigned int iteration, const DoubleVector &x, const DoubleVector &gradient, double alpha, RnFunction *fn) const;
    virtual void print(const DoubleVector &x, const DoubleVector &g, unsigned int iterationNumber) const;
    virtual void project(DoubleVector &x, int index);

    void calculateU(DoubleMatrix &u);
    void calculateP(DoubleMatrix &p, const DoubleMatrix &u);

    double initial(unsigned int n) const;
    double mu(unsigned int i UNUSED_PARAM) const { return 1.0; }

private:
    unsigned int L;
    unsigned int N;
    unsigned int M;
    double hx;
    double ht;
    double h;

    double Ti;
    double Te;
    double alpha;
    double lambda0;
    double lambdal;

    double alpha0;
    double alpha1;
    double alpha2;
    double alpha3;
    double a;

    DoubleVector V;
    const DoubleVector *px;

public:
    static void Main(int argc, char* argv[]);
};

#endif // PROBLEM1M_H

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

class Problem1L2 : protected RnFunction, protected IGradient, public IPrinter, public IProjection
{
public:
    Problem1L2();

protected:
    virtual double fx(const DoubleVector &x) const;
    virtual void gradient(const DoubleVector &x, DoubleVector &g);
    virtual void print(unsigned int i, const DoubleVector &x, const DoubleVector &g, double fx) const;
    virtual void project(DoubleVector &x, int index);

    void calculateU(DoubleMatrix &u) const;
    void calculateP(DoubleMatrix &p, const DoubleMatrix &u);

    double initial(unsigned int n) const;
    double mu(unsigned int i UNUSED_PARAM) const { return 1.0; }

    double integral(const DoubleVector &x) const;
    double norm(const DoubleVector &x0) const;

    void getComponents(DoubleVector &k1, DoubleVector &z1, DoubleVector &e1, const DoubleVector &x) const;
    void qovmaFirstColM(double *a, double *b, double *c, double *d, double *x, unsigned int n, double *e) const;
    void qovmaFirstRowM(double *a, double *b, double *c, double *d, double *x, unsigned int n, double *e) const;

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
    //DoubleVector e;
    //DoubleVector z;
    //DoubleVector k;

    double zmin;
    double zmax;

public:
    static void Main(int argc, char* argv[]);
};

#endif // PROBLEM1L2_H

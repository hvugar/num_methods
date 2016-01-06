#ifndef DISCRETEHYPERBOLIC_H
#define DISCRETEHYPERBOLIC_H

#include <function.h>
#include <projection.h>
#include <printer.h>
#include <doublevector.h>
#include <math.h>

class DiscreteHyperbolic : public RnFunction, public R1Function, public Projection, public Printer
{
public:
    DiscreteHyperbolic();
    ~DiscreteHyperbolic();

    void calculateSettings();

public:
    virtual double fx(const DoubleVector& x);
    virtual void gradient(const DoubleVector& x, DoubleVector& g, double gradient_step=0.000001);

    virtual double fx(double x);

    virtual void project(DoubleVector &x, int index);

    virtual void print(unsigned int iteration, const DoubleVector& x, const DoubleVector &gradient, double alpha, RnFunction* fn) const;

    void calculateU(const DoubleVector& v, DoubleMatrix& u);
    void calculateP(const DoubleMatrix &u, DoubleVector &g);

    virtual double fi1(unsigned int i) const;
    virtual double fi2(unsigned int i) const;
    virtual double m1(unsigned int j) const;
    virtual double m2(unsigned int j) const;
    virtual double f(unsigned int i, unsigned int j) const;

    double t0;
    double t1;
    double x0;
    double x1;
    double ht;
    double hx;
    unsigned int M;
    unsigned int N;
    unsigned int D;
    unsigned int L;

    double U;
    DoubleVector e;

    double a;
    double lamda;
    double R;

    const DoubleVector *pv;
    const DoubleMatrix *pu;

    unsigned int count;

    static void main();
};

#endif // DISCRETEHYPERBOLIC_H

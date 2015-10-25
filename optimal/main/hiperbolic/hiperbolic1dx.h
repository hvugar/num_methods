#ifndef HIPERBOLIC1DX_H
#define HIPERBOLIC1DX_H

#include <function.h>
#include <printer.h>
#include <projection.h>
#include <doublevector.h>

class Hiperbolic1DX : public RnFunction, Projection, Printer
{
public:
    Hiperbolic1DX();
    ~Hiperbolic1DX();

    virtual double fx(const DoubleVector& x);
    virtual void gradient(const DoubleVector& x, DoubleVector& g, double gradient_step=0.000001);
    virtual void project(DoubleVector &x, int index);
    virtual void print(unsigned int iteration, const DoubleVector& x, const DoubleVector &gradient, double alpha, RnFunction* fn) const;

    void calculateU(const DoubleVector& e, DoubleVector& u);
    void calculateP(const DoubleVector& e, const DoubleMatrix& uT, DoubleVector& g);

    double u(double x, double t);
    double fi1(double x);
    double fi2(double x);
    double mu1(double t);
    double mu2(double t);
    double f(double x, double t);
    static void main();

private:
    double t0;
    double t1;
    double x0;
    double x1;
    double ht;
    double hx;
    unsigned int M;
    unsigned int N;
    double a;
    double lamda;
};

#endif // HIPERBOLIC1DX_H

#ifndef HYPERBOLIC1DX_H
#define HYPERBOLIC1DX_H

#include <function.h>
#include <printer.h>
#include <projection.h>
#include <doublevector.h>

class Hiperbolic1DX : public RnFunction, Projection, Printer
{
public:
    Hiperbolic1DX(unsigned int M, unsigned int N);

    virtual double fx(const DoubleVector& x);
    virtual void gradient(const DoubleVector& x, DoubleVector& g, double gradient_step=0.000001);
    virtual void project(DoubleVector &x, int index);
    virtual void print(unsigned int iteration, const DoubleVector& x, const DoubleVector &gradient, double alpha, RnFunction* fn) const;

    void calculateU(const DoubleVector& e, DoubleVector& u);
    void calculateP(const DoubleVector& e, DoubleVector& g);
    void calculateG(const DoubleVector& e, const DoubleVector& psi, DoubleVector& g, unsigned int j);
    void calculateG1(const DoubleVector& e, DoubleVector& g);

    double calculateIntegral(const DoubleVector& e);
    double calculateNorm(const DoubleVector& e);
    void psiDerivative(double &psiX, double e, const DoubleVector& psi);

    double u(double x, double t);
    double fi1(double x);
    double fi2(double x);
    double mu1(double t);
    double mu2(double t);
    double f(double x, double t);

    double fxt(unsigned int i, unsigned int j, const DoubleVector& e);

    static void main();

    double pfi1(double x) const;
    double pfi2(unsigned int i) const;
    double pmu1(double t) const;
    double pmu2(double t) const;

    double f1(double t) const { return t; }
    double f2(double t) const { return t*t; }
    double f3(double t) const { return t*t*t; }

    void initialize();
    void test(int j);

private:
    double t0;
    double t1;
    double x0;
    double x1;
    double ht;
    double hx;
    unsigned int M;
    unsigned int N;
    unsigned int L;
    double a;
    double lamda;

    DoubleVector uT;
    DoubleVector U;
};

#endif // HYPERBOLIC1DX_H

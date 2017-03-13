#ifndef HYPERBOLIC1DX_H
#define HYPERBOLIC1DX_H

#include <function.h>
#include <printer.h>
#include <projection.h>
#include <vector2d.h>
#include <matrix2d.h>
#include <hyperbolicequation.h>
#include <gradient_cjt.h>
#include <stdlib.h>
#include <math.h>

class Hyperbolic1DX : public RnFunction, public IGradient, public IHyperbolicEquation, public IBackwardHyperbolicEquation, public IProjection, public IPrinter
{
public:
    Hyperbolic1DX(unsigned int M, unsigned int N);
    virtual ~Hyperbolic1DX() {}

    virtual double fx(const DoubleVector &e) const;
    virtual void gradient(const DoubleVector &e, DoubleVector &g);
    virtual void project(DoubleVector &x, int index);
    virtual void print(unsigned int i, const DoubleVector& f0, const DoubleVector &g, double fx, GradientMethod::MethodResult result) const;

    virtual double initial1(unsigned int i) const;
    virtual double initial2(unsigned int i) const;
    virtual double boundary(Boundary type, unsigned int j) const;
    virtual double f(unsigned int i, unsigned int j) const;

    virtual double binitial1(unsigned int i) const;
    virtual double binitial2(unsigned int i) const;
    virtual double bboundary(Boundary type, unsigned int j) const;
    virtual double bf(unsigned int i, unsigned int j) const;

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
    unsigned int L;
    double a;
    double lamda;

    double v1(double t) const { return 3.0*t; }
    double v2(double t) const { return 2.5*t; }
    double v3(double t) const { return 4.0*t; }

    void calculateG(const DoubleVector& e, const DoubleVector& psi, DoubleVector& g, unsigned int j);
    void psiDerivative(double &psiX, double e, const DoubleVector& psi);

    DoubleVector U;
    const DoubleVector *pe;
    const DoubleVector *pu;
};

#endif // HYPERBOLIC1DX_H

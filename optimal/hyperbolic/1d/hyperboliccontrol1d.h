#ifndef HYPERBOLICCONTROL1D_H
#define HYPERBOLICCONTROL1D_H

#include <function.h>
#include <hyperbolicequation.h>
#include <vector2d.h>
#include <matrix2d.h>
#include <printer.h>
#include <projection.h>

class HyperbolicControl1D : public RnFunction, public IGradient, public IHyperbolicEquation, public IBackwardHyperbolicEquation, public IPrinter, public IProjection
{
public:
    HyperbolicControl1D();
    virtual ~HyperbolicControl1D() {}

    virtual double fx(const DoubleVector &v) const;
    virtual void gradient(const DoubleVector &v, DoubleVector &g);
    virtual void print(unsigned int i, const DoubleVector& f0, const DoubleVector &g, double fx) const;
    virtual void project(DoubleVector &x, int index);

    virtual double initial1(unsigned int i) const;
    virtual double initial2(unsigned int i) const;
    virtual double boundary(Boundary type, unsigned int j) const;
    virtual double f(unsigned int i, unsigned int j) const;

    virtual double binitial1(unsigned int i) const;
    virtual double binitial2(unsigned int i) const;
    virtual double bboundary(Boundary type, unsigned int j) const;
    virtual double bf(unsigned int i, unsigned int j) const;

    static void main();
protected:
    double t0;
    double t1;
    double x0;
    double x1;
    double a;
    double ht;
    double hx;
    double dt;

    unsigned int M;
    unsigned int N;
    double lamda;

    double v1(double t) const { return t*t; }
    double v2(double t) const { return t*t + 1.0; }

    const DoubleVector *pv;
    const DoubleVector *pu;
    DoubleVector U;
};

#endif

#ifndef HEATCONTROL2D_H
#define HEATCONTROL2D_H

#include <function.h>
#include <doublevector.h>
#include <parabolicequation.h>
#include <gradient_cjt.h>
#include <gradient_sd.h>
#include <printer.h>

class HeatControl2D :public RnFunction, public IGradient, public IParabolicEquation2D, public IBackwardParabolicEquation2D, public IPrinter
{
public:
    HeatControl2D(unsigned int M, unsigned int N2, unsigned int N1);
    virtual ~HeatControl2D();

    virtual double fx(const DoubleVector &x);
    virtual void gradient(const DoubleVector &x, DoubleVector &g);
    virtual void print(unsigned int i, const DoubleVector &x, const DoubleVector &s, double a, RnFunction *f) const;

    double fi(unsigned int i, unsigned int j) const;
    double m1(unsigned int j, unsigned int k) const;
    double m2(unsigned int j, unsigned int k) const;
    double m3(unsigned int i, unsigned int k) const;
    double m4(unsigned int i, unsigned int k) const;
    double f(unsigned int i, unsigned j, unsigned int k) const;

    double bfi(unsigned int i, unsigned int j) const;
    double bm1(unsigned int j, unsigned int k) const;
    double bm2(unsigned int j, unsigned int k) const;
    double bm3(unsigned int i, unsigned int k) const;
    double bm4(unsigned int i, unsigned int k) const;
    double bf(unsigned int i, unsigned int j, unsigned int k) const;

    static void main();
private:
    double u(double x1, double x2, double t) const;
    double fxt(double x1, double x2, double t);

    unsigned int N1;
    unsigned int N2;
    unsigned int M;
    double h1;
    double h2;
    double ht;

    double t0;
    double t1;
    double x10;
    double x11;
    double x20;
    double x21;

    double a1;
    double a2;

    DoubleMatrix U;
    const DoubleVector *pf;
    const DoubleMatrix *pu;
};

#endif // HEATCONTROL2D_H

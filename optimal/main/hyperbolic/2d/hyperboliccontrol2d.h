#ifndef HYPERBOLICCONTROL2D_H
#define HYPERBOLICCONTROL2D_H

#include <function.h>
#include <hyperbolicequation.h>
#include <printer.h>
#include <gradient_cjt.h>

class HyperbolicControl2D : public RnFunction, public IGradient, public IHyperbolicEquation2D, public IBackwardHyperbolicEquation2D, public IPrinter
{
public:
    HyperbolicControl2D();
    virtual ~HyperbolicControl2D();

    virtual double fx(const DoubleVector &x);
    virtual void gradient(const DoubleVector &v, DoubleVector &g);

    virtual double fi1(unsigned int i, unsigned int j) const;
    virtual double fi2(unsigned int i, unsigned int j) const;
    virtual double m1(unsigned int j, unsigned int k) const;
    virtual double m2(unsigned int j, unsigned int k) const;
    virtual double m3(unsigned int i, unsigned int k) const;
    virtual double m4(unsigned int i, unsigned int k) const;
    virtual double f(unsigned int i, unsigned int j, unsigned int k) const;

    virtual double bfi1(unsigned int i, unsigned int j) const;
    virtual double bfi2(unsigned int i, unsigned int j) const;
    virtual double bm1(unsigned int j, unsigned int k) const;
    virtual double bm2(unsigned int j, unsigned int k) const;
    virtual double bm3(unsigned int i, unsigned int k) const;
    virtual double bm4(unsigned int i, unsigned int k) const;
    virtual double bf(unsigned int i, unsigned int j, unsigned int k) const;

    virtual void print(unsigned int iteration, const DoubleVector& x, const DoubleVector &gradient, double alpha, RnFunction* fn) const;

    double v1(double t) const;
    double v2(double t) const;
    double v3(double t) const;
    double fxt(unsigned int i, unsigned int j, unsigned int k) const;

    static void main();

private:
    double t0;
    double t1;
    double x10;
    double x11;
    double x20;
    double x21;
    double ht;
    double h1;
    double h2;
    unsigned int M;
    unsigned int N1;
    unsigned int N2;
    unsigned int L;
    double a1;
    double a2;
    double alpha0;
    double alpha1;
    double qamma;

    DoubleMatrix U0;
    DoubleMatrix U1;
    DoubleVector E;
    const DoubleVector *pv;
    const DoubleCube *pu;

    inline double u(unsigned int i, unsigned int j, unsigned int k) const;
    double norm(const DoubleVector& v) const;
};

#endif // HYPERBOLICCONTROL2D_H

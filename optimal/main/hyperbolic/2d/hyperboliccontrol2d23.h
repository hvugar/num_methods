#ifndef HYPERBOLICCONTROL2D23_H
#define HYPERBOLICCONTROL2D23_H

#include <function.h>
#include <hyperbolicequation.h>
#include <printer.h>
#include <gradient_cjt.h>

class HyperbolicControl2D23 : public R1Function, public RnFunction,
        public IHyperbolicEquation2D, public IBackwardHyperbolicEquation2D,
        public IGradient, public IPrinter
{
public:
    HyperbolicControl2D23();
    virtual ~HyperbolicControl2D23() {}

    virtual double fx(double x);
    virtual double fx(const DoubleVector &x);

    virtual void gradient(const DoubleVector &v, DoubleVector &g);
    virtual void print(unsigned int iteration, const DoubleVector& x, const DoubleVector &gradient, double alpha, RnFunction* fn) const;

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

public:
    static void main(int argc, char ** argv);

    double fxt(unsigned int i, unsigned int j, unsigned int k) const;

    double t0;
    double t1;
    double x10;
    double x11;
    double x20;
    double x21;
    double h1;
    double h2;
    double ht;
    unsigned int L;
    unsigned int N1;
    unsigned int N2;
    unsigned int M;

    double alpha0;
    double alpha1;
    double alpha2;
    double alpha3;
    double qamma;
    double a1;
    double a2;

    double U0;
    double U1;

    DoubleVector e;
    DoubleVector x;
    const DoubleVector *pv;
    const DoubleCube *pu;
    DoubleVector v0;

    FILE *file;
};

#endif // HYPERBOLICCONTROL2D23_H

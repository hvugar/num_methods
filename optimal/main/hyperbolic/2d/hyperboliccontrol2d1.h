#ifndef HYPERBOLICCONTROL2D1_H
#define HYPERBOLICCONTROL2D1_H

#define ONLY_POWER
//#define ONLY_COORDINATE
//#define POWER_COORDINATE

#include <function.h>
#include <hyperbolicequation.h>
#include <printer.h>
#include <gradient_cjt.h>

class HyperbolicControl2D1 : public R1Function, public RnFunction, public IGradient, public IHyperbolicEquation2D, public IBackwardHyperbolicEquation2D, public IPrinter, public Projection
{
public:
    HyperbolicControl2D1();
    virtual ~HyperbolicControl2D1() {}

    virtual double fx(double x);
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
    virtual void project(DoubleVector &x, int index);


public:
    static void main();

    double u(double i, double j, double k) const;
    void psiDerivative(double &psiX1, double &psiX2, double x1, double x2, const DoubleMatrix &psi);
    double fxt(unsigned int i, unsigned int j, unsigned int k) const;
    double norm(const DoubleVector& v) const;
    void initialize();

    double v1(double t) const;
    double v2(double t) const;

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
    double a;

    DoubleMatrix U0;
    DoubleMatrix U1;

    DoubleVector e;
    DoubleVector c;
    const DoubleVector *px;
    const DoubleCube *pu;
    double h;

//    double vd;
//    double vu;
};

#endif // HYPERBOLICCONTROL2D1_H

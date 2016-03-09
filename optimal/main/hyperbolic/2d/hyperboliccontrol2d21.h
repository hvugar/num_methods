#ifndef HYPERBOLICCONTROL2D2_H
#define HYPERBOLICCONTROL2D2_H

#include <function.h>
#include <hyperbolicequation.h>
#include <printer.h>
#include <gradient_cjt.h>

class HyperbolicControl2D2 : public R1Function, public RnFunction, public IHyperbolicEquation2D
{
public:
    HyperbolicControl2D2();
    virtual ~HyperbolicControl2D2() {}

    virtual double fx(double x);
    virtual double fx(const DoubleVector &x);

    virtual double fi1(unsigned int i, unsigned int j) const;
    virtual double fi2(unsigned int i, unsigned int j) const;
    virtual double m1(unsigned int j, unsigned int k) const;
    virtual double m2(unsigned int j, unsigned int k) const;
    virtual double m3(unsigned int i, unsigned int k) const;
    virtual double m4(unsigned int i, unsigned int k) const;
    virtual double f(unsigned int i, unsigned int j, unsigned int k) const;

public:
    static void main();

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

    FILE *file;
};

#endif // HYPERBOLICCONTROL2D2_H

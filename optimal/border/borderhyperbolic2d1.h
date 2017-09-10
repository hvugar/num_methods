#ifndef BORDERHYPERBOLIC21_H
#define BORDERHYPERBOLIC21_H

#include <pde_old/hyperbolicequation.h>

class BorderHyperbolic21 : public IHyperbolicEquation2D
{
public:
    BorderHyperbolic21();
    virtual ~BorderHyperbolic21();

    virtual double initial1(unsigned int i, unsigned int j) const;
    virtual double initial2(unsigned int i, unsigned int j) const;
    virtual double boundary(unsigned int i, unsigned int j, unsigned int k) const;
    virtual double f(unsigned int i, unsigned int j, unsigned int k) const;

    double fxt(unsigned int i, unsigned int j, unsigned int k) const;

    double x10;
    double x11;
    double x20;
    double x21;
    double t0;
    double t1;
    double h1;
    double h2;
    double ht;
    unsigned int N1;
    unsigned int N2;
    unsigned int M;
    double a1;
    double a2;
    double qamma;
    double U;

    double alpha0;
    double alpha1;
    double alpha2;
    double alpha3;

    DoubleVector e;

    static void main();
};

#endif // BORDERHYPERBOLIC21_H

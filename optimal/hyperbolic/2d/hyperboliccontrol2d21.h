#ifndef HYPERBOLICCONTROL2D21_H
#define HYPERBOLICCONTROL2D21_H

#include <function.h>
#include <hyperbolicequation.h>
#include <printer.h>
#include <gradient_cjt.h>

class MINIMUMSHARED_EXPORT HyperbolicControl2D21 : public R1Function, public RnFunction, public IHyperbolicEquation2D
{
public:
    HyperbolicControl2D21();
    virtual ~HyperbolicControl2D21() {}

    virtual double fx(double x) const;
    virtual double fx(const DoubleVector &x) const;

    virtual double initial1(unsigned int i, unsigned int j) const;
    virtual double initial2(unsigned int i, unsigned int j) const;
    virtual double boundary(unsigned int i, unsigned int j, unsigned int k) const;
    virtual double f(unsigned int i, unsigned int j, unsigned int k) const;

public:
    static void Main(int argc, char *argv[]);

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

#endif // HYPERBOLICCONTROL2D21_H

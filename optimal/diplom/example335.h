#ifndef PARABOLICCONTROL2D335_H
#define PARABOLICCONTROL2D335_H

#include <math.h>
#include <stdlib.h>

#include <function.h>
#include <parabolicequation.h>
#include <printer.h>
#include <projection.h>
#include <tomasmethod.h>
#include <gradient_cjt.h>
#include <gradient_sd.h>

class Parabolic1DControl335 : public RnFunction, public IGradient, public IParabolicEquation2D, public IBackwardParabolicEquation2D, public IPrinter, public Projection
{
public:
    Parabolic1DControl335();
    virtual ~Parabolic1DControl335() {}

    virtual double fx(const DoubleVector &v);
    virtual void gradient(const DoubleVector &v, DoubleVector &g);
    virtual void print(unsigned int iteration, const DoubleVector &v, const DoubleVector &g, double alpha, RnFunction* fn) const;
    virtual void project(DoubleVector &v, int index);

    virtual double fi(unsigned int i, unsigned int j) const;
    virtual double m1(unsigned int j, unsigned int k) const;
    virtual double m2(unsigned int j, unsigned int k) const;
    virtual double m3(unsigned int i, unsigned int k) const;
    virtual double m4(unsigned int i, unsigned int k) const;
    virtual double f(unsigned int i, unsigned int j, unsigned int k) const;

    virtual double bfi(unsigned int i, unsigned int j) const;
    virtual double bm1(unsigned int j, unsigned int k) const;
    virtual double bm2(unsigned int j, unsigned int k) const;
    virtual double bm3(unsigned int i, unsigned int k) const;
    virtual double bm4(unsigned int i, unsigned int k) const;
    virtual double bf(unsigned int i, unsigned int j, unsigned int k) const;

    static void main(int argc, char *argv[]);
private:
    double u(double x1, double x2, double t) const { return x1*x1 + x2*x2 + t*t; }
    double norm(const DoubleVector& v) const;

    double t0;
    double t1;

    double x10;
    double x11;
    double x20;
    double x21;

    unsigned int M;
    unsigned int N1;
    unsigned int N2;
    unsigned int L;

    double a1;
    double a2;

    double h1;
    double h2;
    double ht;

    double alpha;
    double gause_a;
    double gause_b;

    DoubleMatrix U;
    const DoubleVector *pv;
    const DoubleMatrix *pu;
    DoubleVector E;
};

#endif // PARABOLICCONTROL2D335_H

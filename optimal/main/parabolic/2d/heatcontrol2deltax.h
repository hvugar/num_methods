#ifndef HEATCONTROL2DELTAX_H
#define HEATCONTROL2DELTAX_H

#include <math.h>
#include <stdlib.h>

#include <function.h>
#include <parabolicequation.h>
#include <printer.h>
#include <projection.h>
#include <tomasmethod.h>
#include <gradient_cjt.h>
#include <gradient_sd.h>

class HeatControl2DeltaX : public RnFunction, public IGradient, public IParabolicEquation2D, public IBackwardParabolicEquation2D, public IPrinter, public Projection
{
public:
    HeatControl2DeltaX(unsigned int M, unsigned int N2, unsigned int N1);
    virtual ~HeatControl2DeltaX() {}

    virtual double fx(const DoubleVector &x);
    virtual void gradient(const DoubleVector &x, DoubleVector &g);
    virtual void print(unsigned int iteration, const DoubleVector &x, const DoubleVector &g, double alpha, RnFunction* fn) const;
    virtual void project(DoubleVector &x, int index);

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

    void calculateGX(const DoubleVector &x, const DoubleMatrix& psi, DoubleVector& g, unsigned int k);

    static void main();
private:
    double u(double x1, double x2, double t) const { return x1*x1 + x2*x2 + t*t; }
    double norm(const DoubleVector& v) const;

    inline double v1(double t) const { return 10*t; }
    inline double v2(double t) const { return 20*t; }
    inline double v3(double t) const { return 30*t; }

    inline void psiDerivative(double &psiX1, double &psiX2, double e1, double e2, const DoubleMatrix &psi);

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
    const DoubleVector *px;
    const DoubleMatrix *pu;
};

#endif // HEATCONTROL2DELTAX_H

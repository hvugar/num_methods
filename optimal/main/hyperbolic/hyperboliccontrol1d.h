#ifndef HYPERBOLICCONTROL1D_H
#define HYPERBOLICCONTROL1D_H

#include <function.h>
#include <hyperbolicequation.h>
#include <doublevector.h>
#include <printer.h>
#include <projection.h>

class HyperbolicControl1D : public RnFunction, IPrinter, Projection
{
public:
    HyperbolicControl1D();
    virtual ~HyperbolicControl1D() {}

    virtual double fx(const DoubleVector& x);
    virtual void gradient(const DoubleVector& x, DoubleVector& g, double gradient_step=0.000001);

    virtual double fi1(unsigned int i) const;
    virtual double fi2(unsigned int i) const;
    virtual double m1(unsigned int j) const;
    virtual double m2(unsigned int j) const;
    virtual double f(unsigned int i, unsigned int j) const;

    double pfi1(unsigned int i) const { return 0.0; }
    double pfi2(unsigned int i) const { return 2.0*(uT[i] - U[i]); }
    double pmu1(unsigned int j) const { return 0.0; }
    double pmu2(unsigned int j) const { return 0.0; }

    virtual void print(unsigned int iteration, const DoubleVector& x, const DoubleVector &gradient, double alpha, RnFunction* fn) const;

    void project(DoubleVector &x, int index) {}

    void calculateU(const DoubleVector& v, DoubleVector &u);
    void calculareP(const DoubleVector &u, DoubleVector &g);
    void calculateG(const DoubleVector& psi, DoubleVector& g, unsigned int j);
    void initialize();

    double g1(double t) const { return t*t; }
    double g2(double t) const { return t*t + 1.0; }

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
    unsigned int DM;
    double lamda;
    DoubleVector U;
    DoubleVector uT;

    const DoubleVector *pv;
};

#endif

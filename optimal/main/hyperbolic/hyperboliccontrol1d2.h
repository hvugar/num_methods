#ifndef HYPERBOLICCONTROL1D2_H
#define HYPERBOLICCONTROL1D2_H

#include <function.h>
#include <doublevector.h>
#include <gradient_cjt.h>
#include <projection.h>
#include <printer.h>
#include <hyperbolicequation.h>
#include <tomasmethod.h>
#include <stdlib.h>

class HyperbolicControl1D2 : public RnFunction, public Printer, public Projection
{
public:
    HyperbolicControl1D2(double t0, double t1);
    virtual ~HyperbolicControl1D2();

    virtual double fx(const DoubleVector& x);
    virtual void gradient(const DoubleVector& x, DoubleVector& g, double gradient_step=0.000001);

    virtual double fi1(unsigned int i) const { return i*hx*i*hx; }
    virtual double fi2(unsigned int i) const { return 0.0; }
    virtual double m1(unsigned int j) const {  return (*pv)[j]; }
    virtual double m2(unsigned int j) const { return (*pv)[M+1+DM + j]; }
    virtual double f(unsigned int i, unsigned int j) const { return 0.0; }

    double pfi1(unsigned int x) const { return 0.0; }
    double pfi2(unsigned int i) const { return 0.0; }
    double pmu1(unsigned int t) const { return 0.0; }
    double pmu2(unsigned int t) const { return 0.0; }

    virtual void print(unsigned int iteration, const DoubleVector& x, const DoubleVector &gradient, double alpha, RnFunction* fn) const;

    void project(DoubleVector &x, int index) {}

    void calculateU(const DoubleVector& v, DoubleMatrix &u);
    void calculareP(const DoubleMatrix &u, DoubleVector &g);
    void calculateG(const DoubleVector& psi, DoubleVector& g, unsigned int j);
    void initialize();

    double g1(double t) const { return t*t; }
    double g2(double t) const { return t*t + 1.0; }

    static void main();

public:
    double t0;
    double t1;
    double x0;
    double x1;
    double a;
    double ht;
    double hx;
    //double dt;

    unsigned int M;
    unsigned int N;
    unsigned int DM;
    double lamda;
    DoubleMatrix U;
    DoubleMatrix uT;

    const DoubleVector *pv;
    double R;

    double UT;
};

#endif

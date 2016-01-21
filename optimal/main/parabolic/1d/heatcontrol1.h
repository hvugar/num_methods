#ifndef HEATCONTROL1_H
#define HEATCONTROL1_H

#include <function.h>
#include <doublevector.h>
#include <gradient_cjt.h>
#include <parabolicequation.h>
#include <printer.h>

struct HeatControl1 : public RnFunction, public IGradient, public IParabolicEquation, public IPrinter
{
    HeatControl1();
    virtual ~HeatControl1() {}

    virtual double fx(const DoubleVector& x);
    virtual void gradient(const DoubleVector& x, DoubleVector& g);

    virtual void print(unsigned int iteration, const DoubleVector& x, const DoubleVector &gradient, double alpha, RnFunction* fn) const;

    virtual double fi(unsigned int i) const;
    virtual double m1(unsigned int j) const;
    virtual double m2(unsigned int j) const;
    virtual double f(unsigned int i, unsigned int j) const;

//    virtual double pfi(unsigned int i) const { return 0.0; }
    virtual double pm1(unsigned int j) const { return 0.0; }
    virtual double pm2(unsigned int j) const { return 0.0; }
//    virtual double pf(unsigned int i, unsigned int j) const { return 0.0; }

    double fxt(double x, double t) { return 2.0*t - 2.0*a; }

    void calculateP(const DoubleVector &u, DoubleVector &g);
    void calculateG(const DoubleVector &f, const DoubleVector& psi, DoubleVector& g, unsigned int j);

    static void main();

private:
    double t0;
    double t1;
    double x0;
    double x1;
    unsigned int N;
    unsigned int M;
    double hx;
    double ht;
    double a;

    const DoubleVector *pf;
    DoubleVector U;
};

#endif // HEATCONTROL1_H

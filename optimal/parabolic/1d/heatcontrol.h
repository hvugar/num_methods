#ifndef HEATCONTROL_H
#define HEATCONTROL_H

#include <math.h>
#include <stdlib.h>
#include <function.h>
#include <parabolicequation.h>
#include <printer.h>
#include <gradient_cjt.h>
#include <gradient_sd.h>

#include <grid/pibvp.h>
#include <grid/bpibvp.h>

/**
 * @brief The HeatControl struct
 * du/dt = a(d^2u/dx^2) + f(x,t);
 * u(x,0) = fi(x);
 * u(0,t) = m1(t);
 * u(l,t) = m2(t);
 */
class MINIMUMSHARED_EXPORT HeatControl : public RnFunction, public IGradient,
        public IParabolicEquation, public IBackwardParabolicEquation,
        public IPrinter
{
public:
    HeatControl();
    virtual ~HeatControl() {}

    virtual double fx(const DoubleVector& f);
    virtual void gradient(const DoubleVector& f, DoubleVector &g);

    virtual double initial(unsigned int i) const;
    virtual double boundary(Boundary type, unsigned int j) const;
    inline virtual double f(unsigned int i, unsigned int j) const;

    virtual double binitial(unsigned int i) const;
    virtual double bboundary(Boundary type, unsigned int j) const;
    virtual double bf(unsigned int i, unsigned int j) const;

    virtual void print(unsigned int i, const DoubleVector& f0, const DoubleVector &s, double a, RnFunction* f) const;

private:
    double ht;
    double hx;
    unsigned int N;
    unsigned int M;
    double a;

    double u(double x, double t) const;
    double fxt(double x, double t);
    const DoubleVector* pf;
    const DoubleVector* pu;
    DoubleVector U;


public:
    static void Main(int argc, char *argv[]);
};

#endif // HEATCONTROL_H

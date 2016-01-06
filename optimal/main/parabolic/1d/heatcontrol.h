#ifndef HEATCONTROL_H
#define HEATCONTROL_H

#include <math.h>
#include <stdlib.h>
#include <function.h>
#include <parabolicequation.h>
#include <printer.h>
#include <gradient_cjt.h>
#include <gradient_sd.h>

/**
 * @brief The HeatControl struct
 * du/dt = a(d^2u/dx^2) + f(x,t);
 * u(x,0) = fi(x);
 * u(0,t) = m1(t);
 * u(l,t) = m2(t);
 */

struct HeatControl : public RnFunction, public Printer
{
public:
    HeatControl();
    virtual ~HeatControl() {}

    virtual double fx(const DoubleVector& u);
    virtual void gradient(const DoubleVector& f, DoubleVector &g, double gradient_step);

    unsigned int N;
    unsigned int M;
    unsigned int C;

protected:
    void initializeU();

private:
    double t0;
    double t1;

    double x0;
    double x1;

    double ht;
    double hx;

    double a1;

    inline double u(double x, double t) const { return x*x+t*t; }

    inline double fi(double x) { return u(x, 0.0); }
    inline double m1(double t) { return u(0.0, t); }
    inline double m2(double t) { return u(1.0, t); }

    inline double pfi(double x) { return 0.0; }
    inline double pm1(double t) { return 0.0; }
    inline double pm2(double t) { return 0.0; }

    inline double fxt(double x, double t) { return 2.0*t - 2.0*a1; }

    inline void calculateU(const DoubleVector& f, DoubleVector& u);
    inline void calculateU1(const DoubleVector &f);
    inline void calculateP(const DoubleVector& f, DoubleVector& g);
    inline void calculateG(const DoubleVector& f, const DoubleVector& p, DoubleVector& g, unsigned int j);

    DoubleVector uT;
    DoubleVector U;

    virtual void print(unsigned int i, const DoubleVector& f0, const DoubleVector &s, double a, RnFunction* f) const;

public:
    static void main();
};

#endif // HEATCONTROL_H

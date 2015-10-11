#ifndef HEATCONTROL_H
#define HEATCONTROL_H

#include <math.h>
#include <stdlib.h>
#include <function.h>
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

struct HeatControl : public RnFunction
{
public:
    HeatControl();

protected:
    virtual double fx(const DoubleVector& u);
    virtual void gradient(double gradient_step, const DoubleVector& f, DoubleVector &g);

private:
    double t0;
    double t1;

    double x0;
    double x1;

    double dt;
    double h;

    double a1;

    unsigned int N;
    unsigned int M;
    unsigned int C;

    double u(double x, double t) const;
//    double U(double x) const;

    double fi(double x);
    double m1(double t);
    double m2(double t);

    double pfi(double x);
    double pm1(double t);
    double pm2(double t);

    double fxt(double x, double t);

    void calculateU(const DoubleVector& f);
    void calculateP(const DoubleVector& f, DoubleVector& g);

    DoubleVector uT;
    DoubleVector U;
    DoubleVector mf;

public:
    static void main();
};

struct HeatControlPrinter : public Printer
{
    virtual void print(unsigned int i, const DoubleVector& f0, const DoubleVector &s, double a, RnFunction* f) const;
};

#endif // HEATCONTROL_H

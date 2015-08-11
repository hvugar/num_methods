#ifndef HEATCONTROL_H
#define HEATCONTROL_H

#include <math.h>
#include <function.h>
#include <printer.h>
#include <gradient_cjt.h>
#include <gradient_sd.h>

struct HeatControl : public RnFunction
{
public:
    HeatControl(double t0, double t1, double x0, double x1, double dt, double dx);

protected:
    virtual double fx(const DoubleVector& u);
    virtual void gradient(double gradient_step, const DoubleVector& u, DoubleVector &g);

private:
    double t0;
    double t1;
    double x0;
    double x1;
    double dt;
    double dx;
    unsigned int n;
    unsigned int m;

    double u(double x, double t);
    double U(double x) const;
    double F(double x, double t);
    double m1(double t);
    double m2(double t);
    double fi(double x);
    double fxt1(double x, double t);
    void calculate_u();

    DoubleVector mf;

public:
    static void main();
};

#endif // HEATCONTROL_H

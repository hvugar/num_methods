#ifndef HEATCONTROL1_H
#define HEATCONTROL1_H

#include "gradient.h"
#include "gradient_sd.h"
#include "gradient_cjt.h"
#include "gridmethod.h"

class HeatControl1;

class HeatGradientMethod1 : public SteepestDescentGradient
{
public:
    HeatControl1* heatControl;
    virtual void calculateGradient();
};

class HeatControl1 : public RnFunction, public SteepestDescentGradient
{
public:
    HeatControl1();
    ~HeatControl1();

    double u(double x, double t);

    double U(double x) const;
    double F(double x, double t);
    double m1(double t);
    double m2(double t);
    double fi(double x);

    double fxt1(double x, double t);


    void calculate_u();

    void calculate();

    virtual double fx(const DoubleVector& x) const;
    virtual void gradient(DoubleVector &g) const;

public:
    HeatGradientMethod1 gradient1;

    DoubleVector mx;
    DoubleVector mt;

    DoubleVector mu;
    DoubleVector mg;
    DoubleVector mf;

    double x0; // length start
    double x1; // length end
    double t0; // time start
    double t1; // time end

    double dx; // length grid delta
    double dt; // time grid delta

    double n;  // length grid size
    double m;  // time grid size

    friend class HeatGradientMethod;
};

#endif // HEATCONTROL1_H

#ifndef HEATCONTROL_H
#define HEATCONTROL_H

#include "gradient.h"
#include "gradient_sd.h"
#include "gradient_cjt.h"
#include "gridmethod.h"

#include <vector>

class HeatControl;

class HeatGradientMethod : public SteepestDescentGradient
{
public:
    HeatControl* heatControl;
    virtual void calculateGradient();
};

class HeatControl : public RnFunction
{
public:
    HeatControl();
    ~HeatControl();

    double u(double x, double t);
    double f(double x, double t);

    double U(double x);
    double fxt1(double x, double t);

    double fi(double x);
    double m1(double t);
    double m2(double t);

    void calculate_u();
    void calculate();
    virtual double fx(const std::vector<double>& x);

private:
    HeatGradientMethod gradient;

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

#endif // HEATCONTROL_H

#ifndef HEATCONTROL_H
#define HEATCONTROL_H

#include "gradient.h"
#include "gradient_sd.h"
#include "gradient_cjt.h"

#include <vector>

typedef std::vector<double> DoubleVector;

class HeatControl : public SteepestDescentGradient
{
public:
    HeatControl();
    ~HeatControl();

    double _y(double x, double t);
    double _f(double x, double t);
    double _u(double x, double t);

    double _fi(double x);
    double _m1(double t);
    double _m2(double t);

    double JSum();

protected:
    virtual void calculateGradient();

private:
    Gradient* gradient;

    std::vector<double> x;
    std::vector<double> t;
    std::vector<DoubleVector> u;
    //std::vector<DoubleVector> f;

    double x0; // length start
    double x1; // length end
    double t0; // time start
    double t1; // time end

    double dx; // length grid delta
    double dt; // time grid delta

    double n; // length grid size
    double m; // time grid size
};

#endif // HEATCONTROL_H

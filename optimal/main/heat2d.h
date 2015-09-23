#ifndef HEAT2D_H
#define HEAT2D_H

#include <function.h>

struct Heat2DControl : public RnFunction
{
public:
    Heat2DControl();
    ~Heat2DControl();

    virtual double fx(const DoubleVector& x) = 0;
    virtual void gradient(double step, const DoubleVector& x, DoubleVector& g) = 0;

    double f(double x1, double x2, double t);
    double fi(double x1, double x2);
    double m1(double x2, double t);
    double m2(double x2, double t);
    double m3(double x1, double t);
    double m4(double x1, double t);

private:
    double t0;
    double t1;
    double x11;
    double x12;
    double x21;
    double x22;

    double N;
    double M;
};

#endif // HEAT2D_H

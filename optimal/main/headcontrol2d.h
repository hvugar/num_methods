#ifndef HEADCONTROL2D_H
#define HEADCONTROL2D_H

#include <function.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

class HeadControl2D : public RnFunction
{
public:
    HeadControl2D();

    virtual double fx(const DoubleVector& x);
    virtual void gradient(double step, const DoubleVector& x, DoubleVector& g);

    void calculateU();

protected:
    double u(double x1, double x2, double t);
    double fi(double x1, double x2);
    double m1(double x1, double t);
    double m2(double x1, double t);
    double m3(double x2, double t);
    double m4(double x2, double t);
    double f(double x1, double x2, double t);

private:
    double t0;
    double t1;
    double x10;
    double x11;
    double x20;
    double x21;

    double a1;
    double a2;

    DoubleVector e1;
    DoubleVector e2;
    int L;

    unsigned int N1;
    unsigned int N2;
    unsigned int M;

    double h1;
    double h2;
    double ht;

    DoubleMatrix U;
};

#endif // HEADCONTROL2D_H

#ifndef HEATEQUATION_H
#define HEATEQUATION_H

#include "global.h"
#include "tomasmethod.h"
#include "printer.h"

struct MINIMUMSHARED_EXPORT HeatEquation
{
    HeatEquation();
    virtual ~HeatEquation();

    void setTimeInterval(double t0, double t1);
    void setLengthInterval(double x0, double x1);
    void setPartNumber(unsigned int M, unsigned int N);

    virtual double fi(double x, unsigned int i) const = 0;
    virtual double m1(double t, unsigned int j) const = 0;
    virtual double m2(double t, unsigned int j) const = 0;
    virtual double f(double x, unsigned int i, double t, unsigned int j) const = 0;

    void calculate_u(DoubleVector& u);
    void calculate_u1(DoubleVector& u);

    /* number of parts of time */
    unsigned int M;
    /* number of parts of length */
    unsigned int N;
    unsigned int C;

    double t0;
    double t1;
    double x0;
    double x1;
    double a;
    double ht;
    double hx;
};

struct MINIMUMSHARED_EXPORT HeatEquation2D
{
    HeatEquation2D();
    virtual double fi(double x, double y) const = 0;
    virtual double m1(double y, double t) const = 0;
    virtual double m2(double y, double t) const = 0;
    virtual double m3(double x, double t) const = 0;
    virtual double m4(double x, double t) const = 0;
    virtual double f(double x, double y, double t) const = 0;

    void setBorders(double t0, double t1, double x10, double x11, double x20, double x21);
    void setPartNumbers(unsigned int N1, unsigned int N2, unsigned M);

    void calculate(DoubleMatrix& u);
    void calculateBack(DoubleMatrix& u);

    unsigned int N1;
    unsigned int N2;
    unsigned int M;

    double t0;
    double t1;
    double x10;
    double x11;
    double x20;
    double x21;
    double a1;
    double a2;
    double h1;
    double h2;
    double ht;
};

#endif // HEATEQUATION_H

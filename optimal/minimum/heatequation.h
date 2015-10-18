#ifndef HEATEQUATION_H
#define HEATEQUATION_H

#include "global.h"
#include "tomasmethod.h"

class MINIMUMSHARED_EXPORT HeatEquation
{
public:
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

    /* number of parts of time */
    unsigned int M;
    /* number of parts of length */
    unsigned int N;
    unsigned int C;

public:
    double t0;
    double t1;
    double x0;
    double x1;
    double a;
    double ht;
    double hx;
};

#endif // HEATEQUATION_H

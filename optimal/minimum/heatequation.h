#ifndef HEATEQUATION_H
#define HEATEQUATION_H

#include "tomasmethod.h"

class HeatEquation
{
public:
    HeatEquation();
    ~HeatEquation();

    void setTimeInterval(double t0, double t1);
    void setLengthInterval(double x0, double x1);
    void setPartNumber(unsigned int M, unsigned int N);

    virtual double fi(unsigned int i) const = 0;
    virtual double m1(unsigned int j) const = 0;
    virtual double m2(unsigned int j) const = 0;
    virtual double f(unsigned int j, unsigned int i) = 0;

    void calculateU(DoubleVector& u);

private:
    double t0;
    double t1;
    double x0;
    double x1;
    double a;
    double ht;
    double hx;
    /* number of parts of time */
    unsigned int M;
    /* number of parts of length */
    unsigned int N;
};

#endif // HEATEQUATION_H

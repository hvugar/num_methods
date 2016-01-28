#ifndef HEATDELTACENTER_H
#define HEATDELTACENTER_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <parabolicequation.h>
#include <printer.h>

class HeatDeltaCenter : public IParabolicEquation2D
{
public:
    HeatDeltaCenter();
    virtual ~HeatDeltaCenter() {}

    virtual double fi(unsigned int i, unsigned int j) const;
    virtual double m1(unsigned int j, unsigned int k) const;
    virtual double m2(unsigned int j, unsigned int k) const;
    virtual double m3(unsigned int i, unsigned int k) const;
    virtual double m4(unsigned int i, unsigned int k) const;
    virtual double f(unsigned int i, unsigned int j, unsigned int k) const;

    double v1(double t) const;

    static void main();
private:
    double t0;
    double t1;
    double x0;
    double x1;
    double y0;
    double y1;
    double ht;
    double hx;
    double hy;
    unsigned int M;
    unsigned int N1;
    unsigned int N2;
};

#endif // HEATDELTACENTER_H

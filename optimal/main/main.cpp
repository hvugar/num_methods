#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <function.h>
#include <r1minimize.h>
#include <gradient.h>
#include <gradient_sd.h>
#include <gradient_cjt.h>
#include <gridmethod.h>
#include <rungekutta.h>
#include <doublevector.h>
#include <heatequation.h>

#include "cfunction1.h"
#include "cfunction2.h"
#include "cfunction3.h"
#include "cfunction.h"
#include "rosenbrock.h"
#include "bealesfunction.h"
#include "boothfunction.h"
#include "pointcontrol.h"
#include "pointcontrol1.h"
#include "pointcontrol2.h"
#include "utils.h"

#include "heat/1d/heatcontrol.h"
#include "heat/1d/heatcontroldeltaf.h"
#include "heat/1d/heatcontroldeltax.h"
#include "heat/2d/heatcontrol2d.h"
#include "heat/2d/heatcontrol2delta.h"
#include "heat/2d/heatcontrol2deltaf.h"
#include "heat/2d/heatcontrol2deltax.h"

struct HeatEquation1 : public HeatEquation
{
    double u(double x, double t) const;
    virtual double fi(double x, unsigned int i) const;
    virtual double m1(double t, unsigned int j) const;
    virtual double m2(double t, unsigned int j) const;
    virtual double f(double x, unsigned int i, double t, unsigned int j) const;
};

double HeatEquation1::u(double x, double t) const
{
    return x*x + t*t;
}

double HeatEquation1::fi(double x, unsigned int i) const
{
    return u(x, 1.0);
}

double HeatEquation1::m1(double t, unsigned int j) const
{
    return u(0.0, t);
}

double HeatEquation1::m2(double t, unsigned int j) const
{
    return u(1.0, t);
}

double HeatEquation1::f(double x, unsigned int i, double t, unsigned int j) const
{
    return 2.0 * t - 2.0*a;
}


struct HeatEquation2D1 : public HeatEquation2D
{
    double u(double x1, double x2, double t) const { return 10*x1*x1 + 20*x2*x2 + 10*t*t + 10*x1*x2; }

    virtual double fi(double i, double j) const { return u(i*h1, j*h2, t1); }
    virtual double m1(double j, double k) const { return u(x10, j*h2, k*ht); }
    virtual double m2(double j, double k) const { return u(x11, j*h2, k*ht); }
    virtual double m3(double i, double k) const { return u(i*h1, x20, k*ht); }
    virtual double m4(double i, double k) const { return u(i*h1, x21, k*ht); }
    virtual double f(double i, double j, double k) const { return 20.0 * (k*ht) - 20.0*a1 - 40.0*a2; }
};

int main()
{
    HeatControl2Delta::main();

//    DoubleVector u;
//    HeatEquation1 he;
//    printf("%f\n", he.a);
//    he.setTimeInterval(0.0, 1.0);
//    he.setLengthInterval(0.0, 1.0);
//    he.setPartNumber(1000, 100);
//    he.calculate_u1(u);
//    Printer::printVector(u, he.N/10);

//    DoubleMatrix m;
//    HeatEquation2D1 he2;
//    he2.setBorders(0.0, 1.0, 0.0, 1.0, 0.0, 1.0);
//    he2.setPartNumbers(100, 100, 1000);
//    he2.calculateBack(m);
//    Printer::printMatrix(m, he2.N2/10, he2.N1/10);

    return 0;
}




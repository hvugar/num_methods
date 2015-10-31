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
#include <parabolicequation.h>
#include <hyperbolicequation.h>

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
#include "heat/1d/heatcontrol1.h"
#include "heat/1d/heatcontroldeltaf.h"
#include "heat/1d/heatcontroldeltax.h"
#include "heat/2d/heatcontrol2d.h"
#include "heat/2d/heatcontrol2delta.h"
#include "heat/2d/heatcontrol2deltaf.h"
#include "heat/2d/heatcontrol2deltax.h"

#include "hyperbolic/hyperbolic1dx.h"

struct HyperbolicEquation2D1 : public HyperbolicEquation2D
{
    HyperbolicEquation2D1(double t0, double t1, double x10, double x11, double x20, double x21, double a1, double a2, unsigned int M, unsigned int N1, unsigned int N2);
    virtual ~HyperbolicEquation2D1();

    virtual double fi1(unsigned int i, unsigned int j) const { return u(i*h1, j*h2, 0.0); }
    virtual double fi2(unsigned int i, unsigned int j) const { return 0.0; }
    virtual double m1(unsigned int j, unsigned int k) const { return u(x10, j*h2, k*ht); }
    virtual double m2(unsigned int j, unsigned int k) const { return u(x11, j*h2, k*ht); }
    virtual double m3(unsigned int i, unsigned int k) const { return u(i*h1, x20, k*ht); }
    virtual double m4(unsigned int i, unsigned int k) const { return u(i*h1, x21, k*ht); }
    virtual double f(unsigned int i, unsigned int j, unsigned int k) const { return 0.0; }

    double u(double x1, double x2, double t) const { return x1*x1*x1 + x2*x2*x2 + t*t*t; }
};

int main()
{
//    HeatControl2DeltaX::main();
//    Hyperbolic1DX::main();
//    HeatControl1 hc;
//    DoubleVector u;
//    hc.calculateU(u);
//    Printer::printVector(u);

    return 0;
}




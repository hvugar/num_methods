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

struct HyperbolicEquation1 : public HyperbolicEquation
{
    HyperbolicEquation1(double t0=0.0, double t1=1.0, double x0=0.0, double x1=1.0, unsigned int M=1000, unsigned int N=1000, double a=1.0)
        : HyperbolicEquation(t0, t1, x0, x1, M, N, a) {}
    virtual ~HyperbolicEquation1() {}

    virtual double fi1(unsigned int i) const { return u(i*hx, 0.0); }
    virtual double fi2(unsigned int i) const { return 0.0; }
    virtual double m1(unsigned int j) const { return u(x0, j*ht); }
    virtual double m2(unsigned int j) const { return u(x1, j*ht); }
    virtual double f(unsigned int i, unsigned int j) const { return 0.0; }

    double u(double x, double t) const { return x*x + t*t; }
};

int main()
{
//    HeatControl2DeltaX::main();
//    Hyperbolic1DX::main();
//    HeatControl1 hc;
    DoubleVector u;
//    hc.calculateU(u);
//    Printer::printVector(u);

    HyperbolicEquation1 he;
    he.calculate(u);
    Printer::printVector(u);

    return 0;
}




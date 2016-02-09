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
#include <r1minimize.h>

#include <iostream>
#include <stdexcept>

#include "control/cfunction.h"
#include "control/cfunction1.h"
#include "control/cfunction2.h"
#include "control/cfunction3.h"

#include "rnfunction/rosenbrock.h"
#include "rnfunction/bealesfunction.h"
#include "rnfunction/boothfunction.h"

#include "parabolic/1d/heatcontrol.h"
#include "parabolic/1d/heatcontroldeltaf.h"
#include "parabolic/1d/heatcontroldeltax.h"

#include "parabolic/2d/heatcontrol2d.h"
#include "parabolic/2d/heatcontrol2delta.h"
#include "parabolic/2d/heatcontrol2deltaf.h"
#include "parabolic/2d/heatcontrol2deltax.h"

#include "hyperbolic/1d/hyperbolic1dx.h"
#include "hyperbolic/1d/hyperboliccontrol1d.h"
#include "hyperbolic/1d/hyperboliccontrolx.h"
#include "hyperbolic/1d/hyperboliccontrolh.h"

#include "point/pointcontrol11.h"
#include "point/pointcontrol.h"
#include "point/pointcontrol1.h"
#include "point/pointcontrol2.h"

#include "discrete/discreteheat.h"
#include "discrete/discretehyperbolic.h"
#include "discrete/discretehyperbolic1.h"

class A : public IHyperbolicEquation2D
{
public:
    virtual double fi1(unsigned int i, unsigned int j) const { return (i*hx1)*(i*hx1) + (j*hx2)*(j*hx2); }
    virtual double fi2(unsigned int i, unsigned int j) const { return 0.0; }
    virtual double m1(unsigned int j, unsigned int k) const { return (j*hx2)*(j*hx2) + (k*ht)*(k*ht); }
    virtual double m2(unsigned int j, unsigned int k) const { return 1.0 + (j*hx2)*(j*hx2) + (k*ht)*(k*ht); }
    virtual double m3(unsigned int i, unsigned int k) const { return (i*hx1)*(i*hx1) + (k*ht)*(k*ht); }
    virtual double m4(unsigned int i, unsigned int k) const { return 1.0 + (i*hx1)*(i*hx1) + (k*ht)*(k*ht); }
    virtual double f(unsigned int i, unsigned int j, unsigned int k) const { return 2.0 - 2.0*a1 - 2.0*a2; }

    double ht;
    double hx1;
    double hx2;
    unsigned int M;
    unsigned int N1;
    unsigned int N2;
    double t0;
    double t1;
    double x10;
    double x11;
    double x20;
    double x21;
    double a1;
    double a2;
};

class P1 : public IParabolicEquation2D
{
public:
    virtual double fi(unsigned int i, unsigned int j) const { return (i*hx1)*(i*hx1) + (j*hx2)*(j*hx2); }
    virtual double m1(unsigned int j, unsigned int k) const { return (j*hx2)*(j*hx2) + (k*ht)*(k*ht); }
    virtual double m2(unsigned int j, unsigned int k) const { return 1.0 + (j*hx2)*(j*hx2) + (k*ht)*(k*ht); }
    virtual double m3(unsigned int i, unsigned int k) const { return (i*hx1)*(i*hx1) + (k*ht)*(k*ht); }
    virtual double m4(unsigned int i, unsigned int k) const { return 1.0 + (i*hx1)*(i*hx1) + (k*ht)*(k*ht); }
    virtual double f(unsigned int i, unsigned int j, unsigned int k) const
    {
        return 2.0*(0.5*k*ht) - 2.0*a1 - 2.0*a2;
    }

    double ht;
    double hx1;
    double hx2;
    unsigned int M;
    unsigned int N1;
    unsigned int N2;
    double t0;
    double t1;
    double x10;
    double x11;
    double x20;
    double x21;
    double a1;
    double a2;
};

int main()
{
    P1 a;
    a.x10 = a.x20 = a.t0 = 0.0;
    a.x11 = a.x21 = a.t1 = 1.0;
    a.M = a.N2 = a.N1 = 100;
    a.ht = a.hx1 = a.hx2 = 0.01;
    a.a1 = a.a2 = 1.0;
    DoubleMatrix u;
    a.calculateU(u, a.hx1, a.hx2, a.ht, a.N1, a.N2, a.M, a.a1, a.a2);
    IPrinter::printMatrix(u);
    puts("---");
    DoubleMatrix u1;
    a.calculateU1(u1, a.hx1, a.hx2, a.ht, a.N1, a.N2, a.M, a.a1, a.a2);
    IPrinter::printMatrix(u1);


    //    A a;
    //    DoubleMatrix u;
    //    a.calculateN(u, a.hx, a.ht, a.N, a.M);
    //    IPrinter::printMatrix(u);
    //    HeatControl2DeltaX::main();
    //    HeatControlDeltaX::main();
    //    DiscreteHyperbolic1::main();
    //    HyperbolicControlH::main();
    return 0;
}

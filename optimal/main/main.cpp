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

class H1 : public IHyperbolicEquation2D
{
public:
    double u(unsigned int i, unsigned int j, unsigned int k) const
    {
        double x1 = i*hx1;
        double x2 = j*hx2;
        double t  = k*ht;
        return x1*x1*x1 + x2*x2*x2 + t*t*t;
    }

    virtual double fi1(unsigned int i, unsigned int j) const { return u(j, j, 0); }
    virtual double fi2(unsigned int i, unsigned int j) const { return 0.0; }
    virtual double m1(unsigned int j, unsigned int k) const { return u(0, j, k); }
    virtual double m2(unsigned int j, unsigned int k) const { return u(N1, j, k); }
    virtual double m3(unsigned int i, unsigned int k) const { return u(i, 0, k); }
    virtual double m4(unsigned int i, unsigned int k) const { return u(i, N2, k); }
    virtual double f(unsigned int i, unsigned int j, unsigned int k) const
    {
        double x1 = i*hx1;
        double x2 = j*hx2;
        double t  = k*ht;
        return 6.0*t - 6.0*x1*a1*a1 - 6.0*x2*a2*a2;
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

class P1 : public IParabolicEquation2D
{
public:
    double u(unsigned int i, unsigned int j, unsigned int k) const
    {
        double x1 = i*hx1;
        double x2 = j*hx2;
        double t  = 0.5*k*ht;
        return x1*x1 + x2*x2 + t*t;
    }

    virtual double fi(unsigned int i, unsigned int j) const { return u(j, j, 0); }
    virtual double m1(unsigned int j, unsigned int k) const { return u(0, j, k); }
    virtual double m2(unsigned int j, unsigned int k) const { return u(N1, j, k); }
    virtual double m3(unsigned int i, unsigned int k) const { return u(i, 0, k); }
    virtual double m4(unsigned int i, unsigned int k) const { return u(i, N2, k); }
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
    H1 a;
    a.x10 = a.x20 = a.t0 = 0.0;

    a.x11 = 1.0;
    a.x21 = 1.0;
    a.t1  = 1.0;

    a.N1 = 100;
    a.N2 = 100;
    a.M  = 100;

    a.hx1 = 0.01;
    a.hx2 = 0.01;
    a.ht  = 0.01;
    a.a1 = a.a2 = 1.0;

    DoubleMatrix u;
    a.calculateU(u, a.hx1, a.hx2, a.ht, a.N1, a.N2, a.M, a.a1, a.a2);
    IPrinter::printMatrix(u);

    //    puts("---");
    //    DoubleCube c;
    //    a.calculateU(c, a.hx1, a.hx2, a.ht, a.N1, a.N2, a.M, a.a1, a.a2);
    //    IPrinter::printMatrix(c[c.size()-1]);



    //    A a;
    //    DoubleMatrix u;
    //    a.calculateN(u, a.hx, a.ht, a.N, a.M);
    //    IPrinter::printMatrix(u);
    //        HeatControl2DeltaX::main();
    //    HeatControlDeltaX::main();
    //    DiscreteHyperbolic1::main();
    //    HyperbolicControlH::main();
    return 0;
}

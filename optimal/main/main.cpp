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
#include "parabolic/1d/heatcontrol1.h"
#include "parabolic/1d/heatcontroldeltaf.h"
#include "parabolic/1d/heatcontroldeltax.h"
#include "parabolic/2d/heatcontrol2d.h"
#include "parabolic/2d/heatcontrol2delta.h"
#include "parabolic/2d/heatcontrol2deltaf.h"
#include "parabolic/2d/heatcontrol2deltax.h"

#include "hyperbolic/hyperbolic1dx.h"
#include "hyperbolic/hyperboliccontrol1d4.h"
#include "hyperbolic/hyperboliccontrolx.h"

#include "point/pointcontrol11.h"
#include "point/pointcontrol.h"
#include "point/pointcontrol1.h"
#include "point/pointcontrol2.h"

#include "discrete/discreteheat.h"
#include "discrete/discretehyperbolic.h"
#include "discrete/discretehyperbolic1.h"

struct A : public IHyperbolicEquation
{
    double fi1(unsigned int i) const
    {
        double x = i*hx;
        return x*x*x;
    }

    double fi2(unsigned int i) const
    {
        return 0.0;
    }

    double m1(unsigned int j) const
    {
        double t = j*ht;
        return t*t*t;
    }

    double m2(unsigned int j) const
    {
        double t = j*ht;
        return t*t*t+1.0;
    }

    double f(unsigned int i, unsigned int j) const
    {
        double x = i*hx;
        double t = j*ht;
        return 6.0*t - 6.0*x*a*a;
    }

    unsigned int N;
    unsigned int M;
    double hx;
    double ht;
    double a;
};


struct B : public IBackwardHyperbolicEquation
{
    double bfi1(unsigned int i) const
    {
        double x = i*hx;
        return x*x*x+0.9*0.9*0.9;
    }

    double bfi2(unsigned int i) const
    {
        return 3.0*0.9*0.9;
    }

    double bm1(unsigned int j) const
    {
        double t = j*ht;
        return t*t*t;
    }

    double bm2(unsigned int j) const
    {
        double t = j*ht;
        return t*t*t+1.0;
    }

    double bf(unsigned int i, unsigned int j) const
    {
        double x = i*hx;
        double t = j*ht;
        return 6.0*t - 6.0*x*a*a;
    }

    unsigned int N;
    unsigned int M;
    double hx;
    double ht;
    double a;
    double t0;
};

int main()
{
    B a;
    a.a = -1.0;
    a.N = 100;
    a.M = 90;
    a.hx = 1.0/a.N;
    a.ht = 0.9/a.M;
    //a.t0 = 0.1;
    DoubleMatrix u;
    a.calculateU(u, a.hx, a.ht, a.M, a.N, a.a);
    IPrinter::printMatrix(u);
    DoubleVector v;
    a.calculateU(v, a.hx, a.ht, a.M, a.N, a.a);
    IPrinter::printVector(v);

//    DiscreteHeat::main();
//    DiscreteHyperbolic1::main();
//    HyperbolicControl1D4::main();
//    HyperbolicControlX::main();
//    HeatControl2D::main();
//    Rosenbrock::main();
    return 0;
}




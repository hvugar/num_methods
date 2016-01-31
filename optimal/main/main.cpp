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

#include "hyperbolic/hyperbolic1dx.h"
#include "hyperbolic/hyperboliccontrol1d.h"
#include "hyperbolic/hyperboliccontrolx.h"
#include "hyperbolic/hyperboliccontrolh.h"

#include "point/pointcontrol11.h"
#include "point/pointcontrol.h"
#include "point/pointcontrol1.h"
#include "point/pointcontrol2.h"

#include "discrete/discreteheat.h"
#include "discrete/discretehyperbolic.h"
#include "discrete/discretehyperbolic1.h"

struct A : public IParabolicEquation
{
    A() : hx(0.01), ht(0.01), N(100), M(100) { }
    virtual double fi(unsigned int i) const { return hx*i*hx*i; }
    virtual double m1(unsigned int j) const { return 0.0; }
    virtual double m2(unsigned int j) const { return 2.0; }
    virtual double f(unsigned int i, unsigned int j) const { return 2.0*j*ht - 2.0; }

    double hx;
    double ht;
    unsigned int N;
    unsigned int M;
};

int main()
{
//    A a;
//    DoubleMatrix u;
//    a.calculateN(u, a.hx, a.ht, a.N, a.M);
//    IPrinter::printMatrix(u);
//    HeatControl2DeltaX::main();
//    DiscreteHyperbolic1::main();
    HyperbolicControlH::main();
    return 0;
}

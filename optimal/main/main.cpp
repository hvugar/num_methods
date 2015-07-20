#include <stdio.h>
#include <function.h>
#include <r1minimize.h>
#include <gradient.h>
#include <gradient_sd.h>
#include <gradient_cjt.h>
#include <gradient_prj.h>
#include <methods.h>
#include <gridmethod.h>

#include "cfunction1.h"
#include "cfunction2.h"
#include "samplecontrol.h"
#include "rosenbrock.h"
#include "heatcontrol.h"

struct Function1 : public RnFunction
{
    virtual double fx(const std::vector<double> &x);
};

double Function1::fx(const std::vector<double> &x)
{
    return x[0]*x[0] + 2.0*x[1]*x[1]*x[1]*x[1] + x[2]*x[2];
}

struct M1 : public R1Function { virtual double fx(double t) { return t*t; } };
struct M2 : public R1Function { virtual double fx(double t) { return t*t+3.0; } };
struct Fi : public R1Function { virtual double fx(double x) { return x*x+2*x; } };
struct F  : public R2Function { virtual double fx(double x, double t) { return 2.0*t - 2.0; } };

int main()
{
//    Rosenbrock::Main();
//    SampleControl::Main();
    HeatControl hc;
    hc.x0 = 0.0;
    hc.x1 = 1.0;
    hc.t0 = 0.0;
    hc.t1 = 1.0;
    hc.dx = 0.001;
    hc.dt = 0.001;
    hc.n = 1001;
    hc.m = 1001;

    hc.calculate_u();

//    GridMethod gm;
//    gm.setLengthInterval(0.0, 1.0);
//    gm.setTimeInterval(0.0, 1.0);
//    gm.setLengthTimeStep(0.001, 0.001);
//    gm.setLengthTimeStepCount(1001, 1001);
//    gm.setM1(new M1);
//    gm.setM2(new M2);
//    gm.setFi(new Fi);
//    gm.setF(new F);
//    gm.implicitDifferenceScheme();

    return 0;
}


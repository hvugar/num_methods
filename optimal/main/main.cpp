#include <stdio.h>
#include <function.h>
#include <r1minimize.h>
#include <gradient.h>
#include <gradient_sd.h>
#include <gradient_cjt.h>
#include <gridmethod.h>
#include <rungekutta.h>

#include "cfunction1.h"
#include "cfunction2.h"
#include "cfunction.h"
#include "rosenbrock.h"
#include "bealesfunction.h"
#include "boothfunction.h"
#include "heatcontrol.h"
#include "pointcontrol.h"
#include "pointcontrol1.h"
#include "utils.h"

double dx(double t, double x)
{
    double sum = 2.0*t - x + t*t + 1.0;
    if (fabs(t-0.2) < 0.00000001)
    {
        printf("%20f\n", t);
        sum += 0.2;
    }
    return sum;
}

struct A : public R1Function
{
    virtual double fx(double p)
    {
        return fabs(p - 1.08986579);
    }
};

int main()
{
    double dt = 0.0000001;
    double t0 = 0.0;
    double t1 = 0.0+dt;
    double x0 = 1.0;
    double x1 = 0.0;
    unsigned int n = (unsigned int)(ceil((t1 - t0)/dt))+1;
    DoubleVector x(n);

    x[0] = x0;
    for (unsigned int i=1; i<n; i++)
    {
        double k1 = dx(t0,        x0);
        double k2 = dx(t0+dt/2.0, x0+(dt/2.0)*k1);
        double k3 = dx(t0+dt/2.0, x0+(dt/2.0)*k2);
        double k4 = dx(t0+dt,     x0+dt*k3);
        x0 = x0 + (dt/6.0) * (k1 + 2.0*k2 + 2.0*k3 + k4);
        t0 = t0 + dt;
        x[i] = x0;
    }

    printf("%.20f %.20f\n", x[0], x[n-1]);
    //printX("x", x);

//    double alpha0 = 0.0;
//    double a,b,alpha;
//    R1Minimize::StranghLineSearch(alpha0, min_step, a, b, &r1X);
//    R1Minimize::GoldenSectionSearch(a, b, alpha, &r1X, min_epsilon);
//    if (r1X.fx(alpha) > r1X.fx(alpha0)) alpha = alpha0;
//    return alpha;

    //    PointControl1::main();
    //    PointControl::main();
//        Rosenbrock::main();
    //    BealesFunction::main();
    //    BoothFunction::main();
    //    puts("*****************************************************************");
    //    CFunction1::main();
    //    CFunction2::main();
    //    ControlFunction::main();
    //    HeatControl::main();

    return 0;
}


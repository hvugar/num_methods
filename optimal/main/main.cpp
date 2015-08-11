#include <stdio.h>
#include <function.h>
#include <r1minimize.h>
#include <gradient.h>
#include <gradient_sd.h>
#include <gradient_cjt.h>
#include <gridmethod.h>

#include "cfunction1.h"
#include "cfunction2.h"
#include "cfunction.h"
#include "rosenbrock.h"
#include "bealesfunction.h"
#include "boothfunction.h"
#include "heatcontrol.h"
#include "pointcontrol.h"

struct R1Function1 : public R1Function
{
    virtual double fx(double x);
};

double R1Function1::fx(double x)
{
    return 2.0*x*x - 12.0*x;
}



int main()
{
    PointControl::main();
//    Rosenbrock::main();
//    BealesFunction::main();
//    BoothFunction::main();
//    puts("*****************************************************************");
//    CFunction1::main();
//    CFunction2::main();
//    ControlFunction::main();
//    HeatControl::main();

    return 0;
}


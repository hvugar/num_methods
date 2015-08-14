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
#include "pointcontrol2.h"
#include "utils.h"

typedef double (*TFunc1)(double, double);

struct TFunc2 //: public TFunc1
{
    static double x;
    double y;
    double func2(double a, double b)
    {
        return a + b + x + y;
    }

    static double func3(double a, double b)
    {
        return x+a+b;
    }
};

double TFunc2::x = 10;

double func1(double x, double y)
{
     return x + y;
}

void a(TFunc1 f, double a, double b)
{
    printf("%f\n", f(a, b));
}

int main()
{
    TFunc2 f2;
    f2.x = 10;
    f2.y = 20;
    TFunc1 f1 = func1;

    a(TFunc2::func3, 1.0, 2.0);
    //    PointControl::main();
    //    PointControl1::main();
    //    PointControl2::main();
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


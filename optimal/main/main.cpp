#include <stdio.h>
#include <function.h>
#include <r1minimize.h>
#include <gradient.h>
#include <gradient_sd.h>
#include <gradient_cjt.h>
#include <gradient_prj.h>
#include <gridmethod.h>

#include "cfunction1.h"
#include "cfunction2.h"
#include "cfunction.h"
#include "rosenbrock.h"
#include "bealesfunction.h"

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
//    Rosenbrock::main();
    BealesFunction::main();
//    puts("*****************************************************************");
//    CFunction1::main();
//    CFunction2::main();
//    ControlFunction::main();

    R1Function1 f;
    double x = 5.0;
    double a = NAN;
    double b = NAN;
    double c = NAN;

    double step = 5.0;
    double epsilon = 0.15;

    unsigned int n = 9;
    double step1 = 0.2;

    R1Minimize::Swann(x, step, a, b, &f);
    printf("[%.6f, %.6f]\n", a, b);

    R1Minimize::FibonachiMethod(a, b, c, step1, epsilon, &f);
    //printf("[%.6f, %.6f] %.6f\n", a, b, (a+b)/2.0);

    return 0;
}


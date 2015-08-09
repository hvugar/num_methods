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

struct R1Function1 : public R1Function
{
    virtual double fx(double x);
};

double R1Function1::fx(double x)
{
    return (x-5)*(x-5);
}

int main()
{
//    Rosenbrock::main();
//    puts("*****************************************************************");
//    CFunction1::main();
//    CFunction2::main();
//    ControlFunction::main();

    R1Function1 f;
    double a, b;
    a = b = 0.0;
    R1Minimize::Swann(1.0, 1.0, a, b, &f);
    printf("%f %f\n", a, b);

    return 0;
}


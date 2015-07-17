#include <stdio.h>
#include <function.h>
#include <r1minimize.h>
#include <gradient.h>
#include <gradient_sd.h>
#include <gradient_cjt.h>
#include <gradient_prj.h>
#include <methods.h>

#include "cfunction1.h"
#include "cfunction2.h"
#include "samplecontrol.h"

#include "rosenbrock.h"

struct Function1 : public RnFunction
{
    virtual double fx(const std::vector<double> &x);
};

double Function1::fx(const std::vector<double> &x)
{
    return x[0]*x[0] + 2.0*x[1]*x[1]*x[1]*x[1] + x[2]*x[2];
}

int main()
{
    Rosenbrock::Main();
//    SampleControl::Main();
    return 0;
}


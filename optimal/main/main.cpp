#include <stdio.h>
#include <function.h>
#include <r1minimize.h>
#include <gradient.h>
#include <fpgradient.h>
#include <methods.h>

struct Rosenbrock : public RnFunction
{
    virtual double fx(std::vector<double> x);
};

double Rosenbrock::fx(std::vector<double> x)
{
    double x1 = x[0];
    double x2 = x[1];
    return ((1 - x1) * (1 - x1)) + 100 * (x2 - x1 * x1) * (x2 - x1 * x1);
}

int main()
{
    Rosenbrock r;
    std::vector<double> x;
    x.push_back(-1.0);
    x.push_back(+1.2);

    FastProximalGradient fg;
//    fg.setF(&r);
//    fg.setEpsilon(0.000001);
//    fg.setX(x);
//    fg.calculate();
}


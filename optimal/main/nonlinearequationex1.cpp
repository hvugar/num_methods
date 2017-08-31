#include "nonlinearequationex1.h"
#include <math.h>

void NonLinearEquationEx1::Main(int agrc UNUSED_PARAM, char *argv[] UNUSED_PARAM)
{
    NonLinearEquationEx1 nlfs;

    //    DoubleVector x0;
    //    x0 << 3.5 << 2.2;
    //    DoubleVector xs;
    //    nlfs.calculateSimpleIdetartion(x0, xs, 0.001);
    //    IPrinter::print(xs, xs.size());

    DoubleVector x1;
    x1 << 2.0 << 2.0;
    IPrinter::print(x1, x1.size(), 10, 4);
    DoubleVector rx;
    nlfs.calculateNewtonMethodMod2(x1, rx, 0.0001, 0.0001);
    IPrinter::print(rx, rx.size(), 10, 4);

}

double NonLinearEquationEx1::fx(const DoubleVector &x, unsigned int num) const
{
    double x1 = x[0];
    double x2 = x[1];
    //    if (num == 0) return sqrt(x1 + 3.0*log10(x1));
    //    if (num == 1) return sqrt((x1*(x2+5.0)-1)/2.0);

    if (num == 0) return x1+3.0*log10(x1)-x2*x2;
    if (num == 1) return 2.0*x1*x1 - x1*x2 - 5.0*x1 + 1.0;

    return NAN;
}


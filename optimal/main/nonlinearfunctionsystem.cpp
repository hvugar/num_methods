#include "nonlinearfunctionsystem.h"
#include <math.h>

//void NonLinearFunctionSystem::Main(int agrc UNUSED_PARAM, char *argv[] UNUSED_PARAM)
//{
//    NonLinearFunctionSystem nlfs;
//    DoubleVector x0;
//    x0 << 3.5 << 2.2;
//    DoubleVector xs;
//    nlfs.calculate(x0, xs);
//}

//double NonLinearFunctionSystem::fx(const DoubleVector &x, unsigned int num) const
//{
//    double x1 = x[0];
//    double x2 = x[1];
//    if (num == 0) return sqrt((x1*(x2+5.0)-1)/2.0);
//    if (num == 1) return sqrt(x1 + 3.0*log10(x1));
//    return NAN;
//}

void NonLinearFunctionSystem::calculate(const DoubleVector &x0, DoubleVector &x, double epsilon)
{
    unsigned int n = x0.size();
//IPrinter::print(x0,x0.size());
    x = x0;
    DoubleVector x1 = x0;
    DoubleVector dx(n);
    while (true)
    {
        for (unsigned int i=0; i<n; i++)
        {
            x1[i] = fx(x, i);
            dx[i] = x1[i] - x[i];
        }
        x = x1;
        //IPrinter::print(x,x.size());
        //IPrinter::print(dx,dx.size());
        //printf("%12.10f\n", dx.LInfNorm());
        if (dx.L2Norm() <= epsilon) break;
    }
}

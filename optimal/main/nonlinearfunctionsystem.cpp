#include "nonlinearfunctionsystem.h"
#include <math.h>

void NonLinearFunction::Main(int agrc UNUSED_PARAM, char *argv[] UNUSED_PARAM)
{
    NonLinearFunction nlfs;

//    DoubleVector x0;
//    x0 << 3.5 << 2.2;
//    DoubleVector xs;
//    nlfs.calculateSimpleIdetartion(x0, xs, 0.001);
//    IPrinter::print(xs, xs.size());

    DoubleVector x1;
    x1 << 3.4 << 2.2;
    IPrinter::print(x1, x1.size(), 10, 4);
    DoubleVector rx;
    nlfs.calculateNewtonMethod(x1, rx, 0.0001, 0.0001);
    IPrinter::print(rx, rx.size(), 10, 4);

}

double NonLinearFunction::fx(const DoubleVector &x, unsigned int num) const
{
    double x1 = x[0];
    double x2 = x[1];
//    if (num == 0) return sqrt(x1 + 3.0*log10(x1));
//    if (num == 1) return sqrt((x1*(x2+5.0)-1)/2.0);

    if (num == 0) return x1+3.0*log10(x1)-x2*x2;
    if (num == 1) return 2.0*x1*x1 - x1*x2 - 5.0*x1 + 1.0;

    return NAN;
}

void INonLinearFunction::calculateSimpleIdetartion(const DoubleVector &x0, DoubleVector &x, double epsilon)
{
    unsigned int n = x0.size();
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
        if (dx.L2Norm() <= epsilon) break;
    }
}

void INonLinearFunction::calculateNewtonMethod(const DoubleVector &x0, DoubleVector &rx, double diffEspilon, double epsilon)
{
    unsigned int n = x0.size();
    rx = x0;
    DoubleMatrix W(n, n);
    DoubleVector e(x0.size());

    do {
        for (unsigned int row=0; row<n; row++)
        {
            for (unsigned int col=0; col<n; col++)
            {
                DoubleVector x2 = rx;
                DoubleVector x1 = rx;
                x2[col] = x2[col] + diffEspilon;
                x1[col] = x1[col] - diffEspilon;
                W.at(row, col) = (fx(x2, row) - fx(x1, row))/(2.0*diffEspilon);
            }
        }
        W.inverse();
        for (unsigned int row=0; row<n; row++)
        {
            e[row] = 0.0;
            for (unsigned int col=0; col<n; col++)
            {
                e[row] += W[row][col]*fx(rx,col);
            }
            rx[row] -= e[row];
        }
    } while (e.LInfNorm() > epsilon);
}

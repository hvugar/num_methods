#include "cauchyproblemmex.h"

void CauchyProblemMEx::Main(int agrc UNUSED_PARAM, char *argv[] UNUSED_PARAM)
{
    Dimension grid(0.01, 100);
    CauchyProblemMEx exp1(grid);
    DoubleVector y0(2);
    y0[0] = 4.0;
    y0[1] = 0.0;

    DoubleVector yN(2);
    yN[0] = 5.0;
    yN[1] = 2.0;

    {
        DoubleMatrix ry;
        exp1.calculate(0.0, y0, ry, RK2);
        IPrinter::printVector(14, 10, ry.row(0));
        IPrinter::printVector(14, 10, ry.row(1));
    }
    {
        DoubleMatrix ry;
        exp1.calculate(1.0, yN, ry, RK2, R2L);
        IPrinter::printVector(14, 10, ry.row(0));
        IPrinter::printVector(14, 10, ry.row(1));
    }
    IPrinter::printSeperatorLine();
    {
        DoubleMatrix ry;
        exp1.calculate(0.0, y0, ry, RK4);
        IPrinter::printVector(14, 10, ry.row(0));
        IPrinter::printVector(14, 10, ry.row(1));
    }
    {
        DoubleMatrix ry;
        exp1.calculate(1.0, yN, ry, RK4, R2L);
        IPrinter::printVector(14, 10, ry.row(0));
        IPrinter::printVector(14, 10, ry.row(1));
    }
    IPrinter::printSeperatorLine();
    {
        DoubleMatrix ry;
        exp1.calculate(0.0, y0, ry, EULER);
        IPrinter::printVector(14, 10, ry.row(0));
        IPrinter::printVector(14, 10, ry.row(1));
    }
    {
        DoubleMatrix ry;
        exp1.calculate(1.0, yN, ry, EULER, R2L);
        IPrinter::printVector(14, 10, ry.row(0));
        IPrinter::printVector(14, 10, ry.row(1));
    }
}

CauchyProblemMEx::CauchyProblemMEx(const Dimension &grid) : CauchyProblemM(grid)
{}

double CauchyProblemMEx::f(double x, const DoubleVector &y, unsigned int i) const
{
    double y1 = y[0];
    double y2 = y[1];
    if (i == 0)
    {
        return x*y1 - y2 - x;
    }
    if (i == 1)
    {
        return (2.0*x+3.0)*y1 - 2.0*y2 - 6.0*x - 11.0;
    }
    return 0.0;
}

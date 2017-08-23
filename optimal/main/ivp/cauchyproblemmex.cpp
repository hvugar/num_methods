#include "cauchyproblemmex.h"
#include <stdlib.h>
#include <float.h>
#include <math.h>

void CauchyProblemMEx::Main(int agrc UNUSED_PARAM, char *argv[] UNUSED_PARAM)
{
    CauchyProblemMEx exp1;
    exp1.setGrid(ODEGrid(Dimension(0.01, 100, 0)));
    DoubleVector y0(3);
#ifdef SAMPLE_1
    y0[0] = 4.0;//250;
    y0[1] = 0.0;//625;
#endif
#ifdef SAMPLE_2
    y0[0] = 0.0;
    y0[1] = 0.0;
#endif
#ifdef SAMPLE_3
    y0[0] = 7.0;
    y0[1] = 4.0;
    y0[2] = 2.0;
#endif
    double x0 = 1.0;

    //    DoubleVector yN(2);
    //    yN[0] = 5.0;
    //    yN[1] = 2.0;

    {
        //        DoubleMatrix ry;
        //        exp1.calculateCP(x0, y0, ry, RK2);
        //        IPrinter::printVector(14, 10, ry.row(0));
        //        IPrinter::printVector(14, 10, ry.row(1));
    }
    //    {
    //        DoubleMatrix ry;
    //        exp1.calculate(1.0, yN, ry, RK2, R2L);
    //        IPrinter::printVector(14, 10, ry.row(0));
    //        IPrinter::printVector(14, 10, ry.row(1));
    //    }
    IPrinter::printSeperatorLine();
    {
        std::vector<DoubleVector> ry;
        exp1.cauchyProblem(x0, y0, ry, RK4, R2L);
        IPrinter::printVector(14, 10, ry[0]);
        IPrinter::printVector(14, 10, ry[1]);
        IPrinter::printVector(14, 10, ry[2]);

        DoubleVector py;
        exp1.cauchyProblem(x0, y0, py, RK4, R2L);
        //IPrinter::print(py,py.size());
    }
    //    {
    //        DoubleMatrix ry;
    //        exp1.calculate(1.5, yN, ry, RK4, R2L);
    //        IPrinter::printVector(14, 10, ry.row(0));
    //        IPrinter::printVector(14, 10, ry.row(1));
    //    }
    //    IPrinter::printSeperatorLine();
    //    {
    //        DoubleMatrix ry;
    //        exp1.calculate(0.0, y0, ry, EULER);
    //        IPrinter::printVector(14, 10, ry.row(0));
    //        IPrinter::printVector(14, 10, ry.row(1));
    //    }
    //    {
    //        DoubleMatrix ry;
    //        exp1.calculate(1.0, yN, ry, EULER, R2L);
    //        IPrinter::printVector(14, 10, ry.row(0));
    //        IPrinter::printVector(14, 10, ry.row(1));
    //    }
}

double CauchyProblemMEx::f(double x, const DoubleVector &y, unsigned int k UNUSED_PARAM, unsigned int i) const
{
    double y1 = y[0];
    double y2 = y[1];
    double y3 = y[2];
#ifdef SAMPLE_1
    if (i == 0) return x*y1 - y2 - x;
    if (i == 1) return (2.0*x+3.0)*y1 - 2.0*y2 - 6.0*x - 11.0;
#endif
#ifdef SAMPLE_2
    if (i == 0) return 2.0*y1 + 3.0*y2 + 2.0 - 13.0*x;
    if (i == 1) return 5.0*y1 + 4.0*y2 + 3.0 - 22.0*x;
#endif
#ifdef SAMPLE_3
    if (i==0) return (+2.0*y1 - 3.0*y2 + 1.0*y3) + (+3.0     - (2.0*(3.0*x+4.0) - 3.0*(4.0*x*x) + 1.0*(x*x+x)));
    if (i==1) return (+3.0*y1 + 1.0*y2 - 2.0*y3) + (+8.0*x   - (3.0*(3.0*x+4.0) + 1.0*(4.0*x*x) - 2.0*(x*x+x)));
    if (i==2) return (+1.0*y1 - 5.0*y2 - 3.0*y3) + (+2.0*x+1 - (1.0*(3.0*x+4.0) - 5.0*(4.0*x*x) - 3.0*(x*x+x)));
#endif
    return NAN;
}

double CauchyProblemMEx::y(double x, unsigned int k UNUSED_PARAM, unsigned int i) const
{
#ifdef SAMPLE_1
    if (i==0) return 2.0*x;
    if (i==1) return 3.0*x;
#endif
#ifdef SAMPLE_2
    if (i==0) return 2.0*x;
    if (i==1) return 3.0*x;
#endif
#ifdef SAMPLE_3
    if (i==0) return 3.0*x+4.0;
    if (i==1) return 4.0*x*x;
    if (i==2) return x*x+x;
#endif
    return NAN;
}

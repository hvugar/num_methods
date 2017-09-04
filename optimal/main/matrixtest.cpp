#include "matrixtest.h"
#include <iomanip>

void MatrixTest::Main(int agrc UNUSED_PARAM, char *argv[] UNUSED_PARAM)
{   
    DoubleMatrix m1(3,4);

    m1[0][0] = 5.0; m1[0][1] = 2.0; m1[0][2] = 3.0; m1[0][3] = 1.0;
    m1[1][0] = 4.0; m1[1][1] = 1.0; m1[1][2] = 2.0; m1[1][3] = 0.0;
    m1[2][0] = 1.0; m1[2][1] = 5.0; m1[2][2] = 7.0; m1[2][3] = 3.0;

    DoubleVector v1(4);

    v1[0] = 2.0; v1[1] = 3.0; v1[2] = 0.0; v1[3] = 1.0;

    IPrinter::print(m1);
    IPrinter::printSeperatorLine();

    IPrinter::print(v1);
    IPrinter::printSeperatorLine();

    DoubleMatrix m2 = m1*v1;
    IPrinter::print(m2);
    IPrinter::printSeperatorLine();

    m1 *= v1;
    IPrinter::print(m1);
    IPrinter::printSeperatorLine();

    DoubleMatrix m3(5, 6, 1.0);
    IPrinter::print(m3);
    IPrinter::printSeperatorLine();
    m3.resize(3, 3, 2.0);
    IPrinter::print(m3,m3.rows(),m3.cols());

    IPrinter::printSeperatorLine();
}

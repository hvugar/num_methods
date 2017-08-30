#include "matrixtest.h"

void MatrixTest::Main(int agrc UNUSED_PARAM, char *argv[] UNUSED_PARAM)
{
    DoubleMatrix m1(4,3);
    m1.randomData();

    DoubleMatrix m2(4,3);
    m2.randomData();

    DoubleMatrix m3(4,3);
    m3.randomData();

    IPrinter::printSeperatorLine();
    IPrinter::print(m1, 0, 0);
    IPrinter::printSeperatorLine();

    IPrinter::printSeperatorLine();
    IPrinter::print(m2, 0, 0);
    IPrinter::printSeperatorLine();

    m2 += m1;

    IPrinter::printSeperatorLine();
    IPrinter::print(m1, 0, 0);
    IPrinter::printSeperatorLine();

    IPrinter::printSeperatorLine();
    IPrinter::print(m2, 0, 0);
    IPrinter::printSeperatorLine();
}

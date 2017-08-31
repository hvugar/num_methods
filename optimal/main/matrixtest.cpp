#include "matrixtest.h"

void MatrixTest::Main(int agrc UNUSED_PARAM, char *argv[] UNUSED_PARAM)
{
    DoubleMatrix m1(4,3);

    DoubleMatrix m2(4,3);

    DoubleMatrix m3(4,3);

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

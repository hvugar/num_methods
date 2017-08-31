#include "matrixtest.h"
#include <iomanip>

void MatrixTest::Main(int agrc UNUSED_PARAM, char *argv[] UNUSED_PARAM)
{
    DoubleVector v1;
    v1 << 1.0 << 3.0 << 5.0;
    IPrinter::print(v1);

    DoubleVector v2;
    v2 << 4.0 << 8.0 << 2.0;
    IPrinter::print(v2);

    DoubleVector v3 = v1*v2;
    IPrinter::print(v1);
    IPrinter::print(v3);
}

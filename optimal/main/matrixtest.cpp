#include "matrixtest.h"
#include <iomanip>
#include <utils/random.h>

void MatrixTest::Main(int agrc UNUSED_PARAM, char *argv[] UNUSED_PARAM)
{
    DoubleVector vec(5, 0.0);
    Random::fillVector(vec, 0, 5, 5);
    DoubleMatrix m = DoubleMatrix::IdentityMatrix(5);
    //DoubleMatrix m = DoubleMatrix::DiagonalMatrix(vec);
    IPrinter::print(m);
    m.inverse();
    IPrinter::printSeperatorLine();
    IPrinter::print(m);
}

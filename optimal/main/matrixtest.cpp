#include "matrixtest.h"

void MatrixTest::Main(int agrc UNUSED_PARAM, char *argv[] UNUSED_PARAM)
{
    DoubleMatrix m1(4,3);
    m1.randomData();

    DoubleMatrix m2(4,5);
    m2.randomData();

    //m2 = m1;

    IPrinter::printSeperatorLine();
    IPrinter::print(m1, 0, 0);
    IPrinter::printSeperatorLine();

    IPrinter::printSeperatorLine();
    IPrinter::print(m2, 0, 0);
    IPrinter::printSeperatorLine();

    if (m1 == m2)
        puts("not equals");
    else
        puts("equals");

    DoubleMatrix im = DoubleMatrix::IdentityMatrix(5);
    IPrinter::printSeperatorLine();
    IPrinter::print(im);
    IPrinter::printSeperatorLine();
    im.inverse();
    IPrinter::print(im);
    IPrinter::printSeperatorLine();

    DoubleMatrix hm = DoubleMatrix::HilbertMatrix(10, 10);
    IPrinter::print(hm);
    IPrinter::printSeperatorLine();
    printf("%.10f\n", hm.determinant());
}

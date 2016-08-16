#include "headers.h"
#include <time.h>
#include "gaussianelimination.h"
#include "problem1.h"
#include "problem2.h"

#include <QtGui/QGuiApplication>
#include <imaging.h>

int main(int argc, char *argv[])
{
    C_UNUSED(argc);
    C_UNUSED(argv);
    srand(time(NULL));

//    Problem1 p;
//    p.calculate3();

    DoubleMatrix m;

    Problem2 p2;
    p2.calculate1(m, p2.ht, p2.hx, p2.M, p2.N, p2.lambdaM, p2.lambdaL, p2.lambdaR, p2.a);

    IPrinter::printMatrix(m);

    return 0;
}

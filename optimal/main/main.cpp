#include "headers.h"
#include <time.h>
#include "gaussianelimination.h"
#include "problem1.h"

#include <QtGui/QGuiApplication>

int main(int argc, char *argv[])
{
    C_UNUSED(argc);
    C_UNUSED(argv);
    srand(time(NULL));

//    QGuiApplication app(argc, argv);

    Problem1 p;
    p.calculate4();

    //GaussianEliminationTester::main(argc, argv);
    //HeatControl2D::main(argc, argv);
    return 0;
}

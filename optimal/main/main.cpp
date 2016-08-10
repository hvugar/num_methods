#include "headers.h"
#include <time.h>
#include "gaussianelimination.h"
#include "problem1.h"

int main(int argc, char *argv[])
{
    C_UNUSED(argc);
    C_UNUSED(argv);
    srand(time(NULL));

    Problem1 p;
    p.calculate();

    //GaussianEliminationTester::main(argc, argv);
    //HeatControl2D::main(argc, argv);
    return 0;
}

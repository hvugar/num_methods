#include "problem1.h"
#include "problem1k.h"
#include "problem1z.h"
#include "problem1x.h"
#include "problem1x1.h"
#include "problem1kz.h"
#include "problem1kzx.h"

#include "problem3.h"
#include "loadedsystems.h"
#include "example1.h"
#include <cmethods.h>

#include <QtGui/QGuiApplication>
#include <imaging.h>

#include <float.h>

int main(int argc, char *argv[])
{
    C_UNUSED(argc);
    C_UNUSED(argv);

    Problem1KZX::Main(argc, argv);

    //srand(time(NULL));
    //Problem1K::Main(argc, argv);
    //    LoadedSystems ls;

    //    Example1 ex;
    //    ex.calculate();

//    unsigned int n = 7;
//    double a[] = {0.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0};
//    double b[] = {1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0};
//    double c[] = {2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 0.0};
//    double d[] = {5.0, 8.0, 8.0, 10.0, 7.0, 8.0, 5.0};
//    double *x = (double*)malloc(sizeof(double)*n);
//    double e[] = {1.0, 2.0, 0.0, 3.0, 0.0, 2.0, 0.0};
////    unsigned int E[] = {3, 4};
////    qovmaE(a, b, c, d, x, n, e, E, 2);
//    qovma2(a, b, c, d, x, n, e);
//    printf("%10.6f %10.6f %10.6f %10.6f %10.6f %10.6f %10.6f\n", x[0], x[1], x[2], x[3], x[4], x[5], x[6]);

    return 0;
}

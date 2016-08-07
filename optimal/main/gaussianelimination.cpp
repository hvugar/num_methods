#include "gaussianelimination.h"
#include <matrix2d.h>
#include <vector2d.h>
#include <printer.h>
#include <cmethods.h>

void GaussianEliminationTester::main(int argc, char **argv)
{
    unsigned int n = 5;

//    double **a = (double**) malloc(sizeof(double*))

    DoubleMatrix m(n,n);
    m.randomData();
    DoubleVector x(n);
    x.randomData();
    DoubleVector b(n);

    puts("x");
    for (unsigned int i=0; i<x.size(); i++)
    {
        printf("%10.6f ", x[i]);
    }
    puts("\nm");
    //m.at(0,0) = 0.0;
    m.print();
    for (unsigned int j=0; j<n; j++)
    {
        b[j] = 0;
        for (unsigned int i=0; i<n; i++)
        {
            b[j] += m[j][i]*x[i];
        }
    }
    puts("b");
    for (unsigned int i=0; i<b.size(); i++)
    {
        printf("%10.6f ", b[i]);
    }
    puts("");
    puts("--------------");

    gaussianElimination(m.data(), b.data(), x.data(), n);

}

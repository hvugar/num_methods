#include "gaussianelimination.h"
#include <matrix2d.h>
#include <vector2d.h>
#include <printer.h>
#include <cmethods.h>

void GaussianEliminationTester::main(int argc __attribute__ ((unused)), char **argv __attribute__ ((unused)))
{
    unsigned int n = 10;
    DoubleMatrix m(n,n);
    m.randomData();
    //for (unsigned int i=0; i<n-1; i++) m.at(0,i) = 0.0;
    //for (unsigned int i=1; i<n; i++) m.at(1,i) = 0.0;
    m.print();
    puts("");

    DoubleVector x(n, 0.0);
    x.randomData();
    for (unsigned int i=0; i<x.size(); i++) x[i] = x[i]*0.01;
    x.print();

    DoubleVector b(n);
    for (unsigned int j=0; j<n; j++)
    {
        b[j] = 0;
        for (unsigned int i=0; i<n; i++)
        {
            b[j] += m.at(j,i)*x[i];
        }
    }
    b.print();

    for (unsigned int i=0; i<x.size(); i++) x[i] = 0.0;
    //gaussianElimination(m.data(), b.data(), x.data(), n);
    GaussianElimination(m, b, x);
    x.print();
}

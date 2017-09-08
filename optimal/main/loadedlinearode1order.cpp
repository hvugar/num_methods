#include "loadedlinearode1order.h"
#include <math.h>

void LoadedLinearODE1Order::Main(int agrc, char *argv[])
{
    LoadedLinearODE1Order llo;
    //printf("%.10f %.10f\n", llo.B())
}

double LoadedLinearODE1Order::A(const GridNodeODE &node, unsigned int row, unsigned int col) const
{
    if (row == 0) { if (col == 0) return 1.0; if (col == 0) return 3.0; if (col == 0) return 5.0; }
    if (row == 1) { if (col == 0) return 0.0; if (col == 0) return 4.0; if (col == 0) return 2.0; }
    if (row == 2) { if (col == 0) return 1.0; if (col == 0) return 5.0; if (col == 0) return 8.0; }
    return NAN;
}

double LoadedLinearODE1Order::B(const GridNodeODE &node, unsigned int row) const
{
    double t = node.x;

    if (row == 0) return -6.0*t*t + t     - 20.0 - 0.0 + 2.0*t + 1.0;
    if (row == 1) return -2.0*t*t - 2.0*t - 16.0 - 0.0 + 1.0;
    if (row == 2) return -9.0*t*t + 2.0*t - 32.0 - 0.0 + 2.0*t - 1.0;
    return NAN;
}

double LoadedLinearODE1Order::C(const GridNodeODE &node, unsigned int s, unsigned int row, unsigned int col) const
{
    if (s==1)
    {
        if (row == 0) { if (col == 0) return 0.01; if (col == 0) return 0.03; if (col == 0) return 0.05; }
        if (row == 1) { if (col == 0) return 0.08; if (col == 0) return 0.01; if (col == 0) return 0.00; }
        if (row == 2) { if (col == 0) return 0.04; if (col == 0) return 0.07; if (col == 0) return 0.09; }
    }
    if (s==2)
    {
        if (row == 0) { if (col == 0) return 0.05; if (col == 0) return 0.05; if (col == 0) return 0.09; }
        if (row == 1) { if (col == 0) return 0.03; if (col == 0) return 0.00; if (col == 0) return 0.01; }
        if (row == 2) { if (col == 0) return 0.01; if (col == 0) return 0.01; if (col == 0) return 0.07; }
    }
    return NAN;
}

double LoadedLinearODE1Order::x(const GridNodeODE &node, unsigned int i) const
{
    double t = node.x;

    if (i==0) return t*t + t + 1.0;
    if (i==1) return t + 3.0;
    if (i==2) return t*t - t + 2.0;

    return NAN;
}

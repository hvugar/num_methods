#include "loadedlinearode1order.h"
#include <math.h>

void LoadedLinearODE1Order::Main(int argc, char *argv[])
{
    LoadedLinearODE1Order llo;

//    unsigned int row = 2;
//    GridNodeODE node1(0.3,0);
//    GridNodeODE node2(0.6,0);
//    printf("%.10f\n", llo.C(node1,1,row,0)*llo.x(node1,0) + llo.C(node1,1,row,1)*llo.x(node1,1) + llo.C(node1,1,row,2)*llo.x(node1,2)
//                     +llo.C(node2,2,row,0)*llo.x(node2,0) + llo.C(node2,2,row,1)*llo.x(node2,1) + llo.C(node2,2,row,2)*llo.x(node2,2));
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

    if (row == 0) return -6.0*t*t + 3.0*t - 19.0 - 0.6388000000;
    if (row == 1) return -2.0*t*t - 2.0*t - 15.0 - 0.2206000000;
    if (row == 2) return -9.0*t*t + 4.0*t - 33.0 - 0.6265000000;

    return NAN;
}

double LoadedLinearODE1Order::C(const GridNodeODE &node UNUSED_PARAM, unsigned int s, unsigned int row, unsigned int col) const
{
    if (s==1)
    {
        if (row == 0) { if (col == 0) return 0.01; if (col == 1) return 0.03; if (col == 2) return 0.05; }
        if (row == 1) { if (col == 0) return 0.08; if (col == 1) return 0.01; if (col == 2) return 0.00; }
        if (row == 2) { if (col == 0) return 0.04; if (col == 1) return 0.07; if (col == 2) return 0.09; }
    }
    if (s==2)
    {
        if (row == 0) { if (col == 0) return 0.05; if (col == 1) return 0.05; if (col == 2) return 0.09; }
        if (row == 1) { if (col == 0) return 0.03; if (col == 1) return 0.00; if (col == 2) return 0.01; }
        if (row == 2) { if (col == 0) return 0.01; if (col == 1) return 0.01; if (col == 2) return 0.07; }
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

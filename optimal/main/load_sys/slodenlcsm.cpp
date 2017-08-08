#include "slodenlcsm.h"

void SystemLinearODENonLocalContionsM::Main(int agrc UNUSED_PARAM, char *argv[] UNUSED_PARAM)
{
//    ODEGrid grid(Dimension(0.01, 100, 0));
//    SystemLinearODENonLocalContionsM cpnlcs(grid);

//    cpnlcs.initialize();
//    DoubleVector x;
//    cpnlcs.calculateForward(x);
//    IPrinter::print(x,x.size());
//    DoubleMatrix m;
//    cpnlcs.calculateBackwardCP(x, m);
//    for (unsigned int row=0; row<cpnlcs.systemOrder(); row++) IPrinter::printVector(m.row(row));
}

SystemLinearODENonLocalContionsM::SystemLinearODENonLocalContionsM(const ODEGrid &grid)
    : ISystemLinearODENonLocalContionsM(grid) {}

double SystemLinearODENonLocalContionsM::A(double t, unsigned int k, unsigned int row, unsigned int col) const
{
    if (row == 1)
    {
        if (col == 1) { return 1.0; }
        if (col == 2) { return 3.0; }
        if (col == 3) { return 9.0; }
    }
    if (row == 2)
    {
        if (col == 1) { return 8.0; }
        if (col == 2) { return 4.0; }
        if (col == 3) { return 1.0; }
    }
    if (row == 3)
    {
        if (col == 1) { return 0.0; }
        if (col == 2) { return 2.0; }
        if (col == 3) { return 3.0; }
    }
    return NAN;
}

double SystemLinearODENonLocalContionsM::B(double t, unsigned int k, unsigned int row, unsigned int col) const
{
    if (row == 1)
    {
        if (col == 1) { return dX(t,k,1,1) - (A(t,k,1,1)*X(t,k,1,1) + A(t,k,1,2)*X(t,k,2,1) + A(t,k,1,3)*X(t,k,3,1)); }
        if (col == 2) { return dX(t,k,1,2) - (A(t,k,1,1)*X(t,k,1,2) + A(t,k,1,2)*X(t,k,2,2) + A(t,k,1,3)*X(t,k,3,2)); }
        if (col == 3) { return dX(t,k,1,3) - (A(t,k,1,1)*X(t,k,1,3) + A(t,k,1,2)*X(t,k,2,3) + A(t,k,1,3)*X(t,k,3,3)); }
    }
    if (row == 2)
    {
        if (col == 1) { return dX(t,k,2,1) - (A(t,k,2,1)*X(t,k,1,1) + A(t,k,2,2)*X(t,k,2,1) + A(t,k,2,3)*X(t,k,3,1)); }
        if (col == 2) { return dX(t,k,2,2) - (A(t,k,2,1)*X(t,k,1,2) + A(t,k,2,2)*X(t,k,2,2) + A(t,k,2,3)*X(t,k,3,2)); }
        if (col == 3) { return dX(t,k,2,3) - (A(t,k,2,1)*X(t,k,1,3) + A(t,k,2,2)*X(t,k,2,3) + A(t,k,2,3)*X(t,k,3,3)); }
    }
    if (row == 3)
    {
        if (col == 1) { return dX(t,k,3,1) - (A(t,k,3,1)*X(t,k,1,1) + A(t,k,3,2)*X(t,k,2,1) + A(t,k,3,3)*X(t,k,3,1)); }
        if (col == 2) { return dX(t,k,3,2) - (A(t,k,3,1)*X(t,k,1,2) + A(t,k,3,2)*X(t,k,2,2) + A(t,k,3,3)*X(t,k,3,2)); }
        if (col == 3) { return dX(t,k,3,3) - (A(t,k,3,1)*X(t,k,1,3) + A(t,k,3,2)*X(t,k,2,3) + A(t,k,3,3)*X(t,k,3,3)); }
    }
    return NAN;
}

double SystemLinearODENonLocalContionsM::X(double t, unsigned int k, unsigned int row, unsigned int col) const
{
    if (row == 1)
    {
        if (col == 1) { return t; }
        if (col == 2) { return t*t; }
        if (col == 3) { return 2.0*t; }
    }
    if (row == 2)
    {
        if (col == 1) { return 3.0*t; }
        if (col == 2) { return t*t*t; }
        if (col == 3) { return 5.0*t; }
    }
    if (row == 3)
    {
        if (col == 1) { return -3.0*t; }
        if (col == 2) { return t*t+t; }
        if (col == 3) { return 8.0*t; }
    }
    return NAN;
}

double SystemLinearODENonLocalContionsM::dX(double t, unsigned int k, unsigned int row, unsigned int col) const
{
    if (row == 1)
    {
        if (col == 1) { return 1.0; }
        if (col == 2) { return 2.0*t; }
        if (col == 3) { return 2.0; }
    }
    if (row == 2)
    {
        if (col == 1) { return 3.0; }
        if (col == 2) { return 3.0*t*t; }
        if (col == 3) { return 5.0; }
    }
    if (row == 3)
    {
        if (col == 1) { return -3.0; }
        if (col == 2) { return 2.0*t+1; }
        if (col == 3) { return 8.0; }
    }
    return NAN;
}

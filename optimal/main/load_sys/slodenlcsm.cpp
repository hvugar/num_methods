#include "slodenlcsm.h"
#include <load_sys/islodenlcsv.h>

void SystemLinearODENonLocalContionsM::Main(int agrc UNUSED_PARAM, char *argv[] UNUSED_PARAM)
{
    SystemLinearODENonLocalContionsM cpnlcs;
    cpnlcs.setGrid(ODEGrid(Dimension(0.01, 100, 0)));

    cpnlcs.initialize();
    DoubleMatrix x;
    cpnlcs.calculateForward(x);
    IPrinter::print(x, x.rows(), x.cols());

    std::vector<std::vector<DoubleVector>> m;
    cpnlcs.calculateBackwardCP(x, m);

    IPrinter::printSeperatorLine();

    for (unsigned int r=0; r<m.size(); r++)
    {
        for (unsigned int c=0; c<m[r].size(); c++)
        {
            IPrinter::printVector(14, 10, m[r][c]);
        }
    }
}

void SystemLinearODENonLocalContionsM::initialize()
{
    Dimension dim = grid().dimension();
    unsigned int N = dim.sizeN();

    unsigned int n0 = 3;
    unsigned int n1 = 0;
    unsigned int n2 = 0;
    unsigned int n = n0 + n1 + n2;

    ISystemLinearODENonLocalContionsV::Condition nsc0;
    nsc0.type = ISystemLinearODENonLocalContionsV::NonSeparated;
    nsc0.time = 0.0;
    nsc0.nmbr = 0;
    nsc0.alpha.resize(n0, n);
    for (unsigned int row=0; row<n0; row++) for(unsigned int col=0; col<n; col++) nsc0.alpha[row][col] = (rand() % 1000) / 1000.0;
    addNonSeparatedCondition(nsc0);

    ISystemLinearODENonLocalContionsV::Condition nsc1;
    nsc1.type = ISystemLinearODENonLocalContionsV::NonSeparated;
    nsc1.time = 0.25;
    nsc1.nmbr = N/4;
    nsc1.alpha.resize(n0, n);
    for (unsigned int row=0; row<n0; row++) for(unsigned int col=0; col<n; col++) nsc1.alpha[row][col] = (rand() % 1000) / 1000.0;
    addNonSeparatedCondition(nsc1);

    ISystemLinearODENonLocalContionsV::Condition nsc2;
    nsc2.type = ISystemLinearODENonLocalContionsV::NonSeparated;
    nsc2.time = 0.5;
    nsc2.nmbr = N/2;
    nsc2.alpha.resize(n0, n);
    for (unsigned int row=0; row<n0; row++) for(unsigned int col=0; col<n; col++) nsc2.alpha[row][col] = (rand() % 1000) / 1000.0;
    addNonSeparatedCondition(nsc2);

    ISystemLinearODENonLocalContionsV::Condition nsc3;
    nsc3.type = ISystemLinearODENonLocalContionsV::NonSeparated;
    nsc3.time = 0.75;
    nsc3.nmbr = 3*(N/4);
    nsc3.alpha.resize(n0, n);
    for (unsigned int row=0; row<n0; row++) for(unsigned int col=0; col<n; col++) nsc3.alpha[row][col] = (rand() % 1000) / 1000.0;
    addNonSeparatedCondition(nsc3);

    ISystemLinearODENonLocalContionsV::Condition nsc4;
    nsc4.type = ISystemLinearODENonLocalContionsV::NonSeparated;
    nsc4.time = 1.0;
    nsc4.nmbr = N;
    nsc4.alpha.resize(n0, n);
    for (unsigned int row=0; row<n0; row++) for(unsigned int col=0; col<n; col++) nsc4.alpha[row][col] = (rand() % 1000) / 1000.0;
    addNonSeparatedCondition(nsc4);

    DoubleMatrix betta(n0, n);
    unsigned int L = nonSeparatedConditions().size();

    for (unsigned int row=0; row<n0; row++)
    {
        for (unsigned int col=0; col<n; col++)
        {
            betta[row][col] = 0.0;
            for (unsigned int s=0; s<L; s++)
            {
                const ISystemLinearODENonLocalContionsV::Condition &c = nonSeparatedConditions().at(s);
                for (unsigned int i=0; i<n; i++) betta[row][col] += c.alpha[row][i] * X(c.time, 0, i, col);
            }
        }
    }

    setBetta(betta);

    //    IPrinter::print(nsc0.alpha, nsc0.alpha.rows(), nsc0.alpha.cols());
    //    IPrinter::printSeperatorLine();
    //    IPrinter::print(nsc4.alpha, nsc4.alpha.rows(), nsc4.alpha.cols());
    //    IPrinter::printSeperatorLine();
    //    IPrinter::print(betta, betta.rows(), betta.cols());
}

double SystemLinearODENonLocalContionsM::A(double, unsigned int, unsigned int row, unsigned int col) const
{
    if (row == 0)
    {
        if (col == 0) { return 1.0; }
        if (col == 1) { return 3.0; }
        if (col == 2) { return 9.0; }
    }
    if (row == 1)
    {
        if (col == 0) { return 8.0; }
        if (col == 1) { return 4.0; }
        if (col == 2) { return 1.0; }
    }
    if (row == 2)
    {
        if (col == 0) { return 0.0; }
        if (col == 1) { return 2.0; }
        if (col == 2) { return 3.0; }
    }
    return NAN;
}

double SystemLinearODENonLocalContionsM::B(double t, unsigned int k, unsigned int row, unsigned int col) const
{
    if (row == 0)
    {
        if (col == 0) { return dX(t,k,0,0) - (A(t,k,0,0)*X(t,k,0,0) + A(t,k,0,1)*X(t,k,1,0) + A(t,k,0,2)*X(t,k,2,0)); }
        if (col == 1) { return dX(t,k,0,1) - (A(t,k,0,0)*X(t,k,0,1) + A(t,k,0,1)*X(t,k,1,1) + A(t,k,0,2)*X(t,k,2,1)); }
        if (col == 2) { return dX(t,k,0,2) - (A(t,k,0,0)*X(t,k,0,2) + A(t,k,0,1)*X(t,k,1,2) + A(t,k,0,2)*X(t,k,2,2)); }
    }
    if (row == 1)
    {
        if (col == 0) { return dX(t,k,1,0) - (A(t,k,1,0)*X(t,k,0,0) + A(t,k,1,1)*X(t,k,1,0) + A(t,k,1,2)*X(t,k,2,0)); }
        if (col == 1) { return dX(t,k,1,1) - (A(t,k,1,0)*X(t,k,0,1) + A(t,k,1,1)*X(t,k,1,1) + A(t,k,1,2)*X(t,k,2,1)); }
        if (col == 2) { return dX(t,k,1,2) - (A(t,k,1,0)*X(t,k,0,2) + A(t,k,1,1)*X(t,k,1,2) + A(t,k,1,2)*X(t,k,2,2)); }
    }
    if (row == 2)
    {
        if (col == 0) { return dX(t,k,2,0) - (A(t,k,2,0)*X(t,k,0,0) + A(t,k,2,1)*X(t,k,1,0) + A(t,k,2,2)*X(t,k,2,0)); }
        if (col == 1) { return dX(t,k,2,1) - (A(t,k,2,0)*X(t,k,0,1) + A(t,k,2,1)*X(t,k,1,1) + A(t,k,2,2)*X(t,k,2,1)); }
        if (col == 2) { return dX(t,k,2,2) - (A(t,k,2,0)*X(t,k,0,2) + A(t,k,2,1)*X(t,k,1,2) + A(t,k,2,2)*X(t,k,2,2)); }
    }
    return NAN;
}

double SystemLinearODENonLocalContionsM::X(double t, unsigned int, unsigned int row, unsigned int col) const
{
    if (row == 0)
    {
        if (col == 0) { return t; }
        if (col == 1) { return t*t; }
        if (col == 2) { return 2.0*t; }
    }
    if (row == 1)
    {
        if (col == 0) { return 3.0*t; }
        if (col == 1) { return t*t*t; }
        if (col == 2) { return 5.0*t; }
    }
    if (row == 2)
    {
        if (col == 0) { return -3.0*t; }
        if (col == 1) { return t*t+t; }
        if (col == 2) { return 8.0*t; }
    }
    return NAN;
}

double SystemLinearODENonLocalContionsM::dX(double t, unsigned int, unsigned int row, unsigned int col) const
{
    if (row == 0)
    {
        if (col == 0) { return 1.0; }
        if (col == 1) { return 2.0*t; }
        if (col == 2) { return 2.0; }
    }
    if (row == 1)
    {
        if (col == 0) { return 3.0; }
        if (col == 1) { return 3.0*t*t; }
        if (col == 2) { return 5.0; }
    }
    if (row == 2)
    {
        if (col == 0) { return -3.0; }
        if (col == 1) { return 2.0*t+1.0; }
        if (col == 2) { return 8.0; }
    }
    return NAN;
}

#include "llode1oex1.h"
#include <math.h>
#include <printer.h>

void LoadLinearODE1stOrderEx1::Main(int agrc UNUSED_PARAM, char *argv[] UNUSED_PARAM)
{
    LoadLinearODE1stOrderEx1 cpnlcs;
    cpnlcs.setGrid(ODEGrid(Dimension(0.01, 100, 0)));

    std::vector<Condition> nscs;
    DoubleVector betta;
    std::vector<DoubleVector> x;
    cpnlcs.initialize(nscs, betta);
}

void LoadLinearODE1stOrderEx1::initialize(std::vector<Condition> &nscs, DoubleVector &betta)
{
    Dimension dim = grid().dimension();
    unsigned int N = dim.sizeN();

    unsigned int n = 3;

    Condition nsc0;
    nsc0.time = 0.0;
    nsc0.nmbr = 0;
    nsc0.mtrx.resize(n, n);
    for (unsigned int row=0; row<n; row++) for(unsigned int col=0; col<n; col++) nsc0.mtrx[row][col] = (rand() % 1000) / 1000.0;

    Condition nsc1;
    nsc1.time = 0.25;
    nsc1.nmbr = N/4;
    nsc1.mtrx.resize(n, n);
    for (unsigned int row=0; row<n; row++) for(unsigned int col=0; col<n; col++) nsc1.mtrx[row][col] = (rand() % 1000) / 1000.0;

    Condition nsc2;
    nsc2.time = 0.5;
    nsc2.nmbr = N/2;
    nsc2.mtrx.resize(n, n);
    for (unsigned int row=0; row<n; row++) for(unsigned int col=0; col<n; col++) nsc2.mtrx[row][col] = (rand() % 1000) / 1000.0;

    Condition nsc3;
    nsc3.time = 0.75;
    nsc3.nmbr = 3*(N/4);
    nsc3.mtrx.resize(n, n);
    for (unsigned int row=0; row<n; row++) for(unsigned int col=0; col<n; col++) nsc3.mtrx[row][col] = (rand() % 1000) / 1000.0;

    Condition nsc4;
    nsc4.time = 1.0;
    nsc4.nmbr = N;
    nsc4.mtrx.resize(n, n);
    for (unsigned int row=0; row<n; row++) for(unsigned int col=0; col<n; col++) nsc4.mtrx[row][col] = (rand() % 1000) / 1000.0;

    nscs.push_back(nsc0);
    nscs.push_back(nsc1);
    nscs.push_back(nsc2);
    nscs.push_back(nsc3);
    nscs.push_back(nsc4);

    unsigned int L = nscs.size();

    betta.resize(n);
    for (unsigned int row=0; row<n; row++)
    {
        betta[row] = 0.0;
        for (unsigned int s=0; s<L; s++)
        {
            const Condition &c = nscs.at(s);
            for (unsigned int i=0; i<n; i++) betta[row] += c.mtrx[row][i] * X(c.time, i);
        }
    }
}

double LoadLinearODE1stOrderEx1::A(double t UNUSED_PARAM, unsigned int, unsigned int row, unsigned int col) const
{
    if (row==0) { if (col==0) { return +2.0; } if (col==1) { return -3.0; } if (col==2) { return +1.0; } }
    if (row==1) { if (col==0) { return +3.0; } if (col==1) { return +1.0; } if (col==2) { return -2.0; } }
    if (row==2) { if (col==0) { return +1.0; } if (col==1) { return -5.0; } if (col==2) { return -3.0; } }
    return NAN;
}

double LoadLinearODE1stOrderEx1::B(double t, unsigned int, unsigned int row) const
{
    if (row==0) return 3.0     - (2.0*(3.0*t+4.0) - 3.0*(4.0*t*t) + 1.0*(t*t+t));
    if (row==1) return 8.0*t   - (3.0*(3.0*t+4.0) + 1.0*(4.0*t*t) - 2.0*(t*t+t));
    if (row==2) return 2.0*t+1 - (1.0*(3.0*t+4.0) - 5.0*(4.0*t*t) - 3.0*(t*t+t));
    return NAN;
}

double LoadLinearODE1stOrderEx1::X(double t, int row) const
{
    if (row==0) return 3.0*t+4.0;
    if (row==1) return 4.0*t*t;
    if (row==2) return t*t+t;
    return NAN;
}

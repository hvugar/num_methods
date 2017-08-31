#include "slodenlcsv.h"
#include <math.h>

void SystemLinearODENonLocalContionsV::Main(int agrc UNUSED_PARAM, char *argv[] UNUSED_PARAM)
{
    SystemLinearODENonLocalContionsV cpnlcs;
    //cpnlcs.setGrid(ODEGrid(Dimension(0.01, 100, 0)));
    cpnlcs.initialize();
    DoubleVector x;
    cpnlcs.calculateForward(x);

    IPrinter::print(x,x.length());
    std::vector<DoubleVector> m;
    cpnlcs.calculateBackwardCP(x, m);
    for (unsigned int row=0; row<cpnlcs.systemOrder(); row++) IPrinter::printVector(m.at(row));
}

void SystemLinearODENonLocalContionsV::initialize()
{
    Dimension dim = grid().dimension();
    unsigned int N = dim.sizeN();

    unsigned int n0 = 3;
    unsigned int n1 = 0;
    unsigned int n2 = 0;
    unsigned int n = n0 + n1 + n2;

    Condition nsc0;
    nsc0.type = NonSeparated;
    nsc0.time = 0.0;
    nsc0.nmbr = 0;
    nsc0.mtrx.resize(n0, n);
    for (unsigned int row=0; row<n0; row++) for(unsigned int col=0; col<n; col++) nsc0.mtrx[row][col] = (rand() % 1000) / 1000.0;
    addNonSeparatedCondition(nsc0);

    Condition nsc1;
    nsc1.type = NonSeparated;
    nsc1.time = 0.25;
    nsc1.nmbr = N/4;
    nsc1.mtrx.resize(n0, n);
    for (unsigned int row=0; row<n0; row++) for(unsigned int col=0; col<n; col++) nsc1.mtrx[row][col] = (rand() % 1000) / 1000.0;
    addNonSeparatedCondition(nsc1);

    Condition nsc2;
    nsc2.type = NonSeparated;
    nsc2.time = 0.5;
    nsc2.nmbr = N/2;
    nsc2.mtrx.resize(n0, n);
    for (unsigned int row=0; row<n0; row++) for(unsigned int col=0; col<n; col++) nsc2.mtrx[row][col] = (rand() % 1000) / 1000.0;
    addNonSeparatedCondition(nsc2);

    Condition nsc3;
    nsc3.type = NonSeparated;
    nsc3.time = 0.75;
    nsc3.nmbr = 3*(N/4);
    nsc3.mtrx.resize(n0, n);
    for (unsigned int row=0; row<n0; row++) for(unsigned int col=0; col<n; col++) nsc3.mtrx[row][col] = (rand() % 1000) / 1000.0;
    addNonSeparatedCondition(nsc3);

    Condition nsc4;
    nsc4.type = NonSeparated;
    nsc4.time = 1.0;
    nsc4.nmbr = N;
    nsc4.mtrx.resize(n0, n);
    for (unsigned int row=0; row<n0; row++) for(unsigned int col=0; col<n; col++) nsc4.mtrx[row][col] = (rand() % 1000) / 1000.0;
    addNonSeparatedCondition(nsc4);

    DoubleVector betta(n);
    unsigned int L = nonSeparatedConditions().size();

    for (unsigned int row=0; row<n0; row++)
    {
        betta[row] = 0.0;
        for (unsigned int s=0; s<L; s++)
        {
            const Condition &c = nonSeparatedConditions().at(s);
            for (unsigned int i=0; i<n; i++) betta[row] += c.mtrx[row][i] * x(c.time, i);
        }
    }

    Condition lscs;
    lscs.type = SeparatedLeft;
    lscs.time = 0.0;
    lscs.nmbr = 0;
    lscs.mtrx.resize(n1, n);
    for (unsigned int row=0; row<n1; row++)
    {
        for(unsigned int col=0; col<n; col++) lscs.mtrx[row][col] = (rand() % 1000) / 1000.0;
        betta[row+n0] = 0.0;
        for (unsigned int i=0; i<n; i++) betta[row+n0] += lscs.mtrx[row][i] * x(lscs.time, i);
    }
    setLeftSeparatedCondition(lscs);

    Condition rscs;
    rscs.type = SeparatedRight;
    rscs.time = 1.0;
    rscs.nmbr = N;
    rscs.mtrx.resize(n2, n);
    for (unsigned int row=0; row<n2; row++)
    {
        for(unsigned int col=0; col<n; col++) rscs.mtrx[row][col] = (rand() % 1000) / 1000.0;
        betta[row+n0+n1] = 0.0;
        for (unsigned int i=0; i<n; i++) betta[row+n0+n1] += rscs.mtrx[row][i] * x(rscs.time, i);
    }
    setRightSeparatedCondition(rscs);

    setBetta(betta);
}

double SystemLinearODENonLocalContionsV::A(double t UNUSED_PARAM, unsigned int, unsigned int row, unsigned int col) const
{
#ifdef SAMPLE_1
    if (row==0)
    {
        if (col==1) { return t; }
        if (col==2) { return -1.0; }
    }
    if (row==2)
    {
        if (col==1) { return +2.0*t+3.0; }
        if (col==2) { return -2.0; }
    }
#endif

#ifdef SAMPLE_2
    if (row==1)
    {
        if (col==1) { return +2.0; }
        if (col==2) { return +3.0; }
    }
    if (row==2)
    {
        if (col==1) { return +5.0; }
        if (col==2) { return +4.0; }
    }
#endif

#ifdef SAMPLE_3
    if (row==0)
    {
        if (col==0) { return +2.0; }
        if (col==1) { return -3.0; }
        if (col==2) { return +1.0; }
    }
    if (row==1)
    {
        if (col==0) { return +3.0; }
        if (col==1) { return +1.0; }
        if (col==2) { return -2.0; }
    }
    if (row==2)
    {
        if (col==0) { return +1.0; }
        if (col==1) { return -5.0; }
        if (col==2) { return -3.0; }
    }
#endif

    return NAN;
}

double SystemLinearODENonLocalContionsV::B(double t, unsigned int, unsigned int row) const
{
#ifdef SAMPLE_1
    if (row==1) return -t;
    if (row==2) return -6.0*t-11.0;
#endif

#ifdef SAMPLE_2
    if (row==1) return +2.0 - 13.0*t;
    if (row==2) return +3.0 - 22.0*t;
#endif

#ifdef SAMPLE_3
    if (row==0) return 3.0     - (2.0*(3.0*t+4.0) - 3.0*(4.0*t*t) + 1.0*(t*t+t));
    if (row==1) return 8.0*t   - (3.0*(3.0*t+4.0) + 1.0*(4.0*t*t) - 2.0*(t*t+t));
    if (row==2) return 2.0*t+1 - (1.0*(3.0*t+4.0) - 5.0*(4.0*t*t) - 3.0*(t*t+t));
#endif

    return NAN;
}

double SystemLinearODENonLocalContionsV::x(double t, int i) const
{
    //Dimension dim = grid().dimension();
    //double h = dim.step();
    //double t = k*h;

#ifdef SAMPLE_1
    if (i==0) return t*t + 4.0;
    if (i==1) return t*t*t + t;
#endif

#ifdef SAMPLE_2
    if (i==0) return 2.0*t;
    if (i==1) return 3.0*t;
#endif

#ifdef SAMPLE_3
    if (i==0) return 3.0*t+4.0;
    if (i==1) return 4.0*t*t;
    if (i==2) return t*t+t;
#endif
    return NAN;
}


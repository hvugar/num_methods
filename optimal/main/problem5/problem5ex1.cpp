#include "problem5ex1.h"

void Problem5Ex1::Main(int agrc, char *argv[])
{
    ODEGrid grid(Dimension(0.01, 100, 0));
    Problem5Ex1 prob1(grid);
    prob1.initialize();
}

Problem5Ex1::Problem5Ex1(const ODEGrid &grid) : mgrid(grid)
{
}

void Problem5Ex1::initialize()
{
    Dimension dim = mgrid.dimension();
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

    ISystemLinearODENonLocalContionsV::Condition nsc1;
    nsc1.type = ISystemLinearODENonLocalContionsV::NonSeparated;
    nsc1.time = 0.25;
    nsc1.nmbr = N/4;
    nsc1.alpha.resize(n0, n);
    for (unsigned int row=0; row<n0; row++) for(unsigned int col=0; col<n; col++) nsc1.alpha[row][col] = (rand() % 1000) / 1000.0;

    ISystemLinearODENonLocalContionsV::Condition nsc2;
    nsc2.type = ISystemLinearODENonLocalContionsV::NonSeparated;
    nsc2.time = 0.5;
    nsc2.nmbr = N/2;
    nsc2.alpha.resize(n0, n);
    for (unsigned int row=0; row<n0; row++) for(unsigned int col=0; col<n; col++) nsc2.alpha[row][col] = (rand() % 1000) / 1000.0;

    ISystemLinearODENonLocalContionsV::Condition nsc3;
    nsc3.type = ISystemLinearODENonLocalContionsV::NonSeparated;
    nsc3.time = 0.75;
    nsc3.nmbr = 3*(N/4);
    nsc3.alpha.resize(n0, n);
    for (unsigned int row=0; row<n0; row++) for(unsigned int col=0; col<n; col++) nsc3.alpha[row][col] = (rand() % 1000) / 1000.0;

    ISystemLinearODENonLocalContionsV::Condition nsc4;
    nsc4.type = ISystemLinearODENonLocalContionsV::NonSeparated;
    nsc4.time = 1.0;
    nsc4.nmbr = N;
    nsc4.alpha.resize(n0, n);
    for (unsigned int row=0; row<n0; row++) for(unsigned int col=0; col<n; col++) nsc4.alpha[row][col] = (rand() % 1000) / 1000.0;

    //////////////////////////////////////////////////////////////////////////////

    Zetta0 zett0(mgrid, this);
    zett0.addNonSeparatedCondition(nsc0);
    zett0.addNonSeparatedCondition(nsc1);
    zett0.addNonSeparatedCondition(nsc2);
    zett0.addNonSeparatedCondition(nsc3);
    zett0.addNonSeparatedCondition(nsc4);

    DoubleVector betta(n);
    unsigned int L = zett0.nonSeparatedConditions().size();

    for (unsigned int row=0; row<n0; row++)
    {
        betta[row] = 0.0;
        for (unsigned int s=0; s<L; s++)
        {
            const ISystemLinearODENonLocalContionsV::Condition &c = zett0.nonSeparatedConditions().at(s);
            for (unsigned int i=0; i<n; i++) betta[row] += c.alpha[row][i] * X(c.time, 0, i);
        }
    }
    zett0.setBetta(betta);

    DoubleVector z0;
    zett0.calculateForward(z0);

    std::vector<DoubleVector> zm0;
    zett0.calculateBackwardCP(z0,zm0);
    IPrinter::printVector(zm0.at(0));
    IPrinter::printVector(zm0.at(1));
    IPrinter::printVector(zm0.at(2));

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////

//    Zetta1 zett1(mgrid, this);
//    zett1.addNonSeparatedCondition(nsc0);
//    zett1.addNonSeparatedCondition(nsc1);
//    zett1.addNonSeparatedCondition(nsc2);
//    zett1.addNonSeparatedCondition(nsc3);
//    zett1.addNonSeparatedCondition(nsc4);

//    DoubleVector betta1(n, 0.0);
//    zett1.setBetta(betta1);

//    DoubleMatrix zm1;
//    zett1.calculateForward(zm1);

//    IPrinter::print(zm1, zm1.rows(), zm1.cols());
}

double Problem5Ex1::A(double t, unsigned int, unsigned int row, unsigned int col) const
{
#ifdef SAMPLE_1
    if (row == 0) { if (col == 0) { return +2.0; } if (col == 1) { return t; }      if (col == 2) { return -3.0; } }
    if (row == 1) { if (col == 0) { return +3.0; } if (col == 1) { return -4.0*t; } if (col == 2) { return -8.0; } }
    if (row == 2) { if (col == 0) { return +t; }   if (col == 1) { return +1.0; }   if (col == 2) { return -1.0; } }
#endif
    return NAN;
}

double Problem5Ex1::B(double t, unsigned int k, unsigned int row) const
{
#ifdef SAMPLE_1
    if (row == 0) { return dX(t,k,0) - (A(t,k,0,0)*X(t,k,0) + A(t,k,0,1)*X(t,k,1) + A(t,k,0,2)*X(t,k,2))
                - (C(t,k,0,0,0)*g(t,k,0,0) + C(t,k,0,1,0)*g(t,k,1,0) + C(t,k,0,2,0)*g(t,k,2,0))
                - (C(t,k,0,0,1)*g(t,k,0,1) + C(t,k,0,1,1)*g(t,k,1,1) + C(t,k,0,2,1)*g(t,k,2,1)); }
    if (row == 1) { return dX(t,k,1) - (A(t,k,1,0)*X(t,k,0) + A(t,k,1,1)*X(t,k,1) + A(t,k,1,2)*X(t,k,2))
                - (C(t,k,1,0,0)*g(t,k,0,0) + C(t,k,1,1,0)*g(t,k,1,0) + C(t,k,1,2,0)*g(t,k,2,0))
                - (C(t,k,1,0,1)*g(t,k,0,1) + C(t,k,1,1,1)*g(t,k,1,1) + C(t,k,1,2,1)*g(t,k,2,1)); }
    if (row == 2) { return dX(t,k,2) - (A(t,k,2,0)*X(t,k,0) + A(t,k,2,1)*X(t,k,1) + A(t,k,2,2)*X(t,k,2))
                - (C(t,k,2,0,0)*g(t,k,0,0) + C(t,k,2,1,0)*g(t,k,1,0) + C(t,k,2,2,0)*g(t,k,2,0))
                - (C(t,k,2,0,1)*g(t,k,0,1) + C(t,k,2,1,1)*g(t,k,1,1) + C(t,k,2,2,1)*g(t,k,2,1)); }
#endif
    return NAN;
}

double Problem5Ex1::C(double, unsigned int, unsigned int row, unsigned int col, unsigned int i) const
{
#ifdef SAMPLE_1
    if ( i == 0 )
    {
        if (row == 0) { if (col == 0) { return +2.0; } if (col == 1) { return +5.0; } if (col == 2) { return +3.0; } }
        if (row == 1) { if (col == 0) { return +4.0; } if (col == 1) { return +8.0; } if (col == 2) { return +1.0; } }
        if (row == 2) { if (col == 0) { return +1.0; } if (col == 1) { return +3.0; } if (col == 2) { return +4.0; } }
    }
    if ( i == 1 )
    {
        if (row == 0) { if (col == 0) { return +1.0; } if (col == 1) { return +3.0; } if (col == 2) { return +4.0; } }
        if (row == 1) { if (col == 0) { return +2.0; } if (col == 1) { return +3.0; } if (col == 2) { return +1.0; } }
        if (row == 2) { if (col == 0) { return +5.0; } if (col == 1) { return +2.0; } if (col == 2) { return +8.0; } }
    }
#endif
    return NAN;
}

double Problem5Ex1::g(double t, unsigned int k, unsigned int row, unsigned int i) const
{
#ifdef SAMPLE_1
    if (i == 0)
    {
        if (row == 0) { return X(0.3,k,row)*X(0.3,k,row); }
        if (row == 1) { return X(0.3,k,row)*X(0.3,k,2)*X(0.3,k,row); }
        if (row == 2) { return X(0.3,k,row); }
    }
    if (i == 1)
    {
        if (row == 0) { return X(0.6,k,row)*X(0.6,k,row)*X(0.6,k,row); }
        if (row == 1) { return X(0.6,k,row); }
        if (row == 2) { return X(0.6,k,row)*X(0.6,k,row); }
    }
#endif
    return NAN;
}

double Problem5Ex1::X(double t, unsigned int, unsigned int i) const
{
#ifdef SAMPLE_1
    if (i == 0) { return t*t+t+2.0; }
    if (i == 1) { return 2.0*t-3.0; }
    if (i == 2) { return t*t*t+t; }
#endif
    return NAN;
}

double Problem5Ex1::dX(double t, unsigned int, unsigned int i) const
{
#ifdef SAMPLE_1
    if (i == 0) { return 2.0*t+1.0; }
    if (i == 1) { return 2.0; }
    if (i == 2) { return 3.0*t*t+1.0; }
#endif
    return NAN;
}

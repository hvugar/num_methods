#include "problem5ex1.h"

void Problem5Ex1::Main(int agrc, char *argv[])
{
    ODEGrid grid(Dimension(0.01, 100, 0));
    Problem5Ex1 prob1(grid);
    prob1.initialize();
}

Problem5Ex1::Problem5Ex1(const ODEGrid &grid) : mgrid(grid)
{
//    printf("%14.10f\n", g(0,0));
//    printf("%14.10f\n", g(0,1));
//    printf("%14.10f\n", g(0,2));

//    printf("%14.10f\n", g(1,0));
//    printf("%14.10f\n", g(1,1));
//    printf("%14.10f\n", g(1,2));
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
    nsc3.nmbr = 3*N/4;
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
            for (unsigned int i=0; i<n; i++) betta[row] += c.alpha[row][i] * X(c.time, i);
        }
    }
    zett0.setBetta(betta);

    DoubleVector z0;
    std::vector<DoubleVector> zm0;
    zett0.calculateForward(z0);
    zett0.calculateBackwardCP(z0,zm0);

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////

    Zetta1 zett1(mgrid, this);
    zett1.addNonSeparatedCondition(nsc0);
    zett1.addNonSeparatedCondition(nsc1);
    zett1.addNonSeparatedCondition(nsc2);
    zett1.addNonSeparatedCondition(nsc3);
    zett1.addNonSeparatedCondition(nsc4);
    DoubleVector betta1(n, 0.0);
    zett1.setBetta(betta1);

    DoubleMatrix z1;
    std::vector<std::vector<DoubleVector>> zm1;
    zett1.calculateForward(z1);
    zett1.calculateBackwardCP(z1,zm1);

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////

    Zetta2 zett2(mgrid, this);
    zett2.addNonSeparatedCondition(nsc0);
    zett2.addNonSeparatedCondition(nsc1);
    zett2.addNonSeparatedCondition(nsc2);
    zett2.addNonSeparatedCondition(nsc3);
    zett2.addNonSeparatedCondition(nsc4);
    DoubleVector betta2(n, 0.0);
    zett2.setBetta(betta2);

    DoubleMatrix z2;
    std::vector<std::vector<DoubleVector>> zm2;
    zett2.calculateForward(z2);
    zett2.calculateBackwardCP(z2,zm2);

    DoubleVector x1(N+1);
    DoubleVector x2(N+1);
    DoubleVector x3(N+1);
    for (unsigned int i=0; i<=N; i++)
    {
        x1[i] =  zm0[0][i] + (zm1[0][0][i]*g(0,0) + zm1[0][1][i]*g(0,1) + zm1[0][2][i]*g(0,2)) + (zm2[0][0][i]*g(1,0) + zm2[0][1][i]*g(1,1) + zm2[0][2][i]*g(1,2));
        x2[i] =  zm0[1][i] + (zm1[1][0][i]*g(0,0) + zm1[1][1][i]*g(0,1) + zm1[1][2][i]*g(0,2)) + (zm2[1][0][i]*g(1,0) + zm2[1][1][i]*g(1,1) + zm2[1][2][i]*g(1,2));
        x3[i] =  zm0[2][i] + (zm1[2][0][i]*g(0,0) + zm1[2][1][i]*g(0,1) + zm1[2][2][i]*g(0,2)) + (zm2[2][0][i]*g(1,0) + zm2[2][1][i]*g(1,1) + zm2[2][2][i]*g(1,2));
    }
    IPrinter::printVector(x1);
    IPrinter::printVector(x2);
    IPrinter::printVector(x3);

//    FILE *file = fopen("temp.txt", "w");
//    IPrinter::printVector(14, 10, zm0.at(0),NULL,N,0,0,file);
//    IPrinter::printVector(14, 10, zm0.at(1),NULL,N,0,0,file);
//    IPrinter::printVector(14, 10, zm0.at(2),NULL,N,0,0,file);

//    IPrinter::printVector(14, 10, zm1.at(0).at(0),NULL,N,0,0,file);
//    IPrinter::printVector(14, 10, zm1.at(0).at(1),NULL,N,0,0,file);
//    IPrinter::printVector(14, 10, zm1.at(0).at(2),NULL,N,0,0,file);

//    IPrinter::printVector(14, 10, zm1.at(1).at(0),NULL,N,0,0,file);
//    IPrinter::printVector(14, 10, zm1.at(1).at(1),NULL,N,0,0,file);
//    IPrinter::printVector(14, 10, zm1.at(1).at(2),NULL,N,0,0,file);

//    IPrinter::printVector(14, 10, zm1.at(2).at(0),NULL,N,0,0,file);
//    IPrinter::printVector(14, 10, zm1.at(2).at(1),NULL,N,0,0,file);
//    IPrinter::printVector(14, 10, zm1.at(2).at(2),NULL,N,0,0,file);

//    IPrinter::printVector(14, 10, zm2.at(0).at(0),NULL,N,0,0,file);
//    IPrinter::printVector(14, 10, zm2.at(0).at(1),NULL,N,0,0,file);
//    IPrinter::printVector(14, 10, zm2.at(0).at(2),NULL,N,0,0,file);

//    IPrinter::printVector(14, 10, zm2.at(1).at(0),NULL,N,0,0,file);
//    IPrinter::printVector(14, 10, zm2.at(1).at(1),NULL,N,0,0,file);
//    IPrinter::printVector(14, 10, zm2.at(1).at(2),NULL,N,0,0,file);

//    IPrinter::printVector(14, 10, zm2.at(2).at(0),NULL,N,0,0,file);
//    IPrinter::printVector(14, 10, zm2.at(2).at(1),NULL,N,0,0,file);
//    IPrinter::printVector(14, 10, zm2.at(2).at(2),NULL,N,0,0,file);

//    fclose(file);
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
    if (row == 0)
    {
        return dX(t,0)
                - (A(t,k,0,0)*X(t,0) + A(t,k,0,1)*X(t,1) + A(t,k,0,2)*X(t,2))
                - (C(t,k,0,0,0)*g(0,0) + C(t,k,0,0,1)*g(0,1) + C(t,k,0,0,2)*g(0,2)) - (C(t,k,1,0,0)*g(1,0) + C(t,k,1,0,1)*g(1,1) + C(t,k,1,0,2)*g(1,2));
    }
    if (row == 1)
    {
        return dX(t,1)
                - (A(t,k,1,0)*X(t,0) + A(t,k,1,1)*X(t,1) + A(t,k,1,2)*X(t,2))
                - (C(t,k,0,1,0)*g(0,0) + C(t,k,0,1,1)*g(0,1) + C(t,k,0,1,2)*g(0,2)) - (C(t,k,1,1,0)*g(1,0) + C(t,k,1,1,1)*g(1,1) + C(t,k,1,1,2)*g(1,2));
    }
    if (row == 2)
    { return dX(t,2)
                - (A(t,k,2,0)*X(t,0) + A(t,k,2,1)*X(t,1) + A(t,k,2,2)*X(t,2))
                - (C(t,k,0,2,0)*g(0,0) + C(t,k,0,2,1)*g(0,1) + C(t,k,0,2,2)*g(0,2)) - (C(t,k,1,2,0)*g(1,0) + C(t,k,1,2,1)*g(1,1) + C(t,k,1,2,2)*g(1,2));
    }
#endif
    return NAN;
}

double Problem5Ex1::C(double, unsigned int, unsigned int num, unsigned int row, unsigned int col) const
{
#ifdef SAMPLE_1
    if ( num == 0 )
    {
        if (row == 0) { if (col == 0) { return +2.0; } if (col == 1) { return +5.0; } if (col == 2) { return +3.0; } }
        if (row == 1) { if (col == 0) { return +4.0; } if (col == 1) { return +8.0; } if (col == 2) { return +1.0; } }
        if (row == 2) { if (col == 0) { return +1.0; } if (col == 1) { return +3.0; } if (col == 2) { return +4.0; } }
    }
    if ( num == 1 )
    {
        if (row == 0) { if (col == 0) { return +1.0; } if (col == 1) { return +3.0; } if (col == 2) { return +4.0; } }
        if (row == 1) { if (col == 0) { return +2.0; } if (col == 1) { return +3.0; } if (col == 2) { return +1.0; } }
        if (row == 2) { if (col == 0) { return +5.0; } if (col == 1) { return +2.0; } if (col == 2) { return +8.0; } }
    }
#endif
    return NAN;
}

double Problem5Ex1::g(unsigned int num, unsigned int row) const
{
   // return 0.0;
#ifdef SAMPLE_1
//    if (num == 0)
//    {
//        if (row == 0) { return X(0.3, 0)*X(0.3,0) + 2.0*X(0.3,0); }
//        if (row == 1) { return X(0.3, 1)*X(0.3,1) - X(0.3,1); }
//        if (row == 2) { return X(0.3, 2); }
//    }
//    if (num == 1)
//    {
//        if (row == 0) { return X(0.6,0)*X(0.6,0) + X(0.6,1); }
//        if (row == 1) { return X(0.6,1); }
//        if (row == 2) { return X(0.6,2)*X(0.6,2); }
//    }

        if (num == 0)
        {
            if (row == 0) { return sin(X(0.3, 0)); }
            if (row == 1) { return cos(X(0.3, 1)); }
            if (row == 2) { return tan(X(0.3, 2)); }
        }
        if (num == 1)
        {
            if (row == 0) { return sin(2.0*X(0.6,0)); }
            if (row == 1) { return sin(9.0*X(0.6,0)); }
            if (row == 2) { return cos(X(0.6,2)); }
        }


//    if (num == 0)
//    {
//        if (row == 0) { return 10.4921000000; }
//        if (row == 1) { return  8.1600000000; }
//        if (row == 2) { return  0.3270000000; }
//    }
//    if (num == 1)
//    {
//        if (row == 0) { return  6.9616000000; }
//        if (row == 1) { return -1.8000000000; }
//        if (row == 2) { return  0.6658560000; }
//    }

#endif
    return NAN;
}

double Problem5Ex1::X(double t, unsigned int i) const
{
#ifdef SAMPLE_1
    if (i == 0) { return t*t+t+2.0; }
    if (i == 1) { return 2.0*t-3.0; }
    if (i == 2) { return t*t*t+t; }
#endif
    return NAN;
}

double Problem5Ex1::dX(double t, unsigned int i) const
{
#ifdef SAMPLE_1
    if (i == 0) { return 2.0*t+1.0; }
    if (i == 1) { return 2.0; }
    if (i == 2) { return 3.0*t*t+1.0; }
#endif
    return NAN;
}

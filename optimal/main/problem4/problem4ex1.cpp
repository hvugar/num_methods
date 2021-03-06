#include "problem4ex1.h"
#include <math.h>

void Problem4Ex1::Main(int agrc UNUSED_PARAM, char *argv[] UNUSED_PARAM)
{
    UniformODEGrid grid(0.001, 0, 1000);
    Problem4Ex1 prob1(grid);
    prob1.initialize();

    DoubleVector x0;
    x0 << +0.0 << +0.0 << +0.0 << +0.0 << +0.1 << +0.1;
    DoubleVector x;
    //    prob1.calculateSimpleIdetartion(x0, x, 0.001);
    prob1.calculateNewtonMethodMod(x0, x, 0.001, 0.001);

    //    IPrinter::printSeperatorLine();
    //    prob1.printResult();
    IPrinter::print(x,x.length());
    IPrinter::printSeperatorLine();
    prob1.printResult1(x);
}

Problem4Ex1::Problem4Ex1(const UniformODEGrid &grid) : mgrid(grid)
{
}

void Problem4Ex1::initialize()
{
    unsigned int N = mgrid.sizeN();

    unsigned int n0 = 3;
    unsigned int n1 = 0;
    unsigned int n2 = 0;
    unsigned int n = n0 + n1 + n2;

    //    Condition nsc0;
    //    nsc0.type = ISystemLinearODENonLocalContionsV::NonSeparated;
    //    nsc0.time = 0.0;
    //    nsc0.nmbr = 0;
    //    nsc0.mtrx.resize(n0, n);
    //    for (unsigned int row=0; row<n0; row++) for(unsigned int col=0; col<n; col++) nsc0.mtrx[row][col] = (rand() % 10) / 1.0;

    //    Condition nsc1;
    //    nsc1.type = ISystemLinearODENonLocalContionsV::NonSeparated;
    //    nsc1.time = 0.2;
    //    nsc1.nmbr = N/5;
    //    nsc1.mtrx.resize(n0, n);
    //    for (unsigned int row=0; row<n0; row++) for(unsigned int col=0; col<n; col++) nsc1.mtrx[row][col] = (rand() % 10) / 1.0;

    //    Condition nsc2;
    //    nsc2.type = ISystemLinearODENonLocalContionsV::NonSeparated;
    //    nsc2.time = 0.5;
    //    nsc2.nmbr = N/2;
    //    nsc2.mtrx.resize(n0, n);
    //    for (unsigned int row=0; row<n0; row++) for(unsigned int col=0; col<n; col++) nsc2.mtrx[row][col] = (rand() % 10) / 1.0;

    //    Condition nsc3;
    //    nsc3.type = ISystemLinearODENonLocalContionsV::NonSeparated;
    //    nsc3.time = 0.8;
    //    nsc3.nmbr = 4*(N/5);
    //    nsc3.mtrx.resize(n0, n);
    //    for (unsigned int row=0; row<n0; row++) for(unsigned int col=0; col<n; col++) nsc3.mtrx[row][col] = (rand() % 10) / 1.0;

    //    Condition nsc4;
    //    nsc4.type = ISystemLinearODENonLocalContionsV::NonSeparated;
    //    nsc4.time = 1.0;
    //    nsc4.nmbr = N;
    //    nsc4.mtrx.resize(n0, n);
    //    for (unsigned int row=0; row<n0; row++) for(unsigned int col=0; col<n; col++) nsc4.mtrx[row][col] = (rand() % 10) / 1.0;

    //////////////////////////////////////////////////////////////////////////////


    LinearODE1stOrder::Condition c0;
    c0.time = 0.0;
    c0.nmbr = 0;
    c0.mtrx.resize(n0, n);
    //    c0.mtrx = nsc0.mtrx;
    //    for (unsigned int row=0; row<n0; row++) for(unsigned int col=0; col<n; col++) c0.mtrx[row][col] = (rand() % 10) / 1.0;
    c0.mtrx[0][0] = 1.0000000000; c0.mtrx[0][1] = 7.0000000000; c0.mtrx[0][2] = 4.0000000000;
    c0.mtrx[1][0] = 0.0000000000; c0.mtrx[1][1] = 9.0000000000; c0.mtrx[1][2] = 4.0000000000;
    c0.mtrx[2][0] = 8.0000000000; c0.mtrx[2][1] = 8.0000000000; c0.mtrx[2][2] = 2.0000000000;

    LinearODE1stOrder::Condition c1;
    c1.time = 0.2;
    c1.nmbr = N/5;
    c1.mtrx.resize(n0, n);
    //    c1.mtrx = nsc1.mtrx;
    //    for (unsigned int row=0; row<n0; row++) for(unsigned int col=0; col<n; col++) c1.mtrx[row][col] = (rand() % 10) / 1.0;
    c1.mtrx[0][0] = 4.0000000000; c1.mtrx[0][1] = 5.0000000000; c1.mtrx[0][2] = 5.0000000000;
    c1.mtrx[1][0] = 1.0000000000; c1.mtrx[1][1] = 7.0000000000; c1.mtrx[1][2] = 2.0000000000;
    c1.mtrx[2][0] = 1.0000000000; c1.mtrx[2][1] = 5.0000000000; c1.mtrx[2][2] = 2.0000000000;

    LinearODE1stOrder::Condition c2;
    c2.time = 0.5;
    c2.nmbr = N/2;
    c2.mtrx.resize(n0, n);
    //    c2.mtrx = nsc2.mtrx;
    //    for (unsigned int row=0; row<n0; row++) for(unsigned int col=0; col<n; col++) c2.mtrx[row][col] = (rand() % 10) / 1.0;
    c2.mtrx[0][0] = 7.0000000000; c2.mtrx[0][1] = 6.0000000000; c2.mtrx[0][2] = 1.0000000000;
    c2.mtrx[1][0] = 4.0000000000; c2.mtrx[1][1] = 2.0000000000; c2.mtrx[1][2] = 3.0000000000;
    c2.mtrx[2][0] = 2.0000000000; c2.mtrx[2][1] = 2.0000000000; c2.mtrx[2][2] = 1.0000000000;

    LinearODE1stOrder::Condition c3;
    c3.time = 0.8;
    c3.nmbr = 4*(N/5);
    c3.mtrx.resize(n0, n);
    //    c3.mtrx = nsc3.mtrx;
    //    for (unsigned int row=0; row<n0; row++) for(unsigned int col=0; col<n; col++) c3.mtrx[row][col] = (rand() % 10) / 1.0;
    c3.mtrx[0][0] = 6.0000000000; c3.mtrx[0][1] = 8.0000000000; c3.mtrx[0][2] = 5.0000000000;
    c3.mtrx[1][0] = 7.0000000000; c3.mtrx[1][1] = 6.0000000000; c3.mtrx[1][2] = 1.0000000000;
    c3.mtrx[2][0] = 8.0000000000; c3.mtrx[2][1] = 9.0000000000; c3.mtrx[2][2] = 2.0000000000;

    LinearODE1stOrder::Condition c4;
    c4.time = 1.0;
    c4.nmbr = N;
    c4.mtrx.resize(n0, n);
    //    c4.mtrx = nsc4.mtrx;
    //    for (unsigned int row=0; row<n0; row++) for(unsigned int col=0; col<n; col++) c4.mtrx[row][col] = (rand() % 10) / 1.0;
    c4.mtrx[0][0] = 7.0000000000; c4.mtrx[0][1] = 9.0000000000; c4.mtrx[0][2] = 5.0000000000;
    c4.mtrx[1][0] = 4.0000000000; c4.mtrx[1][1] = 3.0000000000; c4.mtrx[1][2] = 1.0000000000;
    c4.mtrx[2][0] = 2.0000000000; c4.mtrx[2][1] = 3.0000000000; c4.mtrx[2][2] = 3.0000000000;

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    //    Zetta0 zett0(*this);
    //    zett0.setGrid(mgrid);
    //    zett0.addNonSeparatedCondition(nsc0);
    //    zett0.addNonSeparatedCondition(nsc1);
    //    zett0.addNonSeparatedCondition(nsc2);
    //    zett0.addNonSeparatedCondition(nsc3);
    //    zett0.addNonSeparatedCondition(nsc4);

    //    DoubleVector betta(n);
    //    unsigned int L = zett0.nonSeparatedConditions().size();

    //    for (unsigned int row=0; row<n0; row++)
    //    {
    //        betta[row] = 0.0;
    //        for (unsigned int s=0; s<L; s++)
    //        {
    //            const Condition &c = zett0.nonSeparatedConditions().at(s);
    //            for (unsigned int i=0; i<n; i++) betta[row] += c.mtrx[row][i] * X(c.time, i);
    //        }
    //    }
    //    zett0.setBetta(betta);

    //    DoubleVector z0;
    //    zett0.calculateForward(z0);
    //    zett0.calculateBackwardCP(z0,zm0);

    Zetta01 zett01(*this);
    zett01.setGrid(mgrid);
    std::vector<LinearODE1stOrder::Condition> cs;
    cs.push_back(c0);
    cs.push_back(c1);
    cs.push_back(c2);
    cs.push_back(c3);
    cs.push_back(c4);

    DoubleVector betta(n);
    unsigned int L = cs.size();

    for (unsigned int row=0; row<n0; row++)
    {
        betta[row] = 0.0;
        for (unsigned int s=0; s<L; s++)
        {
            const LinearODE1stOrder::Condition &c = cs.at(s);
            for (unsigned int i=0; i<n; i++) betta[row] += c.mtrx[row][i] * X(c.time, i);
        }
    }

    zett01.calculate(cs, betta, zm0);

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////

    //    Zettai zett1(*this, 0);
    //    zett1.setGrid(mgrid);
    //    zett1.addNonSeparatedCondition(nsc0);
    //    zett1.addNonSeparatedCondition(nsc1);
    //    zett1.addNonSeparatedCondition(nsc2);
    //    zett1.addNonSeparatedCondition(nsc3);
    //    zett1.addNonSeparatedCondition(nsc4);
    //    DoubleVector betta1(n, 0.0);
    //    zett1.setBetta(betta1);

    //    DoubleMatrix z1;
    //    zett1.calculateForward(z1);
    //    zett1.calculateBackwardCP(z1,zm1);


    ///////////////////////////////////////////////////////////////////////////////////////////////////////////

    //    Zettai zett2(*this, 1);
    //    zett2.setGrid(mgrid);
    //    zett2.addNonSeparatedCondition(nsc0);
    //    zett2.addNonSeparatedCondition(nsc1);
    //    zett2.addNonSeparatedCondition(nsc2);
    //    zett2.addNonSeparatedCondition(nsc3);
    //    zett2.addNonSeparatedCondition(nsc4);
    //    DoubleVector betta2(n, 0.0);
    //    zett2.setBetta(betta2);

    //    DoubleMatrix z2;
    //    zett2.calculateForward(z2);
    //    zett2.calculateBackwardCP(z2,zm2);

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////

    DoubleMatrix betta0(n, n, 0.0);

    Zettai1 zett1(*this, 0);
    zett1.setGrid(mgrid);
    zett1.calculateM(cs, betta0, zm1);

    Zettai1 zett2(*this, 1);
    zett2.setGrid(mgrid);
    zett2.calculateM(cs, betta0, zm2);
}

void Problem4Ex1::printResult()
{
    unsigned int N = mgrid.sizeN();

    DoubleVector x1(N+1);
    DoubleVector x2(N+1);
    DoubleVector x3(N+1);
    for (unsigned int i=0; i<=N; i++)
    {
        x1[i] =  zm0[0][i] + (zm1[0][0][i]*g(0,0) + zm1[0][1][i]*g(0,1) + zm1[0][2][i]*g(0,2)) + (zm2[0][0][i]*g(1,0) + zm2[0][1][i]*g(1,1) + zm2[0][2][i]*g(1,2));
        x2[i] =  zm0[1][i] + (zm1[1][0][i]*g(0,0) + zm1[1][1][i]*g(0,1) + zm1[1][2][i]*g(0,2)) + (zm2[1][0][i]*g(1,0) + zm2[1][1][i]*g(1,1) + zm2[1][2][i]*g(1,2));
        x3[i] =  zm0[2][i] + (zm1[2][0][i]*g(0,0) + zm1[2][1][i]*g(0,1) + zm1[2][2][i]*g(0,2)) + (zm2[2][0][i]*g(1,0) + zm2[2][1][i]*g(1,1) + zm2[2][2][i]*g(1,2));
    }
    IPrinter::printVector(14,10,x1);
    IPrinter::printVector(14,10,x2);
    IPrinter::printVector(14,10,x3);
}

void Problem4Ex1::printResult1(const DoubleVector &x)
{
    unsigned int N = mgrid.sizeN();

    DoubleVector x1(N+1);
    DoubleVector x2(N+1);
    DoubleVector x3(N+1);
    for (unsigned int i=0; i<=N; i++)
    {
        x1[i] =  zm0[0][i] + (zm1[0][0][i]*g(x,0,0) + zm1[0][1][i]*g(x,0,1) + zm1[0][2][i]*g(x,0,2)) + (zm2[0][0][i]*g(x,1,0) + zm2[0][1][i]*g(x,1,1) + zm2[0][2][i]*g(x,1,2));
        x2[i] =  zm0[1][i] + (zm1[1][0][i]*g(x,0,0) + zm1[1][1][i]*g(x,0,1) + zm1[1][2][i]*g(x,0,2)) + (zm2[1][0][i]*g(x,1,0) + zm2[1][1][i]*g(x,1,1) + zm2[1][2][i]*g(x,1,2));
        x3[i] =  zm0[2][i] + (zm1[2][0][i]*g(x,0,0) + zm1[2][1][i]*g(x,0,1) + zm1[2][2][i]*g(x,0,2)) + (zm2[2][0][i]*g(x,1,0) + zm2[2][1][i]*g(x,1,1) + zm2[2][2][i]*g(x,1,2));
    }
    IPrinter::printVector(14,10,x1);
    IPrinter::printVector(14,10,x2);
    IPrinter::printVector(14,10,x3);
}

double Problem4Ex1::fx(const DoubleVector &x, unsigned int num) const
{
    unsigned int N = mgrid.sizeN();

    if (num == 0)
    {
        unsigned int n = 0.3*N;
        return zm0[0][n] + (zm1[0][0][n]*g(x,0,0) + zm1[0][1][n]*g(x,0,1) + zm1[0][2][n]*g(x,0,2))
                + (zm2[0][0][n]*g(x,1,0) + zm2[0][1][n]*g(x,1,1) + zm2[0][2][n]*g(x,1,2)) - x[0];
    }
    if (num == 1)
    {
        unsigned int n = 0.3*N;
        return zm0[1][n] + (zm1[1][0][n]*g(x,0,0) + zm1[1][1][n]*g(x,0,1) + zm1[1][2][n]*g(x,0,2))
                + (zm2[1][0][n]*g(x,1,0) + zm2[1][1][n]*g(x,1,1) + zm2[1][2][n]*g(x,1,2)) - x[1];
    }
    if (num == 2)
    {
        unsigned int n = 0.3*N;
        return zm0[2][n] + (zm1[2][0][n]*g(x,0,0) + zm1[2][1][n]*g(x,0,1) + zm1[2][2][n]*g(x,0,2))
                + (zm2[2][0][n]*g(x,1,0) + zm2[2][1][n]*g(x,1,1) + zm2[2][2][n]*g(x,1,2)) - x[2];
    }

    if (num == 3)
    {
        unsigned int n = 0.6*N;
        return zm0[0][n] + (zm1[0][0][n]*g(x,0,0) + zm1[0][1][n]*g(x,0,1) + zm1[0][2][n]*g(x,0,2))
                + (zm2[0][0][n]*g(x,1,0) + zm2[0][1][n]*g(x,1,1) + zm2[0][2][n]*g(x,1,2)) - x[3];
    }
    if (num == 4)
    {
        unsigned int n = 0.6*N;
        return zm0[1][n] + (zm1[1][0][n]*g(x,0,0) + zm1[1][1][n]*g(x,0,1) + zm1[1][2][n]*g(x,0,2))
                + (zm2[1][0][n]*g(x,1,0) + zm2[1][1][n]*g(x,1,1) + zm2[1][2][n]*g(x,1,2)) - x[4];
    }
    if (num == 5)
    {
        unsigned int n = 0.6*N;
        return zm0[2][n] + (zm1[2][0][n]*g(x,0,0) + zm1[2][1][n]*g(x,0,1) + zm1[2][2][n]*g(x,0,2))
                + (zm2[2][0][n]*g(x,1,0) + zm2[2][1][n]*g(x,1,1) + zm2[2][2][n]*g(x,1,2)) - x[5];
    }
    return NAN;
}

double Problem4Ex1::A(const GridNodeODE &node, unsigned int row, unsigned int col) const
{
#ifdef SAMPLE_1
    double t = node.x;
    if (row == 0) { if (col == 0) { return +2.0; } if (col == 1) { return t; }      if (col == 2) { return -3.0; } }
    if (row == 1) { if (col == 0) { return +3.0; } if (col == 1) { return -4.0*t; } if (col == 2) { return -8.0; } }
    if (row == 2) { if (col == 0) { return t; }    if (col == 1) { return +1.0; }   if (col == 2) { return -1.0; } }
#endif
#ifdef SAMPLE_2
    if (row == 0) { if (col == 0) return +1.0; if (col == 0) return +3.0; }
    if (row == 1) { if (col == 0) return +2.0; if (col == 0) return +4.0; }
#endif
    return NAN;
}

double Problem4Ex1::C(const GridNodeODE &node UNUSED_PARAM, unsigned int num, unsigned int row, unsigned int col) const
{
    double K = 10.0;
#ifdef SAMPLE_1

    if ( num == 0 )
    {
        if (row == 0) { if (col == 0) { return +0.002*K; } if (col == 1) { return +0.005*K; } if (col == 2) { return +0.003*K; } }
        if (row == 1) { if (col == 0) { return +0.004*K; } if (col == 1) { return +0.008*K; } if (col == 2) { return +0.001*K; } }
        if (row == 2) { if (col == 0) { return +0.001*K; } if (col == 1) { return +0.003*K; } if (col == 2) { return +0.004*K; } }
    }

    if ( num == 1 )
    {
        if (row == 0) { if (col == 0) { return +0.001*K; } if (col == 1) { return +0.003*K; } if (col == 2) { return +0.004*K; } }
        if (row == 1) { if (col == 0) { return +0.002*K; } if (col == 1) { return +0.005*K; } if (col == 2) { return +0.001*K; } }
        if (row == 2) { if (col == 0) { return +0.005*K; } if (col == 1) { return +0.002*K; } if (col == 2) { return +0.008*K; } }
    }

#endif
#ifdef SAMPLE_2
    if ( num == 0 )
    {
        if (row == 0) { if (col == 0) { return +0.002; } if (col == 1) { return +0.005; } }
        if (row == 1) { if (col == 0) { return +0.004; } if (col == 1) { return +0.008; } }
    }
#endif
    return NAN;
}

double Problem4Ex1::B(const GridNodeODE &node UNUSED_PARAM, unsigned int row) const
{
    double K = 10.0;
    double t = node.x;
#ifdef SAMPLE_1
    if (row == 0) return 3.0*t*t*t - 4.0*t*t + 6.0*t - 3.0 - 0.0414023610*K;
    if (row == 1) return 8.0*t*t*t + 5.0*t*t - 7.0*t - 4.0 - 0.0735341450*K;
    if (row == 2) return 2.0*t*t - 3.0*t             + 4.0 - 0.0639390740*K;
#endif
#ifdef SAMPLE_2
    if (row == 0) return 4.0*cos(4.0*t) - sin(4.0*t) - 3.0*(t*t-t) - 0.0069072257;
    if (row == 1) return (2.0*t-1)-2.0*sin(4.0*t) - 4.0*(t*t-t) - 0.0115811592;
#endif
    return NAN;
}

double Problem4Ex1::g(unsigned int num, unsigned int row) const
{
    //    return 0.0;
#ifdef SAMPLE_1
    if (num == 0)
    {
        if (row == 0) { return X(0.3, 0)*X(0.3,0) + X(0.3, 1)*X(0.3,1) - X(0.3,2); }
        if (row == 1) { return X(0.3, 1)*X(0.3,1) - X(0.3, 0)*X(0.3,2); }
        if (row == 2) { return X(0.3, 0) + X(0.3, 1) + X(0.3, 2)*X(0.3, 2); }
    }
    if (num == 1)
    {
        if (row == 0) { return X(0.6,0)*X(0.6,0); }
        if (row == 1) { return X(0.6,1)*X(0.6,1)*X(0.6,1); }
        if (row == 2) { return X(0.6,2)*X(0.6,2); }
    }
#endif

#ifdef SAMPLE_2
    if (num == 0)
    {
        DoubleVector x;
        x << X(0.5, 0) << X(0.5, 1);
        return g(x, num, row);
    }
#endif

    return NAN;
}

double Problem4Ex1::g(const DoubleVector &x, unsigned int num, unsigned int row) const
{
#ifdef SAMPLE_1
    double x1_03 = x[0]; double x1_06 = x[3];
    double x2_03 = x[1]; double x2_06 = x[4];
    double x3_03 = x[2]; double x3_06 = x[5];

    if (num == 0)
    {
        if (row == 0) return x1_03*x1_03 + x2_03*x2_03 - x3_03;
        if (row == 1) return x2_03*x2_03 - x1_03*x3_03;
        if (row == 2) return x1_03 + x2_03 + x3_03*x3_03;
    }
    if (num == 1)
    {
        if (row == 0) { return  x1_06*x1_06; }
        if (row == 1) return x2_06*x2_06*x2_06;
        if (row == 2) return x3_06*x3_06;
    }
#endif
#ifdef SAMPLE_2
    double x1_05 = x[0];
    double x2_05 = x[1];
    if (num == 0)
    {
        if (row == 0) return x1_05*x1_05 + x2_05*x2_05 + x1_05*x2_05;
        if (row == 1) return x1_05*x1_05 + x2_05*x2_05 - x1_05*x2_05;
    }
#endif

    return NAN;
}

double Problem4Ex1::X(double t, unsigned int i) const
{
#ifdef SAMPLE_1
    if (i == 0) { return t*t+t+2.0; }
    if (i == 1) { return 2.0*t-3.0; }
    if (i == 2) { return t*t*t+t; }
#endif
#ifdef SAMPLE_2
    if (i == 0) { return sin(4.0*t); }
    if (i == 1) { return t*t-t; }
#endif
    return NAN;
}

double Problem4Ex1::dX(double t, unsigned int i) const
{
#ifdef SAMPLE_1
    if (i == 0) { return 2.0*t+1.0; }
    if (i == 1) { return 2.0; }
    if (i == 2) { return 3.0*t*t+1.0; }
#endif
    return NAN;
}

#include "problem4ex1.h"

void Problem4Ex1::Main(int agrc, char *argv[])
{
    ODEGrid grid(Dimension(0.1, 10, 0));
    Problem4Ex1 prob1(grid);
    prob1.initialize();

    for (unsigned int i=0; i<=grid.dimension().sizeN(); i++)
    {
        printf("%10.6f", prob1.X(i*grid.dimension().step(),2));
    }
    puts("");

//    DoubleVector x0;
//    x0 << 2.3 << -2.4 << 0.3 << 2.9 << -1.8 << 0.8;
//    DoubleVector x;
//    prob1.calculate(x0, x, 0.001);

//    printf("%14.10f %14.10f %14.10f\n", x[0], x[1], x[2]);
//    printf("%14.10f %14.10f %14.10f\n", x[3], x[4], x[5]);

//    ---------------------------------------------------------------------------------------------------------------------------------------------------------------------
//      2.0140283546   2.1215738596   2.2492993715   2.3972401372   2.5654525414   2.7540097649   2.9629877484   3.1924470796   3.4424174121   3.7128890257   4.0038126419
//     -3.0066739103  -2.8157708146  -2.6222968256  -2.4256427719  -2.2258060430  -2.0233146753  -1.8190334957  -1.6139268416  -1.4088512490  -1.2044257489  -1.0009921403
//      0.0178093036   0.1160809591   0.2199508346   0.3356953265   0.4696003900   0.6279080030   0.8167815913   1.0422953680   1.3104447850   1.6271700977   1.9983835684

    DoubleVector x;
    x << 2.3972401372 << -2.4256427719 << 0.3356953265;
    printf("%14.10f\n", prob1.fx(x,2));

    IPrinter::printSeperatorLine();
    prob1.printResult();
}

Problem4Ex1::Problem4Ex1(const ODEGrid &grid) : mgrid(grid)
{
}

void Problem4Ex1::initialize()
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
    nsc1.time = 0.2;
    nsc1.nmbr = N/5;
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
    nsc3.time = 0.8;
    nsc3.nmbr = 4*(N/5);
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
    zett2.calculateForward(z2);
    zett2.calculateBackwardCP(z2,zm2);

//    DoubleVector x1(N+1);
//    DoubleVector x2(N+1);
//    DoubleVector x3(N+1);
//    for (unsigned int i=0; i<=N; i++)
//    {
//        x1[i] =  zm0[0][i] + (zm1[0][0][i]*g(0,0) + zm1[0][1][i]*g(0,1) + zm1[0][2][i]*g(0,2))
//                           + (zm2[0][0][i]*g(1,0) + zm2[0][1][i]*g(1,1) + zm2[0][2][i]*g(1,2));
//        x2[i] =  zm0[1][i] + (zm1[1][0][i]*g(0,0) + zm1[1][1][i]*g(0,1) + zm1[1][2][i]*g(0,2))
//                           + (zm2[1][0][i]*g(1,0) + zm2[1][1][i]*g(1,1) + zm2[1][2][i]*g(1,2));
//        x3[i] =  zm0[2][i] + (zm1[2][0][i]*g(0,0) + zm1[2][1][i]*g(0,1) + zm1[2][2][i]*g(0,2))
//                           + (zm2[2][0][i]*g(1,0) + zm2[2][1][i]*g(1,1) + zm2[2][2][i]*g(1,2));
//    }
//    IPrinter::printVector(x1);
//    IPrinter::printVector(x2);
//    IPrinter::printVector(x3);
}

void Problem4Ex1::printResult()
{
    Dimension dim = mgrid.dimension();
    unsigned int N = dim.sizeN();

    DoubleVector x1(N+1);
    DoubleVector x2(N+1);
    DoubleVector x3(N+1);
    for (unsigned int i=0; i<=N; i++)
    {
        x1[i] =  zm0[0][i] + (zm1[0][0][i]*g(0,0) + zm1[0][1][i]*g(0,1) + zm1[0][2][i]*g(0,2))
                + (zm2[0][0][i]*g(1,0) + zm2[0][1][i]*g(1,1) + zm2[0][2][i]*g(1,2));
        x2[i] =  zm0[1][i] + (zm1[1][0][i]*g(0,0) + zm1[1][1][i]*g(0,1) + zm1[1][2][i]*g(0,2))
                + (zm2[1][0][i]*g(1,0) + zm2[1][1][i]*g(1,1) + zm2[1][2][i]*g(1,2));
        x3[i] =  zm0[2][i] + (zm1[2][0][i]*g(0,0) + zm1[2][1][i]*g(0,1) + zm1[2][2][i]*g(0,2))
                + (zm2[2][0][i]*g(1,0) + zm2[2][1][i]*g(1,1) + zm2[2][2][i]*g(1,2));
    }
    IPrinter::printVector(10,6,x1);
    IPrinter::printVector(10,6,x2);
    IPrinter::printVector(10,6,x3);
}

double Problem4Ex1::fx(const DoubleVector &x, unsigned int num) const
{
    Dimension dim = mgrid.dimension();
    unsigned int N = dim.sizeN();

    if (num == 0)
    {
        unsigned int n = 0.3*N;
        return zm0[0][n] + (zm1[0][0][n]*g(x,0,0) + zm1[0][1][n]*g(x,0,1) + zm1[0][2][n]*g(x,0,2))
                         + (zm2[0][0][n]*g(x,1,0) + zm2[0][1][n]*g(x,1,1) + zm2[0][2][n]*g(x,1,2));
    }
    if (num == 1)
    {
        unsigned int n = 0.3*N;
        return zm0[1][n] + (zm1[1][0][n]*g(x,0,0) + zm1[1][1][n]*g(x,0,1) + zm1[1][2][n]*g(x,0,2))
                         + (zm2[1][0][n]*g(x,1,0) + zm2[1][1][n]*g(x,1,1) + zm2[1][2][n]*g(x,1,2));
    }
    if (num == 2)
    {
        unsigned int n = 0.3*N;
        return zm0[2][n] + (zm1[2][0][n]*g(x,0,0) + zm1[2][1][n]*g(x,0,1) + zm1[2][2][n]*g(x,0,2))
                         + (zm2[2][0][n]*g(x,1,0) + zm2[2][1][n]*g(x,1,1) + zm2[2][2][n]*g(x,1,2));
    }

    if (num == 3)
    {
        unsigned int n = 0.6*N;
        return zm0[0][n] + (zm1[0][0][n]*g(x,0,0) + zm1[0][1][n]*g(x,0,1) + zm1[0][2][n]*g(x,0,2))
                         + (zm2[0][0][n]*g(x,1,0) + zm2[0][1][n]*g(x,1,1) + zm2[0][2][n]*g(x,1,2));
    }
    if (num == 4)
    {
        unsigned int n = 0.6*N;
        return zm0[1][n] + (zm1[1][0][n]*g(x,0,0) + zm1[1][1][n]*g(x,0,1) + zm1[1][2][n]*g(x,0,2))
                         + (zm2[1][0][n]*g(x,1,0) + zm2[1][1][n]*g(x,1,1) + zm2[1][2][n]*g(x,1,2));
    }
    if (num == 5)
    {
        unsigned int n = 0.6*N;
        return zm0[2][n] + (zm1[2][0][n]*g(x,0,0) + zm1[2][1][n]*g(x,0,1) + zm1[2][2][n]*g(x,0,2))
                         + (zm2[2][0][n]*g(x,1,0) + zm2[2][1][n]*g(x,1,1) + zm2[2][2][n]*g(x,1,2));
    }
    return NAN;
}

double Problem4Ex1::A(double t, unsigned int, unsigned int row, unsigned int col) const
{
#ifdef SAMPLE_1
    if (row == 0) { if (col == 0) { return +2.0; } if (col == 1) { return t; }      if (col == 2) { return -3.0; } }
    if (row == 1) { if (col == 0) { return +3.0; } if (col == 1) { return -4.0*t; } if (col == 2) { return -8.0; } }
    if (row == 2) { if (col == 0) { return +t; }   if (col == 1) { return +1.0; }   if (col == 2) { return -1.0; } }
#endif
    return NAN;
}

double Problem4Ex1::B(double, unsigned int, unsigned int num, unsigned int row, unsigned int col) const
{
#ifdef SAMPLE_1

    if ( num == 0 )
    {
        if (row == 0) { if (col == 0) { return +0.002; } if (col == 1) { return +0.005; } if (col == 2) { return +0.003; } }
        if (row == 1) { if (col == 0) { return +0.004; } if (col == 1) { return +0.008; } if (col == 2) { return +0.001; } }
        if (row == 2) { if (col == 0) { return +0.001; } if (col == 1) { return +0.003; } if (col == 2) { return +0.004; } }
    }

    if ( num == 1 )
    {
        if (row == 0) { if (col == 0) { return +0.001; } if (col == 1) { return +0.003; } if (col == 2) { return +0.004; } }
        if (row == 1) { if (col == 0) { return +0.002; } if (col == 1) { return +0.005; } if (col == 2) { return +0.001; } }
        if (row == 2) { if (col == 0) { return +0.005; } if (col == 1) { return +0.002; } if (col == 2) { return +0.008; } }
    }

#endif
    return NAN;
}

double Problem4Ex1::C(double t, unsigned int k, unsigned int row) const
{
#ifdef SAMPLE_1
    if (row == 0)
    {
        return dX(t,0) - (A(t,k,0,0)*X(t,0) + A(t,k,0,1)*X(t,1) + A(t,k,0,2)*X(t,2))
                - (B(t,k,0,0,0)*g(0,0) + B(t,k,0,0,1)*g(0,1) + B(t,k,0,0,2)*g(0,2))
                - (B(t,k,1,0,0)*g(1,0) + B(t,k,1,0,1)*g(1,1) + B(t,k,1,0,2)*g(1,2));
    }
    if (row == 1)
    {
        return dX(t,1) - (A(t,k,1,0)*X(t,0) + A(t,k,1,1)*X(t,1) + A(t,k,1,2)*X(t,2))
                - (B(t,k,0,1,0)*g(0,0) + B(t,k,0,1,1)*g(0,1) + B(t,k,0,1,2)*g(0,2))
                - (B(t,k,1,1,0)*g(1,0) + B(t,k,1,1,1)*g(1,1) + B(t,k,1,1,2)*g(1,2));
    }
    if (row == 2)
    {
        return dX(t,2) - (A(t,k,2,0)*X(t,0) + A(t,k,2,1)*X(t,1) + A(t,k,2,2)*X(t,2))
                - (B(t,k,0,2,0)*g(0,0) + B(t,k,0,2,1)*g(0,1) + B(t,k,0,2,2)*g(0,2))
                - (B(t,k,1,2,0)*g(1,0) + B(t,k,1,2,1)*g(1,1) + B(t,k,1,2,2)*g(1,2));
    }
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
        if (row == 2) { return X(0.3, 2); }
    }
    if (num == 1)
    {
        if (row == 0) { return X(0.6,0)*X(0.6,0); }
        if (row == 1) { return X(0.6,1); }
        if (row == 2) { return X(0.6,2)*X(0.6,2); }
    }

//        if (num == 0)
//        {
//            if (row == 0) { return sin(X(0.3, 0)); }
//            if (row == 1) { return cos(X(0.3, 1)); }
//            if (row == 2) { return tan(X(0.3, 2)); }
//        }
//        if (num == 1)
//        {
//            if (row == 0) { return cos(X(0.6,0)); }
//            if (row == 1) { return sin(X(0.6,1)); }
//            if (row == 2) { return cos(X(0.6,2)); }
//        }

#endif
    return NAN;
}

double Problem4Ex1::g(const DoubleVector &x, unsigned int num, unsigned int row) const
{
    double x1_03 = x[0]; double x1_06 = x[3];
    double x2_03 = x[1]; double x2_06 = x[4];
    double x3_03 = x[2]; double x3_06 = x[5];


    if (num == 0)
    {
        if (row == 0) return sin(x1_03);
        if (row == 1) return cos(x2_03);
        if (row == 2) return tan(x3_03);
    }
    if (num == 1)
    {
        if (row == 0) return cos(x1_06);
        if (row == 1) return sin(x2_06);
        if (row == 2) return cos(x3_06);
    }
    return NAN;
}

double Problem4Ex1::X(double t, unsigned int i) const
{
#ifdef SAMPLE_1
    if (i == 0) { return t*t+t+2.0; }
    if (i == 1) { return 2.0*t-3.0; }
    if (i == 2) { return t*t*t+t; }
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

#include "cpnlcs.h"

void CauchyProblemNonLocalContions::Main(int agrc UNUSED_PARAM, char *argv[] UNUSED_PARAM)
{
    Dimension grid(0.01, 100);
    CauchyProblemNonLocalContions cpnlcs(grid);
    cpnlcs.calculate1();
}

CauchyProblemNonLocalContions::CauchyProblemNonLocalContions(const Dimension &grid) : CauchyProblemM(grid)
{
    unsigned int N = grid.sizeN();
    double h = grid.step();

    times << 0.0 << 1.0;

    L = times.size();
    alpha.resize(L);

    DoubleMatrix &m0 = alpha[0];
    m0.resize(n,n);
    m0[0][0] = 1.0; m0[0][1] = 3.0;
    m0[1][0] = 2.0; m0[1][1] = 1.0;

    //    DoubleMatrix &m1 = alpha[1];
    //    m1.resize(n,n);
    //    m1[0][0] = 0.0; m1[0][1] = 0.0;
    //    m1[1][0] = 0.0; m1[1][1] = 0.0;

    DoubleMatrix &m2 = alpha[1];
    m2.resize(n,n);
    m2[0][0] = 5.0; m2[0][1] = 4.0;
    m2[1][0] = 8.0; m2[1][1] = 1.0;

    betta.resize(n);
    betta[0] = m0[0][0]*x1(0) + m0[0][1]*x2(0) /*+ m1[0][0]*x1(N/2) + m1[0][1]*x2(N/2)*/ + m2[0][0]*x1(N) + m2[0][1]*x2(N);
    betta[1] = m0[1][0]*x1(0) + m0[1][1]*x2(0) /*+ m1[1][0]*x1(N/2) + m1[1][1]*x2(N/2)*/ + m2[1][0]*x1(N) + m2[1][1]*x2(N);
}

void CauchyProblemNonLocalContions::calculate1()
{
    unsigned int N = grid().sizeN();
    double h = grid().step();

    DoubleMatrix &alpha0 = alpha[0];
    DoubleMatrix &alpha1 = alpha[1];

    DoubleVector x0(4);
    x0[0] = alpha0[0][0];
    x0[1] = alpha0[0][1];
    x0[2] = betta[0];
    x0[3] = 1.0;
    double t0 = 0.0;

    DoubleMatrix rx;
    calculate(t0, x0, rx, RK4);

    DoubleVector alpha0_11 = rx.row(0);
    DoubleVector alpha0_12 = rx.row(1);
    DoubleVector betta__1  = rx.row(2);
    DoubleVector M         = rx.row(3);

    DoubleVector alpha1_11(N+1);
    DoubleVector alpha1_12(N+1);

    for (unsigned int i=0; i<=N; i++)
    {
        alpha1_11[i] = M[i]*alpha1[0][0];
        alpha1_12[i] = M[i]*alpha1[0][1];
    }

    IPrinter::printVector(14, 10, alpha0_11);
    IPrinter::printVector(14, 10, alpha0_12);
    IPrinter::printVector(14, 10, betta__1);
    IPrinter::printVector(14, 10, M);
    IPrinter::printVector(14, 10, alpha1_11);
    IPrinter::printVector(14, 10, alpha1_12);

    double a_1_11 = alpha0_11[N] + alpha1_11[N];
    double a_1_12 = alpha0_12[N] + alpha1_12[N];

    printf("%f %f\n", a_1_11*x1(N) + a_1_12*x2(N), betta__1[N]);
}

double CauchyProblemNonLocalContions::A(unsigned int k, unsigned int i, unsigned int j) const
{
    unsigned int N = grid().sizeN();
    double h = grid().step();

    double t = k*h;
    if (i==1)
    {
        if (j==1) { return t; }
        if (j==2) { return -1.0; }
    }
    if (i==2)
    {
        if (j==1) { return +2.0*t+3.0; }
        if (j==2) { return -2.0; }
    }
    return 0.0;
}

double CauchyProblemNonLocalContions::B(unsigned int k, unsigned int i) const
{
    double h = grid().step();

    double t = k*h;
    if (i==1)
    {
        return -t;
    }
    if (i==2)
    {
        return -6.0*t-11.0;
    }
    return 0.0;
}

double CauchyProblemNonLocalContions::x1(unsigned int k) const
{
    double h = grid().step();

    double t = k*h;
    return t*t + 4.0;
}

double CauchyProblemNonLocalContions::x2(unsigned int k) const
{
    double h = grid().step();

    double t = k*h;
    return t*t*t + t;
}

double CauchyProblemNonLocalContions::f(double t, const DoubleVector &x, unsigned int k, unsigned int i) const
{
    //Dimension gr = grid();
    double _S0 = S0(1,t,x,k);

    double alpha0_11 = x[0];
    double alpha0_12 = x[1];
    double betta___1 = x[2];
    double M         = x[3];

    if (i == 0) return _S0*alpha0_11 - ( A(k,1,1)*alpha0_11 + A(k,2,1)*alpha0_12 );
    if (i == 1) return _S0*alpha0_12 - ( A(k,1,2)*alpha0_11 + A(k,2,2)*alpha0_12 );
    if (i == 2) return _S0*betta___1 + ( B(k,1)*alpha0_11   + B(k,2)*alpha0_12 );
    if (i == 3) return _S0*M;

    return 0.0;
}

double CauchyProblemNonLocalContions::S0(unsigned int i, double t, const DoubleVector &x, unsigned int k) const
{
    double alpha0_11 = x[0];
    double alpha0_12 = x[1];
    double betta___1 = x[2];
    double M         = x[3];

    double s1,s2,m1;
    if (i==1)
    {
        s1 = (alpha0_11*A(k,1,1) + alpha0_12*A(k,2,1))*alpha0_11 + (alpha0_11*A(k,1,2) + alpha0_12*A(k,2,2))*alpha0_12;
        s2 = (alpha0_11*B(k,1) + alpha0_12*B(k,2))*betta___1;
        m1 = alpha0_11*alpha0_11 + alpha0_12*alpha0_12 + betta___1*betta___1;
    }
    if (i==2)
    {
        s1 = (alpha0_11*A(k,1,1) + alpha0_12*A(k,2,1))*alpha0_11 + (alpha0_11*A(k,1,2) + alpha0_12*A(k,2,2))*alpha0_12;
        s2 = (alpha0_11*B(k,1) + alpha0_12*B(k,2))*betta___1;
        m1 = alpha0_11*alpha0_11 + alpha0_12*alpha0_12 + betta___1*betta___1;
    }

    return (s1-s2)/m1;
}

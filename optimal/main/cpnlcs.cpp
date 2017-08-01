#include "cpnlcs.h"

void CauchyProblemNonLocalContions::Main(int agrc UNUSED_PARAM, char *argv[] UNUSED_PARAM)
{
    Dimension grid(0.01, 100);
    CauchyProblemNonLocalContions cpnlcs(grid);
    //cpnlcs.calculate1();
    cpnlcs.initialize();

    cpnlcs.calculateInterval(0,1,0);
    cpnlcs.calculateInterval(1,2,0);

    puts("---");
    cpnlcs.calculateInterval(0,1,1);
    cpnlcs.calculateInterval(1,2,1);
}

CauchyProblemNonLocalContions::CauchyProblemNonLocalContions(const Dimension &grid) : CauchyProblemM1stOrder(grid)
{
}

void CauchyProblemNonLocalContions::initialize()
{
    unsigned int N = grid().sizeN();

    n = 2;

    NonSeparatedCondition cond0;
    cond0.time = 0.0;
    cond0.rows = cond0.cols = n;
    cond0.alpha.resize(cond0.rows,cond0.cols);
    cond0.alpha[0][0] = 1.0; cond0.alpha[0][1] = 3.0;
    cond0.alpha[1][0] = 2.0; cond0.alpha[1][1] = 1.0;
    cond0.n = 0;

    NonSeparatedCondition cond1;
    cond1.time = 0.5;
    cond1.rows = cond1.cols = n;
    cond1.alpha.resize(cond1.rows,cond1.cols);
    cond1.alpha[0][0] = 0.0; cond1.alpha[0][1] = 0.0;
    cond1.alpha[1][0] = 0.0; cond1.alpha[1][1] = 0.0;
    cond1.n = N/2;

    NonSeparatedCondition cond2;
    cond2.time = 1.0;
    cond2.rows = cond2.cols = n;
    cond2.alpha.resize(cond2.rows,cond2.cols);
    cond2.alpha[0][0] = 5.0; cond2.alpha[0][1] = 4.0;
    cond2.alpha[1][0] = 8.0; cond2.alpha[1][1] = 1.0;
    cond2.n = N;

    nscs.resize(3);
    nscs[0] = cond0;
    nscs[1] = cond1;
    nscs[2] = cond2;

    L = nscs.size();

    betta.resize(n);
    betta[0] = cond0.alpha[0][0]*x1(0)   + cond0.alpha[0][1]*x2(0) +
               cond1.alpha[0][0]*x1(N/2) + cond1.alpha[0][1]*x2(N/2) +
               cond2.alpha[0][0]*x1(N)   + cond2.alpha[0][1]*x2(N);

    betta[1] = cond0.alpha[1][0]*x1(0)   + cond0.alpha[1][1]*x2(0) +
               cond1.alpha[1][0]*x1(N/2) + cond1.alpha[1][1]*x2(N/2) +
               cond2.alpha[1][0]*x1(N)   + cond2.alpha[1][1]*x2(N);

    M.resize(N+1);
    M[0] = 1.0;
}

void CauchyProblemNonLocalContions::calculateInterval(unsigned int start, unsigned int end, unsigned int r)
{
    //puts("--------------------");
    unsigned int N = grid().sizeN();

    unsigned int L = nscs.size();

    NonSeparatedCondition &sc = nscs.at(start);
    NonSeparatedCondition &ec = nscs.at(end);

    mgrid = Dimension(grid().step(), ec.n, sc.n);
    //printf("%d %d\n", mgrid.minN(), mgrid.maxN());

    DoubleMatrix a0 = sc.alpha;

    DoubleVector x;
    x << a0[r][0] << a0[r][1] << betta[r] << 1.0;//M[sc.n];
    //printf("%f %f %f %f %d %d\n", x[0], x[1], x[2], x[3], sc.n, ec.n);

    DoubleVector rx(x.size());
    calculateCP(sc.time, x, rx, RK4);

    M[ec.n] = rx[3];
    //printf("%f\n", M[ec.n]);

    for (unsigned int s=end; s<L; s++)
    {
        NonSeparatedCondition &cc = nscs.at(s);
        //for (unsigned int i=0; i<n; i++) printf("%d %d %f %f\n", s, i, ec.alpha[r][i], cc.alpha[r][i]);
        for (unsigned int i=0; i<n; i++) cc.alpha[r][i] *= M[ec.n];
        //for (unsigned int i=0; i<n; i++) printf("%d %d %f %f\n", s, i, ec.alpha[r][i], cc.alpha[r][i]);
    }
    //for (unsigned int i=0; i<n; i++) printf("%f\n", ec.alpha[r][i]);
    for (unsigned int i=0; i<n; i++) ec.alpha[r][i] += rx[i];
    //for (unsigned int i=0; i<n; i++) printf("%f\n", ec.alpha[r][i]);
    betta[r] = rx[2];

    if (start==0)
    {
       printf("%f %f %d %d\n", (nscs[1].alpha[r][0]*x1(50)  + nscs[1].alpha[r][1]*x2(50)) +
                               (nscs[2].alpha[r][0]*x1(100) + nscs[2].alpha[r][1]*x2(100)), betta[r], ec.n, N);
    }
    if (start==1)
    {
        printf("%f %f %d %d\n", nscs[2].alpha[r][0]*x1(100) + nscs[2].alpha[r][1]*x2(100), betta[r], ec.n, N);
        //printf("%f %f %d %d\n", ec.alpha[r][0]*x1(N) + ec.alpha[r][1]*x2(N), betta[r], ec.n, N);
    }
}

void CauchyProblemNonLocalContions::calculate1()
{
    unsigned int N = grid().sizeN();
    double h = grid().step();

    DoubleMatrix &alpha0 = nscs.at(0).alpha;
    DoubleMatrix &alpha1 = nscs.at(1).alpha;

    DoubleVector x0(4);
    x0[0] = alpha0[0][0];
    x0[1] = alpha0[0][1];
    x0[2] = betta[0];
    x0[3] = 1.0;
    double t0 = 0.0;

    DoubleMatrix rx;
    calculateCP(t0, x0, rx, RK4);

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

double CauchyProblemNonLocalContions::A(double, unsigned int k, unsigned int row, unsigned int col) const
{
    unsigned int N = grid().sizeN();
    double h = grid().step();
    double t = k*h;

    if (row==1)
    {
        if (col==1) { return t; }
        if (col==2) { return -1.0; }
    }
    if (row==2)
    {
        if (col==1) { return +2.0*t+3.0; }
        if (col==2) { return -2.0; }
    }
    return 0.0;
}

double CauchyProblemNonLocalContions::B(double, unsigned int k, unsigned int row) const
{
    double h = grid().step();
    double t = k*h;

    if (row==1) return -t;
    if (row==2) return -6.0*t-11.0;
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

    if (i == 0) return _S0*alpha0_11 - ( A(t,k,1,1)*alpha0_11 + A(t,k,2,1)*alpha0_12 );
    if (i == 1) return _S0*alpha0_12 - ( A(t,k,1,2)*alpha0_11 + A(t,k,2,2)*alpha0_12 );
    if (i == 2) return _S0*betta___1 + ( B(t,k,1)*alpha0_11   + B(t,k,2)*alpha0_12 );
    if (i == 3) return _S0*M;

    return 0.0;
}

double CauchyProblemNonLocalContions::S0(unsigned int i, double t, const DoubleVector &x, unsigned int k) const
{
    double alpha0_11 = x[0];
    double alpha0_12 = x[1];
    double betta___1 = x[2];
    double M         = x[3];

    double s1=0.0,s2=0.0,m1=0.0;
    if (i==1)
    {
        s1 = (alpha0_11*A(t,k,1,1) + alpha0_12*A(t,k,2,1))*alpha0_11 + (alpha0_11*A(t,k,1,2) + alpha0_12*A(t,k,2,2))*alpha0_12;
        s2 = (alpha0_11*B(t,k,1) + alpha0_12*B(t,k,2))*betta___1;
        m1 = alpha0_11*alpha0_11 + alpha0_12*alpha0_12 + betta___1*betta___1;
    }
    if (i==2)
    {
        s1 = (alpha0_11*A(t,k,1,1) + alpha0_12*A(t,k,2,1))*alpha0_11 + (alpha0_11*A(t,k,1,2) + alpha0_12*A(t,k,2,2))*alpha0_12;
        s2 = (alpha0_11*B(t,k,1) + alpha0_12*B(t,k,2))*betta___1;
        m1 = alpha0_11*alpha0_11 + alpha0_12*alpha0_12 + betta___1*betta___1;
    }

    return (s1-s2)/m1;
}

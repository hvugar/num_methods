#include "cpnlcs.h"

void CauchyProblemNonLocalContions::Main(int agrc UNUSED_PARAM, char *argv[] UNUSED_PARAM)
{
    Dimension grid(0.0001, 10000);
    CauchyProblemNonLocalContions cpnlcs(grid);

    cpnlcs.initialize();
    cpnlcs.calculate1();
}

CauchyProblemNonLocalContions::CauchyProblemNonLocalContions(const ODEGrid &grid) : SystemLinearODE1stOrder(grid) {}

void CauchyProblemNonLocalContions::initialize()
{
    Dimension dim = grid().dimension();
    unsigned int N = dim.sizeN();

    n = 2;

    NonSeparatedCondition cond0;
    cond0.time = 0.0;
    cond0.nmbr = 0;
    cond0.rows = cond0.cols = n;
    cond0.alpha.resize(cond0.rows,cond0.cols);
    cond0.alpha[0][0] = 1.0; cond0.alpha[0][1] = 3.0;
    cond0.alpha[1][0] = 2.0; cond0.alpha[1][1] = 1.0;
    for (unsigned int i=0; i<n; i++) for(unsigned int j=0; j<n; j++) cond0.alpha[i][j] = (rand() % 100);// / 1000.0;

    NonSeparatedCondition cond1;
    cond1.time = 0.5;
    cond1.nmbr = N/2;
    cond1.rows = cond1.cols = n;
    cond1.alpha.resize(cond1.rows,cond1.cols);
    cond1.alpha[0][0] = 0.0; cond1.alpha[0][1] = 0.0;
    cond1.alpha[1][0] = 0.0; cond1.alpha[1][1] = 0.0;
    for (unsigned int i=0; i<n; i++) for (unsigned int j=0; j<n; j++) cond1.alpha[i][j] = (rand() % 100);// / 1000.0;

    NonSeparatedCondition cond2;
    cond2.time = 1.0;
    cond2.nmbr = N;
    cond2.rows = cond2.cols = n;
    cond2.alpha.resize(cond2.rows,cond2.cols);
    cond2.alpha[0][0] = 5.0; cond2.alpha[0][1] = 4.0;
    cond2.alpha[1][0] = 8.0; cond2.alpha[1][1] = 1.0;
    for (unsigned int i=0; i<n; i++) for (unsigned int j=0; j<n; j++) cond2.alpha[i][j] = (rand() % 100);// / 1000.0;

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
}

void CauchyProblemNonLocalContions::calculateInterval(unsigned int start, unsigned int r)
{
    //unsigned int N = grid().sizeN();

    unsigned int L = nscs.size();

    NonSeparatedCondition &sc = nscs.at(start);
    NonSeparatedCondition &ec = nscs.at(start+1);

    //mgrid = Dimension(grid().step(), ec.nmbr, sc.nmbr);
    //printf("%d %d\n", mgrid.minN(), mgrid.maxN());

    DoubleVector x;
    x << sc.alpha[r][0] << sc.alpha[r][1] << betta[r] << 1.0;
    //printf("%d %f %f\n", 0, sc.alpha[r][0],sc.alpha[r][1]);
    //printf("%f %f %f %f %d %d\n", x[0], x[1], x[2], x[3], sc.nmbr, ec.nmbr);

    DoubleVector rx(x.size());

    Dimension dim(grid().dimension().step(), ec.nmbr, sc.nmbr);
    CauchyProblemM1stOrderA cpa(*this, dim);
    cpa.calculateCP(sc.time, x, rx, InitialValueProblem::RK4);

    nscs.at(start).alpha[r][0] = rx[0];
    nscs.at(start).alpha[r][1] = rx[1];
    betta[r] = rx[2];
    double M = rx[3];

    //printf("%f %f %f %f %d %d\n", rx[0], rx[1], rx[2], rx[3], sc.nmbr, ec.nmbr);
    //printf("%d %f %f\n", 0, rx[0], rx[1]);

    for (unsigned int s=start+1; s<L; s++)
    {
        NonSeparatedCondition &cc = nscs.at(s);
        //printf("%d %d %f %f %f\n", r, s, cc.alpha[r][0], cc.alpha[r][1], M);
        for (unsigned int i=0; i<n; i++) cc.alpha[r][i] *= M;
        //printf("%d %d %f %f %f\n", r, s, cc.alpha[r][0], cc.alpha[r][1], M);
    }

    for (unsigned int i=0; i<n; i++) ec.alpha[r][i] += rx[i];
    //printf("%d %d %f %f %f\n", r, start+1, nscs.at(start+1).alpha[r][0], nscs.at(start+1).alpha[r][1], M);


    if (start==0)
    {
        printf("%f %f %d\n", (nscs[1].alpha[r][0]*x1(50)  + nscs[1].alpha[r][1]*x2(50)) +
                (nscs[2].alpha[r][0]*x1(100) + nscs[2].alpha[r][1]*x2(100)), betta[r], ec.nmbr);
    }
    if (start==1)
    {
        printf("%f %f %d\n", nscs[2].alpha[r][0]*x1(100) + nscs[2].alpha[r][1]*x2(100), betta[r], ec.nmbr);
    }
}

void CauchyProblemNonLocalContions::calculateCondition(unsigned int r UNUSED_PARAM)
{
}

void CauchyProblemNonLocalContions::calculate1()
{
    for (unsigned int s=0; s<L-1; s++)
    {
        for (unsigned int row=0; row<n; row++)
        {
            calculateInterval(s,row);
        }
        puts("---");
    }

    NonSeparatedCondition cond1 = nscs.at(nscs.size()-1);

    DoubleVector x(2);
    GaussianElimination(cond1.alpha, betta, x);
    printf("%f %f\n", x[0], x[1]);

    CauchyProblemM1stOrderB cpb(*this, grid());
    DoubleMatrix m;
    cpb.calculateCP(1.0, x, m, CauchyProblemM1stOrder::RK4, CauchyProblemM1stOrder::R2L);
    IPrinter::printVector(m.row(0));
    IPrinter::printVector(m.row(1));
}

double CauchyProblemNonLocalContions::A(double, unsigned int k, unsigned int row, unsigned int col) const
{
    Dimension dim = grid().dimension();
    double h = dim.step();
    double t = k*h;

#ifdef SAMPLE_1
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

    return NAN;
}

double CauchyProblemNonLocalContions::B(double, unsigned int k, unsigned int row) const
{
    Dimension dim = grid().dimension();
    double h = dim.step();
    double t = k*h;
#ifdef SAMPLE_1
    if (row==1) return -t;
    if (row==2) return -6.0*t-11.0;
#endif
#ifdef SAMPLE_2
    if (row==1) return +2.0 - 13.0*t;
    if (row==2) return +3.0 - 22.0*t;
#endif
    return NAN;
}

double CauchyProblemNonLocalContions::x1(unsigned int k) const
{
    Dimension dim = grid().dimension();
    double h = dim.step();
    double t = k*h;
#ifdef SAMPLE_1
    return t*t + 4.0;
#endif
#ifdef SAMPLE_2
    return 2.0*t;
#endif
}

double CauchyProblemNonLocalContions::x2(unsigned int k) const
{
    Dimension dim = grid().dimension();
    double h = dim.step();
    double t = k*h;
#ifdef SAMPLE_1
    return t*t*t + t;
#endif
#ifdef SAMPLE_2
    return 3.0*t;
#endif
}

double CauchyProblemM1stOrderA::f(double t, const DoubleVector &x, unsigned int k, unsigned int i) const
{
    double _SO = S0(t,x,k);

    double ar1 = x[0];
    double ar2 = x[1];
    double btr = x[2];
    double M   = x[3];

    if (i == 0) return _SO*ar1 - ( p.A(t,k,1,1)*ar1 + p.A(t,k,2,1)*ar2 );
    if (i == 1) return _SO*ar2 - ( p.A(t,k,1,2)*ar1 + p.A(t,k,2,2)*ar2 );
    if (i == 2) return _SO*btr + ( p.B(t,k,1)*ar1   + p.B(t,k,2)*ar2 );
    if (i == 3) return _SO*M;

    return NAN;
}

double CauchyProblemM1stOrderA::S0(double t, const DoubleVector &x, unsigned int k) const
{
    double ar1 = x[0];
    double ar2 = x[1];
    double btr = x[2];
    double M   = x[3];

    double s1 = ( ar1*p.A(t,k,1,1) + ar2*p.A(t,k,2,1) )*ar1 + ( ar1*p.A(t,k,1,2) + ar2*p.A(t,k,2,2) )*ar2;
    double s2 = ( ar1*p.B(t,k,1)   + ar2*p.B(t,k,2) )*btr;
    double m1 = ar1*ar1 + ar2*ar2 + btr*btr;

    return (s1-s2)/m1;
}

double CauchyProblemM1stOrderB::f(double t, const DoubleVector &x, unsigned int k, unsigned int i) const
{
    double x1 = x[0];
    double x2 = x[1];

    if (i == 0) return p.A(t,k,1,1)*x1 + p.A(t,k,1,2)*x2 + p.B(t,k,1);
    if (i == 1) return p.A(t,k,2,1)*x1 + p.A(t,k,2,2)*x2 + p.B(t,k,2);

    return NAN;
}

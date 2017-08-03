#include "cpnlcs.h"

void CauchyProblemNonLocalContions::Main(int agrc UNUSED_PARAM, char *argv[] UNUSED_PARAM)
{
    ODEGrid grid(Dimension(0.001, 1000, 0));
    CauchyProblemNonLocalContions cpnlcs(grid);

    cpnlcs.initialize();
    DoubleVector x;
    cpnlcs.calculateForward(x);
    IPrinter::print(x,x.size());

    CauchyProblemM1stOrderB cpb(cpnlcs, grid);
    DoubleMatrix m;
    cpb.calculateCP(1.0, x, m, CauchyProblemM1stOrder::RK4, CauchyProblemM1stOrder::R2L);
    for (unsigned int i=0; i<cpnlcs.n; i++) IPrinter::printVector(m.row(i));
}


CauchyProblemNonLocalContions::CauchyProblemNonLocalContions(const ODEGrid &grid) : SystemLinearODE1stOrder(grid) {}

void CauchyProblemNonLocalContions::initialize()
{
    Dimension dim = grid().dimension();
    unsigned int N = dim.sizeN();

    n0 = 0;
    n1 = 3;
    n2 = 0;
    n = n0 + n1 + n2;

    Condition nsc0;
    nsc0.type = NonSeparated;
    nsc0.time = 0.0;
    nsc0.nmbr = 0;
    nsc0.alpha.resize(n0, n);
    for (unsigned int row=0; row<n0; row++) for(unsigned int col=0; col<n; col++) nsc0.alpha[row][col] = (rand() % 1000) / 1000.0;

    Condition nsc1;
    nsc1.type = NonSeparated;
    nsc1.time = 0.25;
    nsc1.nmbr = N/4;
    nsc1.alpha.resize(n0, n);
    for (unsigned int row=0; row<n0; row++) for(unsigned int col=0; col<n; col++) nsc1.alpha[row][col] = (rand() % 1000) / 1000.0;

    Condition nsc2;
    nsc2.type = NonSeparated;
    nsc2.time = 0.5;
    nsc2.nmbr = N/2;
    nsc2.alpha.resize(n0, n);
    for (unsigned int row=0; row<n0; row++) for(unsigned int col=0; col<n; col++) nsc2.alpha[row][col] = (rand() % 1000) / 1000.0;

    Condition nsc3;
    nsc3.type = NonSeparated;
    nsc3.time = 0.75;
    nsc3.nmbr = 3*(N/4);
    nsc3.alpha.resize(n0, n);
    for (unsigned int row=0; row<n0; row++) for(unsigned int col=0; col<n; col++) nsc3.alpha[row][col] = (rand() % 1000) / 1000.0;

    Condition nsc4;
    nsc4.type = NonSeparated;
    nsc4.time = 1.0;
    nsc4.nmbr = N;
    nsc4.alpha.resize(n0, n);
    for (unsigned int row=0; row<n0; row++) for(unsigned int col=0; col<n; col++) nsc4.alpha[row][col] = (rand() % 1000) / 1000.0;

    nscs.resize(5);
    nscs[0] = nsc0;
    nscs[1] = nsc1;
    nscs[2] = nsc2;
    nscs[3] = nsc3;
    nscs[4] = nsc4;

    L = nscs.size();

    betta.resize(n);

    for (unsigned int row=0; row<n0; row++)
    {
        betta[row] = 0.0;
        for (unsigned int s=0; s<L; s++)
        {
            Condition &c = nscs.at(s);
            for (unsigned int i=0; i<n; i++) betta[row] += c.alpha[row][i] * x(c.nmbr, i);
        }
    }

    lscs.type = SeparatedLeft;
    lscs.time = 0.0;
    lscs.nmbr = 0;
    lscs.alpha.resize(n1, n);
    for (unsigned int row=0; row<n1; row++)
    {
        for(unsigned int col=0; col<n; col++) lscs.alpha[row][col] = (rand() % 1000) / 1000.0;
        betta[row+n0] = 0.0;
        for (unsigned int i=0; i<n; i++) betta[row+n0] += lscs.alpha[row][i] * x(lscs.nmbr, i);
    }

    rscs.type = SeparatedRight;
    rscs.time = 1.0;
    rscs.nmbr = N;
    rscs.alpha.resize(n2, n);
    for (unsigned int row=0; row<n2; row++)
    {
        for(unsigned int col=0; col<n; col++) rscs.alpha[row][col] = (rand() % 1000) / 1000.0;
        betta[row+n0+n1] = 0.0;
        for (unsigned int i=0; i<n; i++) betta[row+n0+n1] += rscs.alpha[row][i] * x(rscs.nmbr, i);
    }
}

void CauchyProblemNonLocalContions::calculateIntervalF(unsigned int start, unsigned int r)
{
    unsigned int L = nscs.size();
    double h = grid().dimension().step();
    DoubleVector x(n+2);
    DoubleVector rx(n+2);

    Condition &sc = nscs.at(start);
    Condition &ec = nscs.at(start+1);

    for (unsigned int i=0; i<n; i++) x[i] = sc.alpha[r][i]; x[n] = betta[r]; x[n+1] = 1.0;

    Dimension dim(h, ec.nmbr, sc.nmbr);
    CauchyProblemM1stOrderA cpa(*this, ODEGrid(dim));
    cpa.calculateCP(sc.time, x, rx, InitialValueProblem::RK4);

    for (unsigned int i=0; i<n; i++) sc.alpha[r][i] = rx[i];
    betta[r] = rx[n];
    double M = rx[n+1];

    for (unsigned int s=start+1; s<L; s++)
    {
        Condition &cc = nscs.at(s);
        for (unsigned int i=0; i<n; i++) cc.alpha[r][i] *= M;
    }

    for (unsigned int i=0; i<n; i++) ec.alpha[r][i] += rx[i];

    x.clear();
    rx.clear();
}

void CauchyProblemNonLocalContions::calculateForward(DoubleVector &x)
{
    for (unsigned int row=0; row<n0; row++)
    {
        for (unsigned int s=0; s<L-1; s++)
        {
            calculateIntervalF(s,row);
        }
    }

    for (unsigned int row=0; row<n1; row++)
    {
        double h = grid().dimension().step();
        unsigned int minN = grid().dimension().minN();
        unsigned int maxN = grid().dimension().maxN();

        DoubleVector x(n+2);
        DoubleVector rx(n+2);

        for (unsigned int i=0; i<n; i++) x[i] = lscs.alpha[row][i]; x[n] = betta[row+n0]; x[n+1] = 1.0;

        Dimension dim(h, maxN, minN);
        CauchyProblemM1stOrderA cpa(*this, ODEGrid(dim));
        cpa.calculateCP(lscs.time, x, rx, InitialValueProblem::RK4);

        for (unsigned int i=0; i<n; i++) lscs.alpha[row][i] = rx[i];
        betta[row+n0] = rx[n];
        x.clear();
        rx.clear();
    }

    DoubleMatrix A(n, n);
    DoubleVector b(n);
    x.clear();
    x.resize(n);

    Condition c0 = nscs.at(L-1);
    for (unsigned int row=0; row<n0; row++)
    {
        for (unsigned int col=0; col<n; col++)
        {
            A[row][col] = c0.alpha[row][col];
        }
        b[row] = betta[row];
    }

    for (unsigned int row=0; row<n1; row++)
    {
        for (unsigned int col=0; col<n; col++)
        {
            A[row+n0][col] = lscs.alpha[row][col];
        }
        b[row+n0] = betta[row+n0];
    }

    for (unsigned int row=0; row<n2; row++)
    {
        for (unsigned int col=0; col<n; col++)
        {
            A[row+n0+n1][col] = rscs.alpha[row][col];
        }
        b[row+n0+n1] = betta[row+n0+n1];
    }

    GaussianElimination(A, b, x);
}

double CauchyProblemNonLocalContions::A(double t, unsigned int, unsigned int row, unsigned int col) const
{
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
#ifdef SAMPLE_3
    if (row==1)
    {
        if (col==1) { return +2.0; }
        if (col==2) { return -3.0; }
        if (col==3) { return +1.0; }
    }
    if (row==2)
    {
        if (col==1) { return +3.0; }
        if (col==2) { return +1.0; }
        if (col==3) { return -2.0; }
    }
    if (row==3)
    {
        if (col==1) { return +1.0; }
        if (col==2) { return -5.0; }
        if (col==3) { return -3.0; }
    }
#endif

    return NAN;
}

double CauchyProblemNonLocalContions::B(double t, unsigned int, unsigned int row) const
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
    if (row==1) return 3.0     - (2.0*(3.0*t+4.0) - 3.0*(4.0*t*t) + 1.0*(t*t+t));
    if (row==2) return 8.0*t   - (3.0*(3.0*t+4.0) + 1.0*(4.0*t*t) - 2.0*(t*t+t));
    if (row==3) return 2.0*t+1 - (1.0*(3.0*t+4.0) - 5.0*(4.0*t*t) - 3.0*(t*t+t));
#endif
    return NAN;
}

double CauchyProblemNonLocalContions::x(unsigned int k, int i) const
{
    Dimension dim = grid().dimension();
    double h = dim.step();
    double t = k*h;
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

double CauchyProblemM1stOrderA::f(double t, const DoubleVector &x, unsigned int k, unsigned int i) const
{
    double _SO = S0(t,x,k);

    unsigned int n = p.n;

    //double ar1 = x[0];
    //double ar2 = x[1];

    //double btr = x[2];
    //double M   = x[3];

    if (i<n)
    {
        double res = _SO*x[i];
        for (unsigned int j=0; j<n; j++) res -= p.A(t,k,j+1,i+1)*x[j];
        return res;
    }
    else if (i==n)
    {
        double res = _SO*x[n];
        for (unsigned int j=0; j<n; j++) res += p.B(t,k,j+1)*x[j];
        return res;
    }
    else
    {
        return _SO*x[n+1];
    }


//    if (i == 0) return _SO*ar1 - ( p.A(t,k,1,1)*ar1 + p.A(t,k,2,1)*ar2 );
//    if (i == 1) return _SO*ar2 - ( p.A(t,k,1,2)*ar1 + p.A(t,k,2,2)*ar2 );

//    if (i == 2) return _SO*btr + ( p.B(t,k,1)*ar1   + p.B(t,k,2)*ar2 );
//    if (i == 3) return _SO*M;

    return NAN;
}

double CauchyProblemM1stOrderA::S0(double t, const DoubleVector &x, unsigned int k) const
{
    unsigned int n = p.n;

    //double ar1 = x[0];
    //double ar2 = x[1];
    double btr = x[n];
    //double M   = x[n+1];

    //double s1 = ( ar1*p.A(t,k,1,1) + ar2*p.A(t,k,2,1) )*ar1 + ( ar1*p.A(t,k,1,2) + ar2*p.A(t,k,2,2) )*ar2;
    //double s2 = ( ar1*p.B(t,k,1)   + ar2*p.B(t,k,2) )*btr;
    //double m1 = ar1*ar1 + ar2*ar2 + btr*btr;

    double s1 = 0.0;
    for (unsigned int i=1; i<=n; i++)
    {
        double aa = 0.0;
        for (unsigned int j=1; j<=n; j++) aa += x[j-1]*p.A(t,k,j,i);
        s1 += aa*x[i-1];
    }

    double s2 = 0.0;
    for (unsigned int i=1; i<=n; i++)
    {
        s2 += x[i-1]*p.B(t,k,i);
    }
    s2 *= btr;


    double m1 = 0.0;
    for (unsigned int i=1; i<=n; i++)
    {
        m1 += x[i-1]*x[i-1];
    }
    m1 += btr*btr;

    return (s1-s2)/m1;
}

double CauchyProblemM1stOrderB::f(double t, const DoubleVector &x, unsigned int k, unsigned int i) const
{
//    double x1 = x[0];
//    double x2 = x[1];
//    double x3 = x[2];

    unsigned int n = p.n;
    double aa = p.B(t,k,i+1);
    for (unsigned int j=1; j<=n; j++) aa += p.A(t,k,i+1,j)*x[j-1];

//    if (i == 0) return p.A(t,k,1,1)*x1 + p.A(t,k,1,2)*x2 + p.B(t,k,1);
//    if (i == 1) return p.A(t,k,2,1)*x1 + p.A(t,k,2,2)*x2 + p.B(t,k,2);

//    if (i == 0) return p.A(t,k,1,1)*x1 + p.A(t,k,1,2)*x2 + p.A(t,k,1,3)*x3 + p.B(t,k,1);
//    if (i == 1) return p.A(t,k,2,1)*x1 + p.A(t,k,2,2)*x2 + p.A(t,k,2,3)*x3 + p.B(t,k,2);
//    if (i == 2) return p.A(t,k,3,1)*x1 + p.A(t,k,3,2)*x2 + p.A(t,k,3,3)*x3 + p.B(t,k,3);

    return aa;
}

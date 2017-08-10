#include "slodenlcsm.h"

class CauchyProblemM1stOrderAM : public CauchyProblemM1stOrder
{
public:
    CauchyProblemM1stOrderAM(ISystemLinearODENonLocalContionsM &parent, const ODEGrid& grid) : CauchyProblemM1stOrder(grid), p(parent) {}

    unsigned int row;
    unsigned int col;

protected:
    virtual double f(double t, const DoubleVector &x, unsigned int k, unsigned int i) const;
    double S0(double t, const DoubleVector &x, unsigned int k) const;
private:
    ISystemLinearODENonLocalContionsM &p;
};

void SystemLinearODENonLocalContionsM::Main(int agrc UNUSED_PARAM, char *argv[] UNUSED_PARAM)
{
    ODEGrid grid(Dimension(0.01, 100, 0));
    SystemLinearODENonLocalContionsM cpnlcs(grid);

    cpnlcs.initialize();
    DoubleMatrix x;
    cpnlcs.calculateForward(x);
//    IPrinter::print(x,x.size());
//    DoubleMatrix m;
//    cpnlcs.calculateBackwardCP(x, m);
//    for (unsigned int row=0; row<cpnlcs.systemOrder(); row++) IPrinter::printVector(m.row(row));
}

SystemLinearODENonLocalContionsM::SystemLinearODENonLocalContionsM(const ODEGrid &grid)
    : ISystemLinearODENonLocalContionsM(grid) {}

void SystemLinearODENonLocalContionsM::initialize()
{
//    bool a,b;
//    if (a and b) {}

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
    nsc0.alpha.resize(n0, n);
    for (unsigned int row=0; row<n0; row++) for(unsigned int col=0; col<n; col++) nsc0.alpha[row][col] = (rand() % 1000) / 1000.0;
    addNonSeparatedCondition(nsc0);

    IPrinter::print(nsc0.alpha, nsc0.alpha.rows(), nsc0.alpha.cols());

//    Condition nsc1;
//    nsc1.type = NonSeparated;
//    nsc1.time = 0.25;
//    nsc1.nmbr = N/4;
//    nsc1.alpha.resize(n0, n);
//    for (unsigned int row=0; row<n0; row++) for(unsigned int col=0; col<n; col++) nsc1.alpha[row][col] = (rand() % 1000) / 1000.0;
//    addNonSeparatedCondition(nsc1);

//    Condition nsc2;
//    nsc2.type = NonSeparated;
//    nsc2.time = 0.5;
//    nsc2.nmbr = N/2;
//    nsc2.alpha.resize(n0, n);
//    for (unsigned int row=0; row<n0; row++) for(unsigned int col=0; col<n; col++) nsc2.alpha[row][col] = (rand() % 1000) / 1000.0;
//    addNonSeparatedCondition(nsc2);

//    Condition nsc3;
//    nsc3.type = NonSeparated;
//    nsc3.time = 0.75;
//    nsc3.nmbr = 3*(N/4);
//    nsc3.alpha.resize(n0, n);
//    for (unsigned int row=0; row<n0; row++) for(unsigned int col=0; col<n; col++) nsc3.alpha[row][col] = (rand() % 1000) / 1000.0;
//    addNonSeparatedCondition(nsc3);

    Condition nsc4;
    nsc4.type = NonSeparated;
    nsc4.time = 1.0;
    nsc4.nmbr = N;
    nsc4.alpha.resize(n0, n);
    for (unsigned int row=0; row<n0; row++) for(unsigned int col=0; col<n; col++) nsc4.alpha[row][col] = (rand() % 1000) / 1000.0;
    addNonSeparatedCondition(nsc4);

    IPrinter::printSeperatorLine();
    IPrinter::print(nsc4.alpha, nsc4.alpha.rows(), nsc4.alpha.cols());
    IPrinter::printSeperatorLine();

    DoubleMatrix betta(n0, n);
    unsigned int L = nonSeparatedConditions().size();

    for (unsigned int row=0; row<n0; row++)
    {
        for (unsigned int col=0; col<n; col++)
        {
            betta[row][col] = 0.0;
            for (unsigned int s=0; s<L; s++)
            {
                const Condition &c = nonSeparatedConditions().at(s);
                for (unsigned int i=0; i<n; i++) betta[row][col] += c.alpha[row][i] * X(c.time, 0, i, col);
            }
        }
    }

//    Condition lscs;
//    lscs.type = SeparatedLeft;
//    lscs.time = 0.0;
//    lscs.nmbr = 0;
//    lscs.alpha.resize(n1, n);
//    for (unsigned int row=0; row<n1; row++)
//    {
//        for(unsigned int col=0; col<n; col++) lscs.alpha[row][col] = (rand() % 1000) / 1000.0;
//        betta[row+n0] = 0.0;
//        for (unsigned int i=0; i<n; i++) betta[row+n0] += lscs.alpha[row][i] * x(lscs.time, i);
//    }
//    setLeftSeparatedCondition(lscs);

//    Condition rscs;
//    rscs.type = SeparatedRight;
//    rscs.time = 1.0;
//    rscs.nmbr = N;
//    rscs.alpha.resize(n2, n);
//    for (unsigned int row=0; row<n2; row++)
//    {
//        for(unsigned int col=0; col<n; col++) rscs.alpha[row][col] = (rand() % 1000) / 1000.0;
//        betta[row+n0+n1] = 0.0;
//        for (unsigned int i=0; i<n; i++) betta[row+n0+n1] += rscs.alpha[row][i] * x(rscs.time, i);
//    }
//    setRightSeparatedCondition(rscs);

    setBetta(betta);

    IPrinter::print(betta, betta.rows(), betta.cols());
}

void SystemLinearODENonLocalContionsM::calculateForward(DoubleMatrix &x)
{
    n = n0 + n1 + n2;
    unsigned int L = nscs.size();

    //for (unsigned int row=0; row<n0; row++)
    {
        //for (unsigned int col=0; col<n; col++)
        {
            for (unsigned int s=0; s<L-1; s++)
            {
                calculateIntervalF(s,0,0);
            }
        }
    }
}

void SystemLinearODENonLocalContionsM::calculateIntervalF(unsigned int start, unsigned int row, unsigned int col)
{
    unsigned int L = nscs.size();
    double h = grid().dimension().step();
    n = 3;
    DoubleVector x(n+2);
    DoubleVector rx(n+2);

    Condition sc = nscs.at(start);
    Condition ec = nscs.at(start+1);

    for (unsigned int i=0; i<n; i++) x[i] = sc.alpha[row][i];
    x[n] = betta[row][col];
    x[n+1] = 1.0;

    IPrinter::printSeperatorLine();
    IPrinter::print(ec.alpha, ec.alpha.rows(), ec.alpha.cols());
    IPrinter::printSeperatorLine();

    Dimension dim(h, ec.nmbr, sc.nmbr);
    CauchyProblemM1stOrderAM cpa(*this, ODEGrid(dim));
    cpa.row = 0;
    cpa.col = 0;
    IPrinter::print(x, x.size());
    cpa.calculateCP(sc.time, x, rx, InitialValueProblem::RK4);

    IPrinter::print(x,x.size());
    IPrinter::print(rx,rx.size());

    for (unsigned int i=0; i<n; i++) sc.alpha[row][i] = rx[i];
    betta[row][col] = rx[n];
    double M = rx[n+1];

    for (unsigned int s=start+1; s<L; s++)
    {
        Condition &cc = nscs.at(s);
        for (unsigned int i=0; i<n; i++) cc.alpha[row][i] *= M;
    }

    for (unsigned int i=0; i<n; i++) ec.alpha[row][i] += rx[i];

    IPrinter::printSeperatorLine();
    IPrinter::print(ec.alpha, ec.alpha.rows(), ec.alpha.cols());
    IPrinter::printSeperatorLine();


    printf("*** %f %f\n", ec.alpha[cpa.row][0] * X(1.0,0, 0, cpa.col) +
                          ec.alpha[cpa.row][1] * X(1.0,0, 1, cpa.col) +
                          ec.alpha[cpa.row][2] * X(1.0,0, 2, cpa.col),
                          betta[cpa.row][cpa.col]);

    x.clear();
    rx.clear();
}

double SystemLinearODENonLocalContionsM::A(double t, unsigned int, unsigned int row, unsigned int col) const
{
    if (row == 0)
    {
        if (col == 0) { return 1.0; }
        if (col == 1) { return 3.0; }
        if (col == 2) { return 9.0; }
    }
    if (row == 1)
    {
        if (col == 0) { return 8.0; }
        if (col == 1) { return 4.0; }
        if (col == 2) { return 1.0; }
    }
    if (row == 2)
    {
        if (col == 0) { return 0.0; }
        if (col == 1) { return 2.0; }
        if (col == 2) { return 3.0; }
    }
    return NAN;
}

double SystemLinearODENonLocalContionsM::B(double t, unsigned int k, unsigned int row, unsigned int col) const
{
    if (row == 0)
    {
        if (col == 0) { return dX(t,k,0,0) - (A(t,k,0,0)*X(t,k,0,0) + A(t,k,0,1)*X(t,k,1,0) + A(t,k,0,2)*X(t,k,2,0)); }
        if (col == 1) { return dX(t,k,0,1) - (A(t,k,0,0)*X(t,k,0,1) + A(t,k,0,1)*X(t,k,1,1) + A(t,k,0,2)*X(t,k,2,1)); }
        if (col == 2) { return dX(t,k,0,2) - (A(t,k,0,0)*X(t,k,0,2) + A(t,k,0,1)*X(t,k,1,2) + A(t,k,0,2)*X(t,k,2,2)); }
    }
    if (row == 1)
    {
        if (col == 0) { return dX(t,k,1,0) - (A(t,k,1,0)*X(t,k,0,0) + A(t,k,1,1)*X(t,k,1,0) + A(t,k,1,2)*X(t,k,2,0)); }
        if (col == 1) { return dX(t,k,1,1) - (A(t,k,1,0)*X(t,k,0,1) + A(t,k,1,1)*X(t,k,1,1) + A(t,k,1,2)*X(t,k,2,1)); }
        if (col == 2) { return dX(t,k,1,2) - (A(t,k,1,0)*X(t,k,0,2) + A(t,k,1,1)*X(t,k,1,2) + A(t,k,1,2)*X(t,k,2,2)); }
    }
    if (row == 2)
    {
        if (col == 0) { return dX(t,k,2,0) - (A(t,k,2,0)*X(t,k,0,0) + A(t,k,2,1)*X(t,k,1,0) + A(t,k,2,2)*X(t,k,2,0)); }
        if (col == 1) { return dX(t,k,2,1) - (A(t,k,2,0)*X(t,k,0,1) + A(t,k,2,1)*X(t,k,1,1) + A(t,k,2,2)*X(t,k,2,1)); }
        if (col == 2) { return dX(t,k,2,2) - (A(t,k,2,0)*X(t,k,0,2) + A(t,k,2,1)*X(t,k,1,2) + A(t,k,2,2)*X(t,k,2,2)); }
    }
    return NAN;
}

double SystemLinearODENonLocalContionsM::X(double t, unsigned int k, unsigned int row, unsigned int col) const
{
    if (row == 0)
    {
        if (col == 0) { return t; }
        if (col == 1) { return t*t; }
        if (col == 2) { return 2.0*t; }
    }
    if (row == 1)
    {
        if (col == 0) { return 3.0*t; }
        if (col == 1) { return t*t*t; }
        if (col == 2) { return 5.0*t; }
    }
    if (row == 2)
    {
        if (col == 0) { return -3.0*t; }
        if (col == 1) { return t*t+t; }
        if (col == 2) { return 8.0*t; }
    }
    return NAN;
}

double SystemLinearODENonLocalContionsM::dX(double t, unsigned int k, unsigned int row, unsigned int col) const
{
    if (row == 0)
    {
        if (col == 0) { return 1.0; }
        if (col == 1) { return 2.0*t; }
        if (col == 2) { return 2.0; }
    }
    if (row == 1)
    {
        if (col == 0) { return 3.0; }
        if (col == 1) { return 3.0*t*t; }
        if (col == 2) { return 5.0; }
    }
    if (row == 2)
    {
        if (col == 0) { return -3.0; }
        if (col == 1) { return 2.0*t+1.0; }
        if (col == 2) { return 8.0; }
    }
    return NAN;
}

double CauchyProblemM1stOrderAM::f(double t, const DoubleVector &x, unsigned int k, unsigned int i) const
{
    unsigned int n = 3;//p.systemOrder();

    double _SO = S0(t,x,k);

    if (i<n)
    {
        double res = _SO*x[i];
        for (unsigned int j=0; j<n; j++) res -= p.A(t,k,j,i)*x[j];
        return res;
    }
    else if (i==n)
    {
        double res = _SO*x[n];
        for (unsigned int j=0; j<n; j++) res += x[j]*p.B(t,k,j,col);
        return res;
    }
    else
    {
        return _SO*x[n+1];
    }

    return NAN;
}

double CauchyProblemM1stOrderAM::S0(double t, const DoubleVector &x, unsigned int k) const
{
    unsigned int n = 3;//p.systemOrder();
    double btr = x[n];

    double s1 = 0.0;
    for (unsigned int i=0; i<n; i++)
    {
        double aa = 0.0;
        for (unsigned int j=0; j<n; j++) aa += x[j]*p.A(t,k,j,i);
        s1 += aa*x[i];
    }

    double s2 = 0.0;
    for (unsigned int i=0; i<n; i++)
    {
        s2 += x[i]*p.B(t,k,i,col);
    }
    s2 *= btr;


    double m1 = 0.0;
    for (unsigned int i=0; i<n; i++)
    {
        m1 += x[i]*x[i];
    }
    m1 += btr*btr;

    return (s1-s2)/m1;
}

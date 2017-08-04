#include "islodenlcs.h"

class CauchyProblemM1stOrderA : public CauchyProblemM1stOrder
{
public:
    CauchyProblemM1stOrderA(ISystemLinearODENonLocalContions &parent, const ODEGrid& grid) : CauchyProblemM1stOrder(grid), p(parent) {}

protected:
    virtual double f(double t, const DoubleVector &x, unsigned int k, unsigned int i) const;
    double S0(double t, const DoubleVector &x, unsigned int k) const;
private:
    ISystemLinearODENonLocalContions &p;
};

class CauchyProblemM1stOrderB : public CauchyProblemM1stOrder
{
public:
    CauchyProblemM1stOrderB(ISystemLinearODENonLocalContions &parent, const ODEGrid& grid) : CauchyProblemM1stOrder(grid), p(parent) {}
protected:
    virtual double f(double t, const DoubleVector &x, unsigned int k, unsigned int i) const;
private:
    ISystemLinearODENonLocalContions &p;
};

ISystemLinearODENonLocalContions::ISystemLinearODENonLocalContions(const ODEGrid &grid) : SystemLinearODE1stOrder(grid) {}

void ISystemLinearODENonLocalContions::setLeftSeparatedCondition(const Condition &lscs)
{
    this->lscs = lscs;
    this->n1 = lscs.alpha.rows();
}

void ISystemLinearODENonLocalContions::setRightSeparatedCondition(const Condition &rscs)
{
    this->rscs = rscs;
    this->n2 = rscs.alpha.rows();
}

void ISystemLinearODENonLocalContions::addNonSeparatedCondition(const Condition &nsc)
{
    this->nscs.push_back(nsc);
    this->n0 = nsc.alpha.rows();
    this->L++;
}

const std::vector<ISystemLinearODENonLocalContions::Condition>& ISystemLinearODENonLocalContions::nonSeparatedConditions() const
{
    return nscs;
}

void ISystemLinearODENonLocalContions::setBetta(const DoubleVector &betta)
{
    this->betta = betta;
}

void ISystemLinearODENonLocalContions::calculateIntervalF(unsigned int start, unsigned int r)
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

void ISystemLinearODENonLocalContions::calculateForward(DoubleVector &x)
{
    n = n0 + n1 + n2;
    L = nscs.size();

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

void ISystemLinearODENonLocalContions::calculateBackward(DoubleVector &x)
{}

void ISystemLinearODENonLocalContions::calculateBackwardCP(const DoubleVector &x, DoubleMatrix &m)
{
    CauchyProblemM1stOrderB cpb(*this, grid());
    cpb.calculateCP(nscs.back().time, x, m, CauchyProblemM1stOrder::RK4, CauchyProblemM1stOrder::R2L);
}

unsigned int ISystemLinearODENonLocalContions::systemOrder() const
{
    return n;
}

double CauchyProblemM1stOrderA::f(double t, const DoubleVector &x, unsigned int k, unsigned int i) const
{
    unsigned int n = p.systemOrder();
    double _SO = S0(t,x,k);

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

    return NAN;
}

double CauchyProblemM1stOrderA::S0(double t, const DoubleVector &x, unsigned int k) const
{
    unsigned int n = p.systemOrder();
    double btr = x[n];

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
    unsigned int n = p.systemOrder();
    double res = p.B(t,k,i+1);
    for (unsigned int j=1; j<=n; j++) res += p.A(t,k,i+1,j)*x[j-1];
    return res;
}

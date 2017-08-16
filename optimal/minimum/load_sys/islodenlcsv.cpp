#include "islodenlcsv.h"

class CauchyProblemM1stOrderA : public CauchyProblemM1stOrder
{
public:
    CauchyProblemM1stOrderA(ISystemLinearODENonLocalContionsV &parent, const ODEGrid& grid) : CauchyProblemM1stOrder(grid), p(parent) {}

protected:
    virtual double f(double t, const DoubleVector &x, unsigned int k, unsigned int i) const;
    double S0(double t, const DoubleVector &x, unsigned int k) const;
private:
    ISystemLinearODENonLocalContionsV &p;
};

class CauchyProblemM1stOrderB : public CauchyProblemM1stOrder
{
public:
    CauchyProblemM1stOrderB(ISystemLinearODENonLocalContionsV &parent, const ODEGrid& grid) : CauchyProblemM1stOrder(grid), p(parent) {}
protected:
    virtual double f(double t, const DoubleVector &x, unsigned int k, unsigned int i) const;
private:
    ISystemLinearODENonLocalContionsV &p;
};

ISystemLinearODENonLocalContionsV::ISystemLinearODENonLocalContionsV(const ODEGrid &grid) : SystemLinearODE1stOrder(grid) {}

void ISystemLinearODENonLocalContionsV::setLeftSeparatedCondition(const Condition &lscs)
{
    this->lscs = lscs;
}

void ISystemLinearODENonLocalContionsV::setRightSeparatedCondition(const Condition &rscs)
{
    this->rscs = rscs;
}

void ISystemLinearODENonLocalContionsV::addNonSeparatedCondition(const Condition &nsc)
{
    this->nscs.push_back(nsc);
}

const std::vector<ISystemLinearODENonLocalContionsV::Condition>& ISystemLinearODENonLocalContionsV::nonSeparatedConditions() const
{
    return nscs;
}

void ISystemLinearODENonLocalContionsV::setBetta(const DoubleVector &betta)
{
    this->betta = betta;
}

void ISystemLinearODENonLocalContionsV::calculateIntervalF(unsigned int start, unsigned int r)
{
    unsigned int L = nscs.size();
    unsigned int n0 = 0;
    if (L > 0) n0 = nscs.at(0).alpha.rows();
    unsigned int n1 = lscs.alpha.rows();
    unsigned int n2 = rscs.alpha.rows();
    unsigned int n = n0 + n1 + n2;

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

void ISystemLinearODENonLocalContionsV::calculateForward(DoubleVector &x)
{
    unsigned int L = nscs.size();

    unsigned int n0 = 0;
    if (L > 0) n0 = nscs.at(0).alpha.rows();

    unsigned int n1 = lscs.alpha.rows();

    unsigned int n2 = rscs.alpha.rows();

    unsigned int n = n0 + n1 + n2;

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

void ISystemLinearODENonLocalContionsV::calculateBackward(DoubleVector &x UNUSED_PARAM)
{}

void ISystemLinearODENonLocalContionsV::calculateBackwardCP(const DoubleVector &x, std::vector<DoubleVector>& m)
{
    CauchyProblemM1stOrderB cpb(*this, grid());
    cpb.calculateCP(nscs.back().time, x, m, CauchyProblemM1stOrder::RK4, CauchyProblemM1stOrder::R2L);
}

unsigned int ISystemLinearODENonLocalContionsV::systemOrder() const
{
    unsigned int L = nscs.size();
    unsigned int n0 = 0;
    if (L > 0) n0 = nscs.at(0).alpha.rows();
    unsigned int n1 = lscs.alpha.rows();
    unsigned int n2 = rscs.alpha.rows();
    unsigned int n = n0 + n1 + n2;
    return n;
}

double CauchyProblemM1stOrderA::f(double t, const DoubleVector &x, unsigned int k, unsigned int i) const
{
    unsigned int n = p.systemOrder();
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
        for (unsigned int j=0; j<n; j++) res += p.B(t,k,j)*x[j];
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
    for (unsigned int i=0; i<n; i++)
    {
        double aa = 0.0;
        for (unsigned int j=0; j<n; j++) aa += x[j]*p.A(t,k,j,i);
        s1 += aa*x[i];
    }

    double s2 = 0.0;
    for (unsigned int i=0; i<n; i++)
    {
        s2 += x[i]*p.B(t,k,i);
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

double CauchyProblemM1stOrderB::f(double t, const DoubleVector &x, unsigned int k, unsigned int i) const
{
    unsigned int n = p.systemOrder();
    double res = p.B(t,k,i);
    for (unsigned int j=0; j<n; j++) res += p.A(t,k,i,j)*x[j];
    return res;
}

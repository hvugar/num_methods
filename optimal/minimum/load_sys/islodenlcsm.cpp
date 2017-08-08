#include "islodenlcsm.h"
#include <grid/cauchyp.h>

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

ISystemLinearODENonLocalContionsM::ISystemLinearODENonLocalContionsM(const ODEGrid& grid) : mgrid(grid)
{
}

void ISystemLinearODENonLocalContionsM::calculateForward(DoubleVector &x)
{
}

void ISystemLinearODENonLocalContionsM::calculateIntervalF(unsigned int start, unsigned int row, unsigned int col)
{
    unsigned int L = nscs.size();
    double h = grid().dimension().step();
    DoubleVector x(n+2);
    DoubleVector rx(n+2);

    Condition &sc = nscs.at(start);
    Condition &ec = nscs.at(start+1);

    for (unsigned int i=0; i<n; i++)
        x[i] = sc.alpha[row][i];
    x[n] = betta[row][col];
    x[n+1] = 1.0;

    Dimension dim(h, ec.nmbr, sc.nmbr);
    CauchyProblemM1stOrderAM cpa(*this, ODEGrid(dim));
    cpa.row = row;
    cpa.col = col;
    cpa.calculateCP(sc.time, x, rx, InitialValueProblem::RK4);

    for (unsigned int i=0; i<n; i++) sc.alpha[row][i] = rx[i];
    betta[row][col] = rx[n];
    double M = rx[n+1];

    for (unsigned int s=start+1; s<L; s++)
    {
        Condition &cc = nscs.at(s);
        for (unsigned int i=0; i<n; i++) cc.alpha[row][i] *= M;
    }

    for (unsigned int i=0; i<n; i++) ec.alpha[row][i] += rx[i];

    x.clear();
    rx.clear();
}

unsigned int ISystemLinearODENonLocalContionsM::systemOrder() const
{
    return n;
}

const ODEGrid& ISystemLinearODENonLocalContionsM::grid() const
{
    return mgrid;
}

double CauchyProblemM1stOrderAM::f(double t, const DoubleVector &x, unsigned int k, unsigned int i) const
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
        for (unsigned int j=0; j<n; j++) res += x[j]*p.B(t,k,j+1,col);
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
        //s2 += x[i-1]*p.B(t,k,i);
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


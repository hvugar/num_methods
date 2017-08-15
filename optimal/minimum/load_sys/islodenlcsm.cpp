#include "islodenlcsm.h"

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

class CauchyProblemM1stOrderBM : public CauchyProblemM1stOrder
{
public:
    CauchyProblemM1stOrderBM(ISystemLinearODENonLocalContionsM &parent, const ODEGrid& grid) : CauchyProblemM1stOrder(grid), p(parent) {}

    unsigned int row;
    unsigned int col;
    unsigned int n;

protected:
    virtual double f(double t, const DoubleVector &x, unsigned int k, unsigned int i) const;
private:
    ISystemLinearODENonLocalContionsM &p;
};

ISystemLinearODENonLocalContionsM::ISystemLinearODENonLocalContionsM(const ODEGrid& grid) : mgrid(grid)
{
}

void ISystemLinearODENonLocalContionsM::calculateForward(DoubleMatrix &x)
{
    unsigned int L = nscs.size();

    unsigned int n0 = 0;
    if (L > 0) n0 = nscs.at(0).alpha.rows();

    unsigned int n1 = lscs.alpha.rows();

    unsigned int n2 = rscs.alpha.rows();

    unsigned int n = n0 + n1 + n2;

    double h = grid().dimension().step();

    DoubleVector ix(n+2);
    DoubleVector ox(n+2);
    DoubleMatrix bt(n0, n);

    DoubleMatrix A(n0*n, n*n);
    DoubleVector b(n0*n);

    for (unsigned int row=0; row<n0; row++)
    {
        for (unsigned int col=0; col<n; col++)
        {
            DoubleMatrix alpha(L, n);
            for (unsigned int i = 0; i<L; i++)
            {
                for (unsigned int j = 0; j<n; j++)
                {
                    alpha[i][j] = nscs.at(i).alpha[row][j];
                }
            }
            bt[row][col] = betta[row][col];

            for (unsigned int start = 0; start<L-1; start++)
            {
                ISystemLinearODENonLocalContions::Condition sc = nscs.at(start);
                ISystemLinearODENonLocalContions::Condition ec = nscs.at(start+1);

                for (unsigned int i=0; i<n; i++) ix[i] = alpha[start][i];
                ix[n+0] = bt[row][col];
                ix[n+1] = 1.0;

                Dimension dim(h, ec.nmbr, sc.nmbr);
                CauchyProblemM1stOrderAM cpa(*this, ODEGrid(dim));
                cpa.row = row;
                cpa.col = col;
                cpa.calculateCP(sc.time, ix, ox, InitialValueProblem::RK4);

                // Assign alpha[0] at point t1
                for (unsigned int i=0; i<n; i++) alpha[start][i] = ox[i];
                // Assign betta at point t1
                bt[row][col] = ox[n];
                double M = ox[n+1];

                // Assign alpha[1,..,L] at point t1
                for (unsigned int s=start+1; s<L; s++)
                {
                    for (unsigned int i=0; i<n; i++) alpha[s][i] *= M;
                }

                // Adding alpha[0] + alpha[1] at point t1
                for (unsigned int i=0; i<n; i++) alpha[start+1][i] += alpha[start][i];
            }

            b[row*n+col] = bt[row][col];

            for (unsigned int i=0; i<n*n; i++)
            {
                if (i%n == (row*n+col)%n)
                {
                    A[row*n+col][i] = alpha[L-1][i/n];
                }
            }
            alpha.clear();
        }
    }

    for (unsigned int row=0; row<n1; row++)
    {
//        unsigned int minN = grid().dimension().minN();
//        unsigned int maxN = grid().dimension().maxN();

//        DoubleVector x(n+2);
//        DoubleVector rx(n+2);

//        for (unsigned int i=0; i<n; i++) x[i] = lscs.alpha[row][i]; x[n] = betta[row+n0]; x[n+1] = 1.0;

//        Dimension dim(h, maxN, minN);
//        CauchyProblemM1stOrderA cpa(*this, ODEGrid(dim));
//        cpa.calculateCP(lscs.time, x, rx, InitialValueProblem::RK4);

//        for (unsigned int i=0; i<n; i++) lscs.alpha[row][i] = rx[i];
//        betta[row+n0] = rx[n];
//        x.clear();
//        rx.clear();
    }

    for (unsigned int row=0; row<n2; row++)
    {}

    DoubleVector m(n*n);
    GaussianElimination(A, b, m);

    x.resize(n, n);
    for (unsigned int i=0; i<n; i++)
    {
        for (unsigned int j=0; j<n; j++)
        {
            x[i][j] = m[i*n+j];
        }
    };

    ix.clear();
    ox.clear();
}

void ISystemLinearODENonLocalContionsM::calculateBackwardCP(DoubleMatrix &x, DoubleVector** m)
{
    unsigned int n = systemOrder();
    DoubleVector y;
    for (unsigned int i=0; i<x.rows(); i++)
    {
        for (unsigned int j=0; j<x.cols(); j++)
        {
            y << x[i][j];
        }
    }
    CauchyProblemM1stOrderBM cpb(*this, grid());
    cpb.n = systemOrder();
    DoubleMatrix m1;
    cpb.calculateCP(1.0, y, m1, CauchyProblemM1stOrder::RK4, CauchyProblemM1stOrder::R2L);

    m = new DoubleVector*[n];
    for (unsigned int i=0; i<n; i++)
    {
        m[i] = new DoubleVector[n];
        for (unsigned int j=0; j<n; j++)
            m[i][j] = m1.row(i*n+j);
    }
}

void ISystemLinearODENonLocalContionsM::calculateBackward(DoubleMatrix &x UNUSED_PARAM)
{}

unsigned int ISystemLinearODENonLocalContionsM::systemOrder() const
{
    unsigned int L = nscs.size();
    unsigned int n0 = 0;
    if (L > 0) n0 = nscs.at(0).alpha.rows();
    unsigned int n1 = lscs.alpha.rows();
    unsigned int n2 = rscs.alpha.rows();
    unsigned int n = n0 + n1 + n2;
    return n;
}

const ODEGrid& ISystemLinearODENonLocalContionsM::grid() const
{
    return mgrid;
}

void ISystemLinearODENonLocalContionsM::setLeftSeparatedCondition(const ISystemLinearODENonLocalContions::Condition &lscs)
{
    C_UNUSED(lscs);
}

void ISystemLinearODENonLocalContionsM::setRightSeparatedCondition(const ISystemLinearODENonLocalContions::Condition &rscs)
{
    C_UNUSED(rscs);
}

void ISystemLinearODENonLocalContionsM::addNonSeparatedCondition(const ISystemLinearODENonLocalContions::Condition &nsc)
{
    this->nscs.push_back(nsc);
}

const std::vector<ISystemLinearODENonLocalContions::Condition>& ISystemLinearODENonLocalContionsM::nonSeparatedConditions() const
{
    return nscs;
}

void ISystemLinearODENonLocalContionsM::setBetta(const DoubleMatrix &betta)
{
    this->betta = betta;
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

double CauchyProblemM1stOrderBM::f(double t, const DoubleVector &x, unsigned int k, unsigned int i) const
{
    int row = i / n;//1
    int col = i % n;//2

    unsigned int n = p.systemOrder();
    double res = p.B(t,k,row,col);
    for (unsigned int j=0; j<n; j++) res += p.A(t,k,row,j)*x[n*j+col];
    return res;
}

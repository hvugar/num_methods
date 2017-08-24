#include "lode1o.h"
#include "nlode1o.h"
#include <math.h>
#include <cmethods.h>

void LinearODE1stOrder::calculate(const std::vector<Condition> &nscs, const DoubleVector &bt, std::vector<DoubleVector> &x)
{
    double h = grid().dimension().step();

    std::vector<Condition> cs = nscs;
    DoubleVector beta = bt;

    unsigned int L = cs.size();
    unsigned int n = 0;
    if (L!=0) n = cs.at(0).mtrx.rows();

    DoubleVector x0(n+2);
    DoubleVector rx(n+2);
    for (unsigned int row=0; row<n; row++)
    {
        for (unsigned int s=0; s<L-1; s++)
        {
            Condition &sc = cs.at(s);
            Condition &ec = cs.at(s+1);

            for (unsigned int i=0; i<n; i++) x0[i] = sc.mtrx[row][i];
            x0[n] = beta[row];
            x0[n+1] = 1.0;

            struct HelperB : public NonLinearODE1stOrder
            {
                LinearODE1stOrder *p;
                unsigned int n;

                virtual double f(double, double, unsigned int) const { return NAN; }

                virtual double f(double t, const DoubleVector &x, unsigned int k, unsigned int i) const
                {
                    //unsigned int n = p->systemOrder();
                    double _SO = S0(t,x,k);

                    if (i<n)
                    {
                        double res = _SO*x[i];
                        for (unsigned int j=0; j<n; j++) res -= p->A(t,k,j,i)*x[j];
                        return res;
                    }
                    else if (i==n)
                    {
                        double res = _SO*x[n];
                        for (unsigned int j=0; j<n; j++) res += p->B(t,k,j)*x[j];
                        return res;
                    }
                    else
                    {
                        return _SO*x[n+1];
                    }

                    return NAN;
                }

                double S0(double t, const DoubleVector &x, unsigned int k) const
                {
                    //unsigned int n = p.systemOrder();
                    double btr = x[n];

                    double s1 = 0.0;
                    for (unsigned int i=0; i<n; i++)
                    {
                        double aa = 0.0;
                        for (unsigned int j=0; j<n; j++) aa += x[j]*p->A(t,k,j,i);
                        s1 += aa*x[i];
                    }

                    double s2 = 0.0;
                    for (unsigned int i=0; i<n; i++)
                    {
                        s2 += x[i]*p->B(t,k,i);
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
            };

            HelperB helper;
            helper.p = this;
            helper.n = n;
            helper.setGrid(ODEGrid(Dimension(h, ec.nmbr, sc.nmbr)));
            helper.cauchyProblem(sc.time, x0, rx, NonLinearODE1stOrder::RK4);

            for (unsigned int i=0; i<n; i++) sc.mtrx[row][i] = rx[i];
            beta[row] = rx[n];
            double M = rx[n+1];

            for (unsigned int j=s+1; j<L; j++)
            {
                Condition &cc = cs.at(j);
                for (unsigned int i=0; i<n; i++) cc.mtrx[row][i] *= M;
            }

            for (unsigned int i=0; i<n; i++) ec.mtrx[row][i] += rx[i];
        }
    }
    x.clear();
    rx.clear();

    // Separated conditions in right side

    DoubleMatrix A(n, n);
    DoubleVector b(n);
    DoubleVector x1(n);

    Condition c0 = cs.at(L-1);
    for (unsigned int row=0; row<n; row++)
    {
        for (unsigned int col=0; col<n; col++)
        {
            A[row][col] = c0.mtrx[row][col];
        }
        b[row] = beta[row];
    }

    GaussianElimination(A, b, x1);

    struct HelperB : public NonLinearODE1stOrder
    {
        LinearODE1stOrder *p;
        unsigned int n;
    protected:
        virtual double f(double t, const DoubleVector &x, unsigned int k, unsigned int i) const
        {
            double res = p->B(t,k,i);
            for (unsigned int j=0; j<n; j++) res += p->A(t,k,i,j)*x[j];
            return res;
        }
    };

    HelperB helper;
    helper.p = this;
    helper.n = n;
    helper.setGrid(grid());
    helper.cauchyProblem(nscs.back().time, x1, x, HelperB::RK4, HelperB::R2L);
}

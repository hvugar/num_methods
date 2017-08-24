#include "lode1o.h"
#include "nlode1o.h"
#include <math.h>
#include <cmethods.h>
#include <float.h>

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

void LinearODE1stOrder::highOder2Accuracy(const std::vector<Condition> &cnds, const DoubleVector& rs)
{
    unsigned int n =  equationsNumber();

    double h = grid().dimension().step();
    unsigned int N = grid().dimension().sizeN();

    if (n == 1)
    {
        DoubleVector betta(N+1);
        for (unsigned int m=0; m<=N; m++) betta[m] = 0.0;

        for (unsigned int i=0; i<cnds.size(); i++)
        {
            const Condition &c = cnds[i];
            double alpha = c.mtrx.at(0,0);
            double time  = c.time;

            for (unsigned int m=0; m<=N; m++)
            {
                double dh = fabs(time - m*h);

                //printf("%f %f %f\n", alpha, time, dh);
                if (dh <= h)
                {
                    betta[m] += alpha*(1.0 - dh/h);
                }
            }
        }
        IPrinter::printVector(14, 10, betta);

        std::vector<unsigned int> indexes;
        for (unsigned int m=1; m<=N; m++) if (betta[m] != 0.0) indexes.push_back(m);

        DoubleVector f(N+1);
        std::vector<DoubleVector> ems(indexes.size()); for (unsigned int i=0; i<ems.size(); i++) ems[i].resize(N+1);
    }
    else
    {
//        DoubleMatrix* alhas = new DoubleMatrix[N+1];
//        for (unsigned int i=0; i<=N; i++) alhas[0].resize(n,n);

//        for (unsigned int s=0; s<cnds.size(); s++)
//        {
//            const Condition &cnd = cnds.at(s);
//            alhas[i] = cnd.m;
//        }

//        unsigned int next = 1;
//        DoubleMatrix betta1 = cnds.at(0).mtrx;
//        DoubleMatrix betta2(n,n);
//        DoubleVector gamma(n);

//        DoubleVector alpha0;
//        DoubleMatrix alpha1;
//        DoubleMatrix alpha2;

//        for (unsigned int i=1; i<N-2; i++)
//        {
//            gamma -= betta1*alpha0;
//            betta2 = betta1*alpha2;
//            betta1 = betta1*alpha1;
//            if (cnds.at(next).n==i)
//            {
//                betta1 += cnds.at(next).m;
//                next++;
//            }
//            if (cnds.at(next).n==i+1)
//            {
//                betta2 += cnds.at(next).m;
//                next++;
//            }
//        }
    }
}

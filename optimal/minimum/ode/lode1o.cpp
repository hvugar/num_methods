#include "lode1o.h"
#include "nlode1o.h"
#include "../matrix2d.h"
#include "../printer.h"
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

    LinearEquation::GaussianElimination(A, b, x1);

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

void LinearODE1stOrder::highOderAccuracy(const std::vector<Condition> &cs, const DoubleVector &rs, std::vector<DoubleVector> &x, unsigned int k, Direction direction UNUSED_PARAM)
{
    switch (k)
    {
    case 2: highOder2Accuracy(cs, rs, x, direction); break;
    case 4: highOder4Accuracy(cs, rs, x, direction); break;
    case 6: highOder6Accuracy(cs, rs, x, direction); break;
    default: break;
    }
}

void LinearODE1stOrder::discretisation(const std::vector<Condition>& cs, double* b) const
{
    double h = grid().dimension().step();
    unsigned int N = grid().dimension().sizeN();

    unsigned int cnd_size = cs.size();
    for (unsigned int s=0; s<cnd_size; s++)
    {
        const Condition &cnd = cs[s];
        double alpha = cnd.mtrx.at(0,0);
        double time  = cnd.time;

        for (unsigned int n=0; n<=N; n++)
        {
            double dh = fabs(time - n*h);

            if (dh <= h)
            {
                b[n] += alpha*(1.0 - dh/h);
            }
        }
    }
}

void LinearODE1stOrder::highOder2Accuracy(const std::vector<Condition> &cnds, const DoubleVector & rs, std::vector<DoubleVector> &x, Direction direction UNUSED_PARAM)
{
    unsigned int en = equationsNumber();

    double h = grid().dimension().step();
    unsigned int N = grid().dimension().sizeN();

    if (en == 1)
    {
        double *p = (double*) malloc(sizeof(double)*N);
        double *q = (double*) malloc(sizeof(double)*N);
        double *r = (double*) malloc(sizeof(double)*N);

        std::vector<unsigned int> ind;
        std::vector<double*> ems;

        /**********************************************************************
         *                  Discretization using trapesium rule
         *********************************************************************/
        double *b = (double*) malloc(sizeof(double)*(N+1));
        for (unsigned int m=0; m<=N; m++) b[m] = 0.0;
        discretisation(cnds, b);

        p[0] = b[0];
        q[0] = b[1];
        for (unsigned int n=2; n<=N; n++)
        {
            if (fabs(b[n]) > DBL_EPSILON)
            {
                ind.push_back(n);
                ems.push_back((double*) malloc(sizeof(double)*N));
                ems[ems.size()-1][0] = b[n];
            }
        }
        free(b);
        r[0] = rs[0];

        unsigned int ind_size = ind.size();
        unsigned int ems_size = ems.size();

        /**********************************************************************
         *                          End of discretization
         *********************************************************************/

        /**********************************************************************
         *                      Finding function at end of grid
         *********************************************************************/

        for (unsigned int n=0; n<=N-2; n++)
        {
            double t = n*h;
            double m = 1.0/(2.0*h*A(t,n)+3.0);
            double alpha1 = +4.0*m;
            double alpha2 = -1.0*m;
            double alpha0 = -2.0*h*B(t,n)*m;

            r[n+1] = r[n] - p[n]*alpha0;
            q[n+1] = p[n]*alpha2;
            p[n+1] = p[n]*alpha1 + q[n];


            for (unsigned int i=0; i<ind_size; i++)
            {
                if (n+2 == ind.at(i))
                {
                    q[n+1] += ems[i][n];
                }
                if (n+2 >= ind.at(i))
                {
                    ems[i][n+1] = 0.0;
                }
                else
                {
                    ems[i][n+1] = ems[i][n];
                }
            }
        }

        DoubleMatrix m(3,3);
        DoubleVector c(3);
        DoubleVector xT(3);

        for (unsigned int row=0; row<3*en; row++)
        {
            for (unsigned int col=0; col<3*en; col++)
            {

            }
        }
        m[0][0] = 0.0;
        m[0][1] = p[N-1];
        m[0][2] = q[N-1];
        c[0] = r[N-1];

        m[1][0] = -1.0;
        m[1][1] = -2.0*h*A((N-1)*h,N-1);
        m[1][2] = +1.0;
        c[1] = 2.0*h*B((N-1)*h,N-1);

        m[2][0] = +1.0;
        m[2][1] = -4.0;
        m[2][2] = 3.0-2.0*h*A(N*h,N);
        c[2] = 2.0*h*B(N*h,N);

        LinearEquation::GaussianElimination(m, c, xT);
        m.clear();
        c.clear();

        /**********************************************************************
         *                      Finding function at end of grid
         *********************************************************************/

        DoubleVector x0(N+1);
        x0[N-0] = xT[2];
        x0[N-1] = xT[1];
        x0[N-2] = xT[0];
        xT.clear();


        for (unsigned int i=N-2; i!=0; i--)
        {
            x0[i-1] = -q[i-1]*x0[i]+r[i-1];
            for (unsigned int s=0; s<ems_size; s++) x0[i-1] -= ems[s][i-1]*x0[ind[s]];
            x0[i-1] /= p[i-1];
        }

        x.clear();
        x.push_back(x0);

        free(p);
        free(q);
        free(r);
        for (unsigned int i=0; i<ems_size; i++) free(ems[i]);
        ems.clear();
        ind.clear();
    }
    else
    {
        DoubleMatrix* p = new DoubleMatrix[N]; for (unsigned int i=0; i<N; i++) p[i].resize(en, en, 0.0);
        DoubleMatrix* q = new DoubleMatrix[N]; for (unsigned int i=0; i<N; i++) q[i].resize(en, en, 0.0);
        DoubleMatrix* r = new DoubleMatrix[N]; for (unsigned int i=0; i<N; i++) r[i].resize(en, 1, 0.0);

        std::vector<unsigned int> ind;
        std::vector<DoubleMatrix*> ems;

        /**********************************************************************
         *                  Discretization using trapesium rule
         *********************************************************************/
        DoubleMatrix* b = new DoubleMatrix[N+1]; for (unsigned int n=0; n<=N; n++) b[n].resize(en, en, 0.0);

        unsigned int cnd_size = cnds.size();
        for (unsigned int s=0; s<cnd_size; s++)
        {
            const Condition &cnd = cnds[s];
            double time  = cnd.time;

            for (unsigned int n=0; n<=N; n++)
            {
                double dh = fabs(time - n*h);
                if (dh <= h)
                {
                    b[n] += (1.0 - dh/h)*cnd.mtrx;
                }
            }
        }

        p[0] = b[0];
        q[0] = b[1];
        for (unsigned int n=2; n<=N; n++)
        {
            if (!b[n].zeroMatrix())
            {
                ind.push_back(n);
                DoubleMatrix *dm = new DoubleMatrix[N];
                for (unsigned int i=0; i<N; i++) dm[i].resize(en,en,0.0);
                ems.push_back(dm);
                ems[ems.size()-1][0] = b[n];
            }
        }
        free(b);
        r[0] = rs;

        unsigned int ind_size = ind.size();
        unsigned int ems_size = ems.size();

        //printf("%d %d %d %d\n", ind[0], ind[1], ind[2], ind[3]);

        /**********************************************************************
         *                          End of discretization
         *********************************************************************/

        /**********************************************************************
         *                      Finding function at end of grid
         *********************************************************************/

        DoubleMatrix m(en, en);
        DoubleMatrix alpha1(en, en, 0.0);
        DoubleMatrix alpha2(en, en, 0.0);
        DoubleMatrix alpha0(en, 1, 0.0);

        for (unsigned int n=0; n<=N-2; n++)
        {
            double t = n*h;

            for (unsigned int row=0; row<en; row++)
            {
                for (unsigned int col=0; col<en; col++)
                {
                    m[row][col] = 2.0*h*A(t,n,row,col);
                    if (row==col) m[row][col] += 3.0;

                    alpha1[row][col] = 0.0;
                    if (row==col) alpha1[row][col] = 4.0;

                    alpha2[row][col] = 0.0;
                    if (row==col) alpha2[row][col] = -1.0;
                }
                alpha0[row][0] = -2.0*h*B(t,n,row);
            }

            m.inverse();

            alpha1 = m*alpha1;
            alpha2 = m*alpha2;
            alpha0 = m*alpha0;

            r[n+1] = r[n] - p[n]*alpha0;
            q[n+1] = p[n]*alpha2;
            p[n+1] = p[n]*alpha1 + q[n];

            for (unsigned int i=0; i<ind_size; i++)
            {
                if (n+2 == ind.at(i))
                {
                    q[n+1] += ems[i][n];
                    //ems[i][n] = DoubleMatrix::ZeroMatrix(en,en);
                }

                if (n+2 >= ind.at(i))
                {
                    ems[i][n+1] = DoubleMatrix::ZeroMatrix(en,en);
                }
                else
                {
                    ems[i][n+1] = ems[i][n];
                }
            }
        }

        m.clear();
        alpha0.clear();
        alpha1.clear();
        alpha2.clear();

        DoubleMatrix M(3*en,3*en);
        DoubleVector C(3*en);
        DoubleVector xT(3*en);

        for (unsigned int row=0*en; row<1*en; row++)
        {
            for (unsigned int col=0*en; col<1*en; col++) M[row][col] = 0.0;
            for (unsigned int col=1*en; col<2*en; col++) M[row][col] = p[N-1][row%en][col%en];
            for (unsigned int col=2*en; col<3*en; col++) M[row][col] = q[N-1][row%en][col%en];
            C[row] = r[N-1][row%en][0];
        }

        for (unsigned int row=1*en; row<2*en; row++)
        {
            for (unsigned int col=0*en; col<1*en; col++) { M[row][col] = +0.0; if (row%en == col%en) M[row][col] += -1.0; }
            for (unsigned int col=1*en; col<2*en; col++) { M[row][col] = -2.0*h*A((N-1)*h, N-1, row%en, col%en); }
            for (unsigned int col=2*en; col<3*en; col++) { M[row][col] = +0.0; if (row%en == col%en) M[row][col] += +1.0; }
            C[row] = 2.0*h*B((N-1)*h, N-1, row%en);
        }

        for (unsigned int row=2*en; row<3*en; row++)
        {
            for (unsigned int col=0*en; col<1*en; col++) { M[row][col] = +0.0; if (row%en == col%en) M[row][col] += +1.0; }
            for (unsigned int col=1*en; col<2*en; col++) { M[row][col] = +0.0; if (row%en == col%en) M[row][col] += -4.0; }
            for (unsigned int col=2*en; col<3*en; col++) { M[row][col] = -2.0*h*A(N*h, N, row%en, col%en); if (row%en == col%en) M[row][col] += +3.0; }
            C[row] = 2.0*h*B(N*h, N, row%en);
        }

        LinearEquation::GaussianElimination(M, C, xT);
        //IPrinter::print(xT, xT.size());

        x.resize(en);
        for (unsigned int i=0; i<en; i++) x[i].resize(N+1);

        std::vector<DoubleMatrix> x0(N+1);
        for (unsigned int i=0; i<=N; i++) x0[i].resize(en,1);

        for (unsigned int row=0; row<en; row++)
        {
            x0[N-2][row][0] = xT[row+0*en];
            x0[N-1][row][0] = xT[row+1*en];
            x0[N-0][row][0] = xT[row+2*en];
        }

        for (unsigned int i=N-2; i!=0; i--)
        {
            x0[i-1] = r[i-1]-q[i-1]*x0[i];
            for (unsigned int s=0; s<ems_size; s++)
            {
                x0[i-1] -= (ems[s][i-1]*x0[ind[s]]);
            }
            p[i-1].inverse();
            x0[i-1] = p[i-1] * x0[i-1];
        }

        for (unsigned int i=0; i<=N; i++)
        {
            for (unsigned int row=0; row<en; row++)
            {
                x[row][i] = x0[i][row][0];
            }
        }

        for (unsigned int i=0; i<N; i++)
        {
            r[i].clear();
            q[i].clear();
            p[i].clear();
        }
        delete [] r;
        delete [] q;
        delete [] p;
        ind.clear();

        for (unsigned int i=0; i<ems_size; i++)
        {
            delete [] ems[i];
        }
        ems.clear();

        for (unsigned int i=0; i<=N; i++)
        {
            x0[i].clear();
        }
        x0.clear();

        M.clear();
        C.clear();
        xT.clear();
    }
}

void LinearODE1stOrder::highOder4Accuracy(const std::vector<Condition> &cnds, const DoubleVector& rs, std::vector<DoubleVector> &x, Direction direction UNUSED_PARAM)
{
    unsigned int en = equationsNumber();

    double h = grid().dimension().step();
    unsigned int N = grid().dimension().sizeN();

    if (en == 1)
    {
        double *p = (double*) malloc(sizeof(double)*(N-2));
        double *q = (double*) malloc(sizeof(double)*(N-2));
        double *v = (double*) malloc(sizeof(double)*(N-2));
        double *u = (double*) malloc(sizeof(double)*(N-2));
        double *r = (double*) malloc(sizeof(double)*(N-2));

        std::vector<unsigned int> ind;
        std::vector<double*> ems;

        /**********************************************************************
         *                  Discretization using trapesium rule
         *********************************************************************/
        double *b = (double*) malloc(sizeof(double)*(N+1));
        for (unsigned int m=0; m<=N; m++) b[m] = 0.0;
        discretisation(cnds, b);

        p[0] = b[0];
        q[0] = b[1];
        v[0] = b[2];
        u[0] = b[3];
        for (unsigned int n=4; n<=N; n++)
        {
            if (fabs(b[n]) > DBL_EPSILON)
            {
                ind.push_back(n);
                ems.push_back((double*) malloc(sizeof(double)*(N-2)));
                ems[ems.size()-1][0] = b[n];
            }
        }
        free(b);
        r[0] = rs[0];

        unsigned int ind_size = ind.size();
        unsigned int ems_size = ems.size();

        /**********************************************************************
         *                          End of discretization
         *********************************************************************/

        /**********************************************************************
         *                      Finding function at end of grid
         *********************************************************************/


        for (unsigned int n=0; n<=N-4; n++)
        {
            double t = n*h;
            double m = +1.0/(-12.0*h*A(t,n)-25.0);
            double alpha1 = -48.0*m;
            double alpha2 = +36.0*m;
            double alpha3 = -16.0*m;
            double alpha4 = +3.0*m;
            double alpha0 = +12.0*h*B(t,n)*m;

            r[n+1] = r[n] - p[n]*alpha0;
            u[n+1] = p[n]*alpha4;
            v[n+1] = p[n]*alpha3 + u[n];
            q[n+1] = p[n]*alpha2 + v[n];
            p[n+1] = p[n]*alpha1 + q[n];

            for (unsigned int i=0; i<ind_size; i++)
            {
                if (n+4 == ind.at(i))
                {
                    u[n+1] += ems[i][n];
                }
                if (n+4 >= ind.at(i))
                {
                    ems[i][n+1] = 0.0;
                }
                else
                {
                    ems[i][n+1] = ems[i][n];
                }
            }
        }

        DoubleMatrix m(5,5);
        DoubleVector c(5);
        DoubleVector xT(5);

        m[0][0] = 0.0;
        m[0][1] = p[N-3];
        m[0][2] = q[N-3];
        m[0][3] = v[N-3];
        m[0][4] = u[N-3];
        c[0] = r[N-3];

        m[1][0] = -3.0;
        m[1][1] = -10.0 - 12.0*h*A((N-3)*h,(N-3));
        m[1][2] = +18.0;
        m[1][3] = -6.0;
        m[1][4] = +1.0;
        c[1] = 12.0*h*B((N-3)*h,N-3);

        m[2][0] = +1.0;
        m[2][1] = -8.0;
        m[2][2] = +0.0 -12.0*h*A((N-2)*h,N-2);
        m[2][3] = +8.0;
        m[2][4] = -1.0;
        c[2] = 12.0*h*B((N-2)*h,N-2);

        m[3][0] = -1.0;
        m[3][1] = +6.0;
        m[3][2] = -18.0;
        m[3][3] = +10.0 - 12.0*h*A((N-1)*h,N-1);
        m[3][4] = +3.0;
        c[3] = 12.0*h*B((N-1)*h,N-1);

        m[4][0] = +3.0;
        m[4][1] = -16.0;
        m[4][2] = +36.0;
        m[4][3] = -48.0;
        m[4][4] = +25.0 - 12.0*h*A((N-0)*h,N-0);
        c[4] = 12.0*h*B((N-0)*h,N-0);

        LinearEquation::GaussianElimination(m, c, xT);
        m.clear();
        c.clear();

        /**********************************************************************
         *                      Finding function at end of grid
         *********************************************************************/
        DoubleVector x0(N+1);
        x0[N-0] = xT[4];
        x0[N-1] = xT[3];
        x0[N-2] = xT[2];
        x0[N-3] = xT[1];
        x0[N-4] = xT[0];
        xT.clear();

        for (unsigned int i=N-4; i!=0; i--)
        {
            x0[i-1] = -q[i-1]*x0[i]-v[i-1]*x0[i+1]-u[i-1]*x0[i+2]+r[i-1];
            for (unsigned int s=0; s<ems_size; s++) x0[i-1] -= ems[s][i-1]*x0[ind[s]];
            x0[i-1] /= p[i-1];
        }

        x.clear();
        x.push_back(x0);

        free(p);
        free(q);
        free(v);
        free(u);
        free(r);
        for (unsigned int i=0; i<ems_size; i++) free(ems[i]);
        ems.clear();
        ind.clear();

    }
    else
    {
        DoubleMatrix* p = new DoubleMatrix[N-2]; for (unsigned int i=0; i<N-2; i++) p[i].resize(en, en, 0.0);
        DoubleMatrix* q = new DoubleMatrix[N-2]; for (unsigned int i=0; i<N-2; i++) q[i].resize(en, en, 0.0);
        DoubleMatrix* v = new DoubleMatrix[N-2]; for (unsigned int i=0; i<N-2; i++) v[i].resize(en, en, 0.0);
        DoubleMatrix* u = new DoubleMatrix[N-2]; for (unsigned int i=0; i<N-2; i++) u[i].resize(en, en, 0.0);
        DoubleMatrix* r = new DoubleMatrix[N-2]; for (unsigned int i=0; i<N-2; i++) r[i].resize(en, 1, 0.0);

        std::vector<unsigned int> ind;
        std::vector<DoubleMatrix*> ems;

        /**********************************************************************
         *                  Discretization using trapesium rule
         *********************************************************************/
        DoubleMatrix* b = new DoubleMatrix[N+1];
        for (unsigned int n=0; n<=N; n++) b[n].resize(en, en, 0.0);

        unsigned int cnd_size = cnds.size();
        for (unsigned int s=0; s<cnd_size; s++)
        {
            const Condition &cnd = cnds[s];
            double time  = cnd.time;

            for (unsigned int n=0; n<=N; n++)
            {
                double dh = fabs(time - n*h);
                if (dh <= h) b[n] += (1.0 - dh/h)*cnd.mtrx;
            }
        }

        p[0] = b[0];
        q[0] = b[1];
        v[0] = b[2];
        u[0] = b[3];
        for (unsigned int n=4; n<=N; n++)
        {
            if (!b[n].zeroMatrix())
            {
                ind.push_back(n);
                DoubleMatrix *dm = new DoubleMatrix[N];
                for (unsigned int i=0; i<N; i++) dm[i].resize(en,en);
                ems.push_back(dm);
                ems[ems.size()-1][0] = b[n];
            }
        }
        free(b);
        r[0] = rs;

        unsigned int ind_size = ind.size();
        unsigned int ems_size = ems.size();

        //printf("%d %d %d %d %d %d %d %d\n", ind[0], ind[1], ind[2], ind[3], ind[4], ind[5], ind[6], ind[7]);

        /**********************************************************************
         *                          End of discretization
         *********************************************************************/

        /**********************************************************************
         *                      Finding function at end of grid
         *********************************************************************/

        DoubleMatrix m(en, en);
        DoubleMatrix alpha1(en, en, 0.0);
        DoubleMatrix alpha2(en, en, 0.0);
        DoubleMatrix alpha3(en, en, 0.0);
        DoubleMatrix alpha4(en, en, 0.0);
        DoubleMatrix alpha0(en, 1, 0.0);

        for (unsigned int n=0; n<=N-4; n++)
        {
            double t = n*h;

            for (unsigned int row=0; row<en; row++)
            {
                for (unsigned int col=0; col<en; col++)
                {
                    m[row][col] = -12.0*h*A(t,n,row,col);
                    if (row==col) m[row][col] -= 25.0;

                    alpha1[row][col] = 0.0;
                    if (row==col) alpha1[row][col] = -48.0;

                    alpha2[row][col] = 0.0;
                    if (row==col) alpha2[row][col] = +36.0;

                    alpha3[row][col] = 0.0;
                    if (row==col) alpha3[row][col] = -16.0;

                    alpha4[row][col] = 0.0;
                    if (row==col) alpha4[row][col] = +3.0;

                }
                alpha0[row][0] = +12.0*h*B(t,n,row);
            }

            m.inverse();

            alpha1 = m*alpha1;
            alpha2 = m*alpha2;
            alpha3 = m*alpha3;
            alpha4 = m*alpha4;
            alpha0 = m*alpha0;

            r[n+1] = r[n] - p[n]*alpha0;
            u[n+1] = p[n]*alpha4;
            v[n+1] = p[n]*alpha3 + u[n];
            q[n+1] = p[n]*alpha2 + v[n];
            p[n+1] = p[n]*alpha1 + q[n];

            for (unsigned int i=0; i<ind_size; i++)
            {
                if (n+4 == ind.at(i))
                {
                    u[n+1] += ems[i][n];
                    //ems[i][n] = DoubleMatrix::ZeroMatrix(en, en);
                }
                if (n+4 >= ind.at(i))
                {
                    ems[i][n+1] = DoubleMatrix::ZeroMatrix(en, en);
                }
                else
                {
                    ems[i][n+1] = ems[i][n];
                }
            }
        }

        m.clear();
        alpha0.clear();
        alpha1.clear();
        alpha2.clear();
        alpha3.clear();
        alpha4.clear();

        DoubleMatrix M(5*en,5*en);
        DoubleVector C(5*en);
        DoubleVector xT(5*en);

        for (unsigned int row=0*en; row<1*en; row++)
        {
            for (unsigned int col=0*en; col<1*en; col++) M[row][col] = 0.0;
            for (unsigned int col=1*en; col<2*en; col++) M[row][col] = p[N-3][row%en][col%en];
            for (unsigned int col=2*en; col<3*en; col++) M[row][col] = q[N-3][row%en][col%en];
            for (unsigned int col=3*en; col<4*en; col++) M[row][col] = v[N-3][row%en][col%en];
            for (unsigned int col=4*en; col<5*en; col++) M[row][col] = u[N-3][row%en][col%en];
            C[row] = r[N-3][row%en][0];
        }

        for (unsigned int row=1*en; row<2*en; row++)
        {
            for (unsigned int col=0*en; col<1*en; col++) { M[row][col] = +0.0; if (row%en == col%en) M[row][col] += -3.0; }
            for (unsigned int col=1*en; col<2*en; col++) { M[row][col] = -12.0*h*A((N-3)*h, N-3, row%en, col%en); if (row%en == col%en) M[row][col] += -10.0; }
            for (unsigned int col=2*en; col<3*en; col++) { M[row][col] = +0.0; if (row%en == col%en) M[row][col] += +18.0; }
            for (unsigned int col=3*en; col<4*en; col++) { M[row][col] = +0.0; if (row%en == col%en) M[row][col] += -6.0; }
            for (unsigned int col=4*en; col<5*en; col++) { M[row][col] = +0.0; if (row%en == col%en) M[row][col] += +1.0; }
            C[row] = 12.0*h*B((N-3)*h, N-3, row%en);
        }

        for (unsigned int row=2*en; row<3*en; row++)
        {
            for (unsigned int col=0*en; col<1*en; col++) { M[row][col] = +0.0; if (row%en == col%en) M[row][col] += +1.0; }
            for (unsigned int col=1*en; col<2*en; col++) { M[row][col] = +0.0; if (row%en == col%en) M[row][col] += -8.0; }
            for (unsigned int col=2*en; col<3*en; col++) { M[row][col] = -12.0*h*A((N-2)*h, N-2, row%en, col%en); }
            for (unsigned int col=3*en; col<4*en; col++) { M[row][col] = +0.0; if (row%en == col%en) M[row][col] += +8.0; }
            for (unsigned int col=4*en; col<5*en; col++) { M[row][col] = +0.0; if (row%en == col%en) M[row][col] += -1.0; }
            C[row] = 12.0*h*B((N-2)*h, N-2, row%en);
        }

        for (unsigned int row=3*en; row<4*en; row++)
        {
            for (unsigned int col=0*en; col<1*en; col++) { M[row][col] = +0.0; if (row%en == col%en) M[row][col] += -1.0; }
            for (unsigned int col=1*en; col<2*en; col++) { M[row][col] = +0.0; if (row%en == col%en) M[row][col] += +6.0; }
            for (unsigned int col=2*en; col<3*en; col++) { M[row][col] = +0.0; if (row%en == col%en) M[row][col] += -18.0; }
            for (unsigned int col=3*en; col<4*en; col++) { M[row][col] = -12.0*h*A((N-1)*h, N-1, row%en, col%en); if (row%en == col%en) M[row][col] += +10.0; }
            for (unsigned int col=4*en; col<5*en; col++) { M[row][col] = +0.0; if (row%en == col%en) M[row][col] += +3.0; }
            C[row] = 12.0*h*B((N-1)*h, N-1, row%en);
        }

        for (unsigned int row=4*en; row<5*en; row++)
        {
            for (unsigned int col=0*en; col<1*en; col++) { M[row][col] = +0.0; if (row%en == col%en) M[row][col] += +3.0; }
            for (unsigned int col=1*en; col<2*en; col++) { M[row][col] = +0.0; if (row%en == col%en) M[row][col] += -16.0; }
            for (unsigned int col=2*en; col<3*en; col++) { M[row][col] = +0.0; if (row%en == col%en) M[row][col] += +36.0; }
            for (unsigned int col=3*en; col<4*en; col++) { M[row][col] = +0.0; if (row%en == col%en) M[row][col] += -48.0; }
            for (unsigned int col=4*en; col<5*en; col++) { M[row][col] = -12.0*h*A((N-0)*h, N-0, row%en, col%en); if (row%en == col%en) M[row][col] += +25.0;}
            C[row] = 12.0*h*B((N-0)*h, N-0, row%en);
        }

        //        unsigned int row = 12;
        //        double tN4 = h*(N-4);
        //        double tN3 = h*(N-3);
        //        double tN2 = h*(N-2);
        //        double tN1 = h*(N-1);
        //        double tN0 = h*(N-0);
        //        printf("%.10f\n%.10f\n", C[row], M[row][0]*X(tN4,0)+M[row][1]*X(tN4,1)+M[row][2]*X(tN4,2)
        //                                     +M[row][3]*X(tN3,0)+M[row][4]*X(tN3,1)+M[row][5]*X(tN3,2)
        //                                     +M[row][6]*X(tN2,0)+M[row][7]*X(tN2,1)+M[row][8]*X(tN2,2)
        //                                     +M[row][9]*X(tN1,0)+M[row][10]*X(tN1,1)+M[row][11]*X(tN1,2)
        //                                     +M[row][12]*X(tN0,0)+M[row][13]*X(tN0,1)+M[row][14]*X(tN0,2));

        LinearEquation::GaussianElimination(M, C, xT);
        //        IPrinter::print(xT, xT.size());

        x.resize(en);
        for (unsigned int i=0; i<en; i++) x[i].resize(N+1);

        std::vector<DoubleMatrix> x0(N+1);
        for (unsigned int i=0; i<=N; i++) x0[i].resize(en,1);

        for (unsigned int row=0; row<en; row++)
        {
            x0[N-4][row][0] = xT[row+0*en];
            x0[N-3][row][0] = xT[row+1*en];
            x0[N-2][row][0] = xT[row+2*en];
            x0[N-1][row][0] = xT[row+3*en];
            x0[N-0][row][0] = xT[row+4*en];
        }

        for (unsigned int i=N-4; i!=0; i--)
        {
            x0[i-1] = r[i-1] - q[i-1]*x0[i] - v[i-1]*x0[i+1] - u[i-1]*x0[i+2];
            for (unsigned int s=0; s<ems_size; s++) x0[i-1] += -1.0*(ems[s][i-1]*x0[ind[s]]);
            p[i-1].inverse();
            x0[i-1] = p[i-1] * x0[i-1];
        }

        for (unsigned int i=0; i<=N; i++)
        {
            for (unsigned int row=0; row<en; row++)
            {
                x[row][i] = x0[i][row][0];
            }
        }

        for (unsigned int i=0; i<N-2; i++)
        {
            r[i].clear();
            u[i].clear();
            v[i].clear();
            q[i].clear();
            p[i].clear();
        }
        delete [] u;
        delete [] v;
        delete [] r;
        delete [] q;
        delete [] p;
        ind.clear();

        for (unsigned int i=0; i<ems_size; i++)
        {
            delete [] ems[i];
        }
        ems.clear();

        for (unsigned int i=0; i<=N; i++)
        {
            x0[i].clear();
        }
        x0.clear();

        M.clear();
        C.clear();
        xT.clear();
    }
}

void LinearODE1stOrder::highOder6Accuracy(const std::vector<Condition> &cnds, const DoubleVector& rs, std::vector<DoubleVector> &x, Direction direction UNUSED_PARAM)
{
    unsigned int en = equationsNumber();

    double h = grid().dimension().step();
    unsigned int N = grid().dimension().sizeN();

    if (en == 1)
    {
        double *p = (double*) malloc(sizeof(double)*(N-4));
        double *q = (double*) malloc(sizeof(double)*(N-4));
        double *v = (double*) malloc(sizeof(double)*(N-4));
        double *u = (double*) malloc(sizeof(double)*(N-4));
        double *w = (double*) malloc(sizeof(double)*(N-4));
        double *z = (double*) malloc(sizeof(double)*(N-4));
        double *r = (double*) malloc(sizeof(double)*(N-4));

        std::vector<unsigned int> ind;
        std::vector<double*> ems;

        /**********************************************************************
         *                  Discretization using trapesium rule
         *********************************************************************/
        double *b = (double*) malloc(sizeof(double)*(N+1));
        for (unsigned int m=0; m<=N; m++) b[m] = 0.0;
        discretisation(cnds, b);

        p[0] = b[0];
        q[0] = b[1];
        v[0] = b[2];
        u[0] = b[3];
        w[0] = b[4];
        z[0] = b[5];
        for (unsigned int n=6; n<=N; n++)
        {
            if (fabs(b[n]) > DBL_EPSILON)
            {
                ind.push_back(n);
                ems.push_back((double*) malloc(sizeof(double)*(N-4)));
                ems[ems.size()-1][0] = b[n];
            }
        }
        free(b);
        r[0] = rs[0];

        unsigned int ind_size = ind.size();
        unsigned int ems_size = ems.size();

        /**********************************************************************
         *                          End of discretization
         *********************************************************************/

        /**********************************************************************
         *                      Finding function at end of grid
         *********************************************************************/


        for (unsigned int n=0; n<=N-6; n++)
        {
            double t = n*h;
            double m = +1.0/(-60.0*h*A(t,n)-147.0);
            double alpha1 = -360.0*m;
            double alpha2 = +450.0*m;
            double alpha3 = -400.0*m;
            double alpha4 = +225.0*m;
            double alpha5 = -72.0*m;
            double alpha6 = +10.0*m;
            double alpha0 = +60.0*h*B(t,n)*m;

            r[n+1] = r[n] - p[n]*alpha0;
            z[n+1] = p[n]*alpha6;
            w[n+1] = p[n]*alpha5 + z[n];
            u[n+1] = p[n]*alpha4 + w[n];
            v[n+1] = p[n]*alpha3 + u[n];
            q[n+1] = p[n]*alpha2 + v[n];
            p[n+1] = p[n]*alpha1 + q[n];

            for (unsigned int i=0; i<ind_size; i++)
            {
                if (n+6 == ind.at(i))
                {
                    z[n+1] += ems[i][n];
                }
                if (n+6 >= ind.at(i))
                {
                    ems[i][n+1] = 0.0;
                }
                else
                {
                    ems[i][n+1] = ems[i][n];
                }
            }
        }

        DoubleMatrix m(7,7);
        DoubleVector c(7);
        DoubleVector xT(7);

        m[0][0] = 0.0;
        m[0][1] = p[N-5];
        m[0][2] = q[N-5];
        m[0][3] = v[N-5];
        m[0][4] = u[N-5];
        m[0][5] = w[N-5];
        m[0][6] = z[N-5];
        c[0] = r[N-5];

        m[1][0] = -10.0;
        m[1][1] = -77.0 - 60.0*h*A((N-5)*h,N-5);
        m[1][2] = +150.0;
        m[1][3] = -100.0;
        m[1][4] = +50.0;
        m[1][5] = -15.0;
        m[1][6] = +2.0;
        c[1] = 60.0*h*B((N-5)*h,N-5);

        m[2][0] = +2.0;
        m[2][1] = -24.0;
        m[2][2] = -35.0 - 60.0*h*A((N-4)*h,N-4);
        m[2][3] = +80.0;
        m[2][4] = -30.0;
        m[2][5] = +8.0;
        m[2][6] = -1.0;
        c[2] = 60.0*h*B((N-4)*h,N-4);

        m[3][0] = -1.0;
        m[3][1] = +9.0;
        m[3][2] = -45.0;
        m[3][3] = -60.0*h*A((N-3)*h,N-3);
        m[3][4] = +45.0;
        m[3][5] = -9.0;
        m[3][6] = +1.0;
        c[3] = 60.0*h*B((N-3)*h,N-3);

        m[4][0] = +1.0;
        m[4][1] = -8.0;
        m[4][2] = +30.0;
        m[4][3] = -80.0;
        m[4][4] = +35.0 - 60.0*h*A((N-2)*h,N-2);
        m[4][5] = +24.0;
        m[4][6] = -2.0;
        c[4] = 60.0*h*B((N-2)*h,N-2);

        m[5][0] = -2.0;
        m[5][1] = +15.0;
        m[5][2] = -50.0;
        m[5][3] = +100.0;
        m[5][4] = -150.0;
        m[5][5] = +77.0 - 60.0*h*A((N-1)*h,N-1);
        m[5][6] = +10.0;
        c[5] = 60.0*h*B((N-1)*h,N-1);

        m[6][0] = +10.0;
        m[6][1] = -72.0;
        m[6][2] = +225.0;
        m[6][3] = -400.0;
        m[6][4] = +450.0;
        m[6][5] = -360.0;
        m[6][6] = +147.0 - 60.0*h*A((N-0)*h,N-0);
        c[6] = 60.0*h*B((N-0)*h,N-0);

        LinearEquation::GaussianElimination(m, c, xT);
        m.clear();
        c.clear();

        /**********************************************************************
         *                      Finding function at end of grid
         *********************************************************************/
        DoubleVector x0(N+1);
        x0[N-0] = xT[6];
        x0[N-1] = xT[5];
        x0[N-2] = xT[4];
        x0[N-3] = xT[3];
        x0[N-4] = xT[2];
        x0[N-5] = xT[1];
        x0[N-6] = xT[0];
        xT.clear();

        for (unsigned int i=N-6; i!=0; i--)
        {
            x0[i-1] = -q[i-1]*x0[i]-v[i-1]*x0[i+1]-u[i-1]*x0[i+2]-w[i-1]*x0[i+3]-z[i-1]*x0[i+4]+r[i-1];
            for (unsigned int s=0; s<ems_size; s++) x0[i-1] -= ems[s][i-1]*x0[ind[s]];
            x0[i-1] /= p[i-1];
        }

        x.clear();
        x.push_back(x0);

        free(p);
        free(q);
        free(v);
        free(u);
        free(r);
        free(w);
        free(z);
        for (unsigned int i=0; i<ems_size; i++) free(ems[i]);
        ems.clear();
        ind.clear();

    }
    else
    {
        DoubleMatrix* p = new DoubleMatrix[N-4]; for (unsigned int i=0; i<N-4; i++) p[i].resize(en, en, 0.0);
        DoubleMatrix* q = new DoubleMatrix[N-4]; for (unsigned int i=0; i<N-4; i++) q[i].resize(en, en, 0.0);
        DoubleMatrix* v = new DoubleMatrix[N-4]; for (unsigned int i=0; i<N-4; i++) v[i].resize(en, en, 0.0);
        DoubleMatrix* u = new DoubleMatrix[N-4]; for (unsigned int i=0; i<N-4; i++) u[i].resize(en, en, 0.0);
        DoubleMatrix* w = new DoubleMatrix[N-4]; for (unsigned int i=0; i<N-4; i++) w[i].resize(en, en, 0.0);
        DoubleMatrix* z = new DoubleMatrix[N-4]; for (unsigned int i=0; i<N-4; i++) z[i].resize(en, en, 0.0);
        DoubleMatrix* r = new DoubleMatrix[N-4]; for (unsigned int i=0; i<N-4; i++) r[i].resize(en, 1, 0.0);

        std::vector<unsigned int> ind;
        std::vector<DoubleMatrix*> ems;

        /**********************************************************************
         *                  Discretization using trapesium rule
         *********************************************************************/
        DoubleMatrix* b = new DoubleMatrix[N+1];
        for (unsigned int n=0; n<=N; n++) b[n].resize(en, en, 0.0);

        unsigned int cnd_size = cnds.size();
        for (unsigned int s=0; s<cnd_size; s++)
        {
            const Condition &cnd = cnds[s];
            double time  = cnd.time;

            for (unsigned int n=0; n<=N; n++)
            {
                double dh = fabs(time - n*h);
                if (dh <= h) b[n] += (1.0 - dh/h)*cnd.mtrx;
            }
        }

        p[0] = b[0];
        q[0] = b[1];
        v[0] = b[2];
        u[0] = b[3];
        w[0] = b[4];
        z[0] = b[5];
        for (unsigned int n=6; n<=N; n++)
        {
            if (!b[n].zeroMatrix())
            {
                ind.push_back(n);
                DoubleMatrix *dm = new DoubleMatrix[N];
                for (unsigned int i=0; i<N; i++) dm[i].resize(en,en);
                ems.push_back(dm);
                ems[ems.size()-1][0] = b[n];
            }
        }
        free(b);
        r[0] = rs;

        unsigned int ind_size = ind.size();
        unsigned int ems_size = ems.size();

        //printf("%d %d %d %d %d %d %d %d\n", ind[0], ind[1], ind[2], ind[3], ind[4], ind[5], ind[6], ind[7]);

        /**********************************************************************
         *                          End of discretization
         *********************************************************************/

        /**********************************************************************
         *                      Finding function at end of grid
         *********************************************************************/

        DoubleMatrix m(en, en);
        DoubleMatrix alpha1(en, en, 0.0);
        DoubleMatrix alpha2(en, en, 0.0);
        DoubleMatrix alpha3(en, en, 0.0);
        DoubleMatrix alpha4(en, en, 0.0);
        DoubleMatrix alpha5(en, en, 0.0);
        DoubleMatrix alpha6(en, en, 0.0);
        DoubleMatrix alpha0(en, 1, 0.0);

        for (unsigned int n=0; n<=N-6; n++)
        {
            double t = n*h;

            for (unsigned int row=0; row<en; row++)
            {
                for (unsigned int col=0; col<en; col++)
                {
                    m[row][col] = -60.0*h*A(t,n,row,col);
                    if (row==col) m[row][col] -= 147.0;

                    alpha1[row][col] = 0.0;
                    if (row==col) alpha1[row][col] = -360.0;

                    alpha2[row][col] = 0.0;
                    if (row==col) alpha2[row][col] = +450.0;

                    alpha3[row][col] = 0.0;
                    if (row==col) alpha3[row][col] = -400.0;

                    alpha4[row][col] = 0.0;
                    if (row==col) alpha4[row][col] = +225.0;

                    alpha5[row][col] = 0.0;
                    if (row==col) alpha5[row][col] = -72.0;

                    alpha6[row][col] = 0.0;
                    if (row==col) alpha6[row][col] = +10.0;
                }
                alpha0[row][0] = +60.0*h*B(t,n,row);
            }

            m.inverse();

            alpha1 = m*alpha1;
            alpha2 = m*alpha2;
            alpha3 = m*alpha3;
            alpha4 = m*alpha4;
            alpha5 = m*alpha5;
            alpha6 = m*alpha6;
            alpha0 = m*alpha0;

            r[n+1] = r[n] - p[n]*alpha0;
            z[n+1] = p[n]*alpha6;
            w[n+1] = p[n]*alpha5 + z[n];
            u[n+1] = p[n]*alpha4 + w[n];
            v[n+1] = p[n]*alpha3 + u[n];
            q[n+1] = p[n]*alpha2 + v[n];
            p[n+1] = p[n]*alpha1 + q[n];

            for (unsigned int i=0; i<ind_size; i++)
            {
                if (n+6 == ind.at(i))
                {
                    z[n+1] += ems[i][n];
                    //ems[i][n] = DoubleMatrix::ZeroMatrix(en, en);
                }
                if (n+6 >= ind.at(i))
                {
                    ems[i][n+1] = DoubleMatrix::ZeroMatrix(en, en);
                }
                else
                {
                    ems[i][n+1] = ems[i][n];
                }
            }
        }

        m.clear();
        alpha0.clear();
        alpha1.clear();
        alpha2.clear();
        alpha3.clear();
        alpha4.clear();
        alpha5.clear();
        alpha6.clear();

        DoubleMatrix M(7*en,7*en);
        DoubleVector C(7*en);
        DoubleVector xT(7*en);

        for (unsigned int row=0*en; row<1*en; row++)
        {
            for (unsigned int col=0*en; col<1*en; col++) M[row][col] = 0.0;
            for (unsigned int col=1*en; col<2*en; col++) M[row][col] = p[N-5][row%en][col%en];
            for (unsigned int col=2*en; col<3*en; col++) M[row][col] = q[N-5][row%en][col%en];
            for (unsigned int col=3*en; col<4*en; col++) M[row][col] = v[N-5][row%en][col%en];
            for (unsigned int col=4*en; col<5*en; col++) M[row][col] = u[N-5][row%en][col%en];
            for (unsigned int col=5*en; col<6*en; col++) M[row][col] = w[N-5][row%en][col%en];
            for (unsigned int col=6*en; col<7*en; col++) M[row][col] = z[N-5][row%en][col%en];
            C[row] = r[N-5][row%en][0];
        }

        for (unsigned int row=1*en; row<2*en; row++)
        {
            for (unsigned int col=0*en; col<1*en; col++) { M[row][col] = +0.0; if (row%en == col%en) M[row][col] += -10.0; }
            for (unsigned int col=1*en; col<2*en; col++) { M[row][col] = -60.0*h*A((N-5)*h, N-5, row%en, col%en); if (row%en == col%en) M[row][col] += -77.0; }
            for (unsigned int col=2*en; col<3*en; col++) { M[row][col] = +0.0; if (row%en == col%en) M[row][col] += +150.0; }
            for (unsigned int col=3*en; col<4*en; col++) { M[row][col] = +0.0; if (row%en == col%en) M[row][col] += -100.0; }
            for (unsigned int col=4*en; col<5*en; col++) { M[row][col] = +0.0; if (row%en == col%en) M[row][col] += +50.0; }
            for (unsigned int col=5*en; col<6*en; col++) { M[row][col] = +0.0; if (row%en == col%en) M[row][col] += -15.0; }
            for (unsigned int col=6*en; col<7*en; col++) { M[row][col] = +0.0; if (row%en == col%en) M[row][col] += +2.0; }
            C[row] = 60.0*h*B((N-5)*h, N-3, row%en);
        }

        for (unsigned int row=2*en; row<3*en; row++)
        {
            for (unsigned int col=0*en; col<1*en; col++) { M[row][col] = +0.0; if (row%en == col%en) M[row][col] += +2.0; }
            for (unsigned int col=1*en; col<2*en; col++) { M[row][col] = +0.0; if (row%en == col%en) M[row][col] += -24.0; }
            for (unsigned int col=2*en; col<3*en; col++) { M[row][col] = -60.0*h*A((N-4)*h, N-4, row%en, col%en); if (row%en == col%en) M[row][col] += -35.0; }
            for (unsigned int col=3*en; col<4*en; col++) { M[row][col] = +0.0; if (row%en == col%en) M[row][col] += +80.0; }
            for (unsigned int col=4*en; col<5*en; col++) { M[row][col] = +0.0; if (row%en == col%en) M[row][col] += -30.0; }
            for (unsigned int col=5*en; col<6*en; col++) { M[row][col] = +0.0; if (row%en == col%en) M[row][col] += +8.0; }
            for (unsigned int col=6*en; col<7*en; col++) { M[row][col] = +0.0; if (row%en == col%en) M[row][col] += -1.0; }
            C[row] = 60.0*h*B((N-4)*h, N-4, row%en);
        }

        for (unsigned int row=3*en; row<4*en; row++)
        {
            for (unsigned int col=0*en; col<1*en; col++) { M[row][col] = +0.0; if (row%en == col%en) M[row][col] += -1.0; }
            for (unsigned int col=1*en; col<2*en; col++) { M[row][col] = +0.0; if (row%en == col%en) M[row][col] += +9.0; }
            for (unsigned int col=2*en; col<3*en; col++) { M[row][col] = +0.0; if (row%en == col%en) M[row][col] += -45.0; }
            for (unsigned int col=3*en; col<4*en; col++) { M[row][col] = -60.0*h*A((N-3)*h, N-3, row%en, col%en);}
            for (unsigned int col=4*en; col<5*en; col++) { M[row][col] = +0.0; if (row%en == col%en) M[row][col] += +45.0; }
            for (unsigned int col=5*en; col<6*en; col++) { M[row][col] = +0.0; if (row%en == col%en) M[row][col] += -9.0; }
            for (unsigned int col=6*en; col<7*en; col++) { M[row][col] = +0.0; if (row%en == col%en) M[row][col] += +1.0; }
            C[row] = 60.0*h*B((N-3)*h, N-3, row%en);
        }

        for (unsigned int row=4*en; row<5*en; row++)
        {
            for (unsigned int col=0*en; col<1*en; col++) { M[row][col] = +0.0; if (row%en == col%en) M[row][col] += +1.0; }
            for (unsigned int col=1*en; col<2*en; col++) { M[row][col] = +0.0; if (row%en == col%en) M[row][col] += -8.0; }
            for (unsigned int col=2*en; col<3*en; col++) { M[row][col] = +0.0; if (row%en == col%en) M[row][col] += +30.0; }
            for (unsigned int col=3*en; col<4*en; col++) { M[row][col] = +0.0; if (row%en == col%en) M[row][col] += -80.0; }
            for (unsigned int col=4*en; col<5*en; col++) { M[row][col] = -60.0*h*A((N-2)*h, N-2, row%en, col%en); if (row%en == col%en) M[row][col] += +35.0; }
            for (unsigned int col=5*en; col<6*en; col++) { M[row][col] = +0.0; if (row%en == col%en) M[row][col] += +24.0; }
            for (unsigned int col=6*en; col<7*en; col++) { M[row][col] = +0.0; if (row%en == col%en) M[row][col] += -2.0; }
            C[row] = 60.0*h*B((N-2)*h, N-2, row%en);
        }

        for (unsigned int row=5*en; row<6*en; row++)
        {
            for (unsigned int col=0*en; col<1*en; col++) { M[row][col] = +0.0; if (row%en == col%en) M[row][col] += -2.0; }
            for (unsigned int col=1*en; col<2*en; col++) { M[row][col] = +0.0; if (row%en == col%en) M[row][col] += +15.0; }
            for (unsigned int col=2*en; col<3*en; col++) { M[row][col] = +0.0; if (row%en == col%en) M[row][col] += -50.0; }
            for (unsigned int col=3*en; col<4*en; col++) { M[row][col] = +0.0; if (row%en == col%en) M[row][col] += +100.0; }
            for (unsigned int col=4*en; col<5*en; col++) { M[row][col] = +0.0; if (row%en == col%en) M[row][col] += -150.0; }
            for (unsigned int col=5*en; col<6*en; col++) { M[row][col] = -60.0*h*A((N-1)*h, N-1, row%en, col%en); if (row%en == col%en) M[row][col] += +77.0; }
            for (unsigned int col=6*en; col<7*en; col++) { M[row][col] = +0.0; if (row%en == col%en) M[row][col] += +10.0; }
            C[row] = 60.0*h*B((N-1)*h, N-1, row%en);
        }

        for (unsigned int row=6*en; row<7*en; row++)
        {
            for (unsigned int col=0*en; col<1*en; col++) { M[row][col] = +0.0; if (row%en == col%en) M[row][col] += +10.0; }
            for (unsigned int col=1*en; col<2*en; col++) { M[row][col] = +0.0; if (row%en == col%en) M[row][col] += -72.0; }
            for (unsigned int col=2*en; col<3*en; col++) { M[row][col] = +0.0; if (row%en == col%en) M[row][col] += +225.0; }
            for (unsigned int col=3*en; col<4*en; col++) { M[row][col] = +0.0; if (row%en == col%en) M[row][col] += -400.0; }
            for (unsigned int col=4*en; col<5*en; col++) { M[row][col] = +0.0; if (row%en == col%en) M[row][col] += +450.0; }
            for (unsigned int col=5*en; col<6*en; col++) { M[row][col] = +0.0; if (row%en == col%en) M[row][col] += -360.0; }
            for (unsigned int col=6*en; col<7*en; col++) { M[row][col] = -60.0*h*A((N-0)*h, N-5, row%en, col%en); if (row%en == col%en) M[row][col] += 147.0; }
            C[row] = 60.0*h*B((N-0)*h, N-0, row%en);
        }

        //unsigned int row = 3;

        //        for (unsigned int row = 0; row < 7*en; row++)
        //        {
        //        double tN6 = h*(N-6);
        //        double tN5 = h*(N-5);
        //        double tN4 = h*(N-4);
        //        double tN3 = h*(N-3);
        //        double tN2 = h*(N-2);
        //        double tN1 = h*(N-1);
        //        double tN0 = h*(N-0);
        //        printf("%4d %.10f %.10f\n", row, C[row],
        //                                     +M[row][0]*X(tN6,0) +M[row][1]*X(tN6,1) +M[row][2]*X(tN6,2)
        //                                     +M[row][3]*X(tN5,0) +M[row][4]*X(tN5,1) +M[row][5]*X(tN5,2)
        //                                     +M[row][6]*X(tN4,0) +M[row][7]*X(tN4,1) +M[row][8]*X(tN4,2)
        //                                     +M[row][9]*X(tN3,0) +M[row][10]*X(tN3,1)+M[row][11]*X(tN3,2)
        //                                     +M[row][12]*X(tN2,0)+M[row][13]*X(tN2,1)+M[row][14]*X(tN2,2)
        //                                     +M[row][15]*X(tN1,0)+M[row][16]*X(tN1,1)+M[row][17]*X(tN1,2)
        //                                     +M[row][18]*X(tN0,0)+M[row][19]*X(tN0,1)+M[row][20]*X(tN0,2));
        //        }
        LinearEquation::GaussianElimination(M, C, xT);
        //IPrinter::print(xT, xT.size());

        x.resize(en);
        for (unsigned int i=0; i<en; i++) x[i].resize(N+1);

        std::vector<DoubleMatrix> x0(N+1);
        for (unsigned int i=0; i<=N; i++) x0[i].resize(en,1);

        for (unsigned int row=0; row<en; row++)
        {
            x0[N-6][row][0] = xT[row+0*en];
            x0[N-5][row][0] = xT[row+1*en];
            x0[N-4][row][0] = xT[row+2*en];
            x0[N-3][row][0] = xT[row+3*en];
            x0[N-2][row][0] = xT[row+4*en];
            x0[N-1][row][0] = xT[row+5*en];
            x0[N-0][row][0] = xT[row+6*en];
        }

        for (unsigned int i=N-6; i!=0; i--)
        {
            x0[i-1] = r[i-1] - q[i-1]*x0[i] - v[i-1]*x0[i+1] - u[i-1]*x0[i+2] - w[i-1]*x0[i+3] - z[i-1]*x0[i+4];
            for (unsigned int s=0; s<ems_size; s++) x0[i-1] += -1.0*(ems[s][i-1]*x0[ind[s]]);
            p[i-1].inverse();
            x0[i-1] = p[i-1] * x0[i-1];
        }

        for (unsigned int i=0; i<=N; i++)
        {
            for (unsigned int row=0; row<en; row++)
            {
                x[row][i] = x0[i][row][0];
            }
        }

        for (unsigned int i=0; i<N-4; i++)
        {
            r[i].clear();
            z[i].clear();
            w[i].clear();
            u[i].clear();
            v[i].clear();
            q[i].clear();
            p[i].clear();
        }
        delete [] z;
        delete [] w;
        delete [] u;
        delete [] v;
        delete [] r;
        delete [] q;
        delete [] p;
        ind.clear();

        for (unsigned int i=0; i<ems_size; i++)
        {
            delete [] ems[i];
        }
        ems.clear();

        for (unsigned int i=0; i<=N; i++)
        {
            x0[i].clear();
        }
        x0.clear();

        M.clear();
        C.clear();
        xT.clear();
    }
}


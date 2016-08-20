#include "problem2.h"

void Problem2::Main(int argc UNUSED_PARAM, char *argv[] UNUSED_PARAM)
{
    Problem2 p;

    DoubleVector k0(p.L);
    k0[0] = 1.0;
    k0[1] = 1.0;
    p.calculateV(k0);
//    IPrinter::printVector(p.V);

    DoubleVector k(p.L);
    k[0] = 2.1;
    k[1] = 2.2;

    puts("Analitic");
    DoubleVector g(p.L,0.0);
    p.gradient(k, g);
    printf("%f %f\n", k[0], k[1]);
    printf("%f %f\n", g[0], g[1]);
    g.L2Normalize();
    printf("%f %f\n", g[0], g[1]);

    puts("Numerical");
    double h = 0.01;
    DoubleVector gn(p.L,0.0);
    IGradient::Gradient(&p, h, k, gn);
    printf("%f %f\n", k[0], k[1]);
    printf("%f %f\n", gn[0], gn[1]);
    gn.L2Normalize();
    printf("%f %f\n", gn[0], gn[1]);
}

Problem2::Problem2()
{
    t0 = 0.0; t1 = 1.0;
    x0 = 0.0; x1 = 1.0;
    hx = 0.01;
    ht = 0.01;
    N = 100;
    M = 100;
    a = 1.0;

    lambdaM = 1.0;
    lambdaL = 1.0;
    lambdaR = 1.0;

    L = 2;

    xi.resize(L);
    xi[0] = 0.4;
    xi[1] = 0.7;

    Xi.resize(L);
    Xi[0] = 40;
    Xi[1] = 70;

    z.resize(L);
    z[0] = 2.89;
    z[1] = 2.81;

    alpha1 = 1.0;
    alpha2 = 0.0;

    Te = 3.0;
}

Problem2::~Problem2()
{}

double Problem2::fx(const DoubleVector &k)
{
    pk = &k;
    DoubleMatrix u;
    calculateU(u, ht, hx, M, N, lambdaM, lambdaL, lambdaR, a);
    pu = &u;
    //IPrinter::printVector(u.row(M));

    double sum = 0.0;
    for (unsigned int n=0; n<=N-1; n++)
    {
        unsigned int n1 = n;
        unsigned int n2 = n + 1;
        double f1 = mu(n1) * (u.at(M, n1) - V[n1])*(u.at(M, n1) - V[n1]);
        double f2 = mu(n2) * (u.at(M, n2) - V[n2])*(u.at(M, n2) - V[n2]);
        sum = sum + (f1 + f2);
    }
    sum = 0.5*hx*sum;

    double norm = 0.0;
    for (unsigned int m=0; m<=M-1; m++)
    {
        unsigned int m1 = m + 0;
        unsigned int m2 = m + 1;
        //norm += (vl(m1) + vl(m2));
        norm = norm + (k[0]*(u.at(m1,Xi[0])-z[0]) + k[1]*(u.at(m1,Xi[1])-z[1]))
                    + (k[0]*(u.at(m2,Xi[0])-z[0]) + k[1]*(u.at(m2,Xi[1])-z[1]));
    }
    norm = 0.5*ht * norm;

    return alpha1*sum + alpha2*norm;
}

void Problem2::gradient(const DoubleVector &k UNUSED_PARAM, DoubleVector &g UNUSED_PARAM)
{
    pk = &k;
    //puts("---");
    DoubleMatrix u;
    calculateU(u, ht, hx, M, N, lambdaM, lambdaL, lambdaR, a);
    pu = &u;
    //IPrinter::printMatrix(u);

    //puts("---");
    DoubleMatrix psi;
    calculateP(psi, ht, hx, M, N, lambdaM, lambdaL, lambdaR, a);
    pp = &psi;
    //IPrinter::printMatrix(psi);

    g[0] = g[1] = 0.0;
    for (unsigned int s = 0; s<L; s++)
    {
        double sum = 0.0;
        for (unsigned int m=0; m<=M-1; m++)
        {
            unsigned int m1 = m + 0;
            unsigned int m2 = m + 1;
            double g1 = (u.at(m1, Xi[s]) - z[s]) * (-lambdaM*a*a*psi.at(m1, 0) + 2.0*alpha2*(k[0]*(u.at(m1,Xi[0]) - z[0]) + k[1]*(u.at(m1,Xi[1]) - z[1])));
            double g2 = (u.at(m2, Xi[s]) - z[s]) * (-lambdaM*a*a*psi.at(m2, 0) + 2.0*alpha2*(k[0]*(u.at(m2,Xi[0]) - z[0]) + k[1]*(u.at(m2,Xi[1]) - z[1])));
            sum = sum + (g1 + g2);
        }
        sum = 0.5 * ht * sum;
        g[s] = sum;
    }
}

double Problem2::initial(unsigned int i UNUSED_PARAM) const { return 1.0; }

double Problem2::vm(unsigned int j UNUSED_PARAM) const
{
    return Te;
}

double Problem2::vl(unsigned int j UNUSED_PARAM) const
{
    const DoubleVector &k = *pk;
    const DoubleMatrix &u = *pu;
    return k[0] * (u.at(j, Xi[0]) - z[0]) + k[1] * (u.at(j, Xi[1]) - z[1]);
}

double Problem2::vr(unsigned int j UNUSED_PARAM) const
{
    return Te;
}

void Problem2::calculateU(DoubleMatrix &u, double ht, double hx, unsigned int M, unsigned int N, double lambdaM, double lambdaL, double lambdaR, double a)
{
    const DoubleVector &k = *pk;

    u.clear();
    u.resize(M+1, N+1);

    DoubleVector rb(N+1);
    DoubleVector rx(N+1);

    for (unsigned int j=0; j<=M; j++)
    {
        if (j==0)
        {
            for (unsigned int i=0; i<=N; i++) u.at(0,i) = initial(i);
        }
        else
        {
            DoubleMatrix ra(N+1, N+1, 0.0);
            ra(0,0) = 1.0 + (a*a*ht)/(hx*hx) + (lambdaL*a*a*ht)/hx + lambdaM*ht;
            ra(0,1) = -(a*a*ht)/(hx*hx);
            rb[0] = u.at(j-1,0) + ((lambdaL*a*a*ht)/(hx))*(-(k[0]*z[0]-k[1]*z[1])) + lambdaM*ht*vm(j);

            ra(0, Xi[0]) = -k[0] * ((lambdaL*a*a*ht)/(hx));
            ra(0, Xi[1]) = -k[1] * ((lambdaL*a*a*ht)/(hx));

            for (unsigned int i=1; i<=N-1; i++)
            {
                ra(i,i-1) = -(a*a*ht)/(hx*hx);
                ra(i,i) = 1.0 + 2.0*((a*a)*ht)/(hx*hx) + lambdaM*ht;
                ra(i,i+1) = -(a*a*ht)/(hx*hx);
                rb[i] = u.at(j-1, i) + lambdaM*ht*vm(j);
            }

            ra(N,N-1) = -(a*a*ht)/(hx*hx);
            ra(N,N)   = 1.0 + (a*a*ht)/(hx*hx) + (lambdaL*a*a*ht)/hx + lambdaM*ht;
            rb[N] = u.at(j-1,N) + ((lambdaR*a*a*ht)/(hx))*vr(j) + lambdaM*ht*vm(j);

            GaussianElimination(ra, rb, rx);

            for (unsigned int i=0; i<=N; i++) u.at(j,i) = rx[i];
        }
    }
}

void Problem2::calculateP(DoubleMatrix &psi, double ht, double hx, unsigned int M, unsigned int N, double lambdaM, double lambdaL, double lambdaR UNUSED_PARAM, double a)
{
    const DoubleVector &k = *pk;
    const DoubleMatrix &u = *pu;

    psi.clear();
    psi.resize(M+1, N+1, 0.0);

    DoubleVector d1(N+1);
    DoubleVector rx(N+1);

    for (unsigned int m1=0; m1<=M; m1++)
    {
        unsigned int m = M-m1;

        if (m==M)
        {
            for (unsigned int n=0; n<=N; n++)
                psi.at(M,n) = -2.0*alpha1*mu(n)*(u.at(M,n) - V[n]);
        }
        else
        {
            DoubleMatrix p(N+1, N+1, 0.0);

            p.at(0,0) = -1.0 - (a*a*ht)/(hx*hx) - (lambdaL*a*a*ht)/hx - lambdaM*ht;
            p.at(0,1) = (a*a*ht)/(hx*hx);
            d1[0]     = -psi.at(m+1,0);

            //p(0, Xi[0]) = -k[0] * ((lambdaL*a*a*ht)/(hx));
            //p(0, Xi[1]) = -k[1] * ((lambdaL*a*a*ht)/(hx));

            for (unsigned int n=1; n<=N-1; n++)
            {
                p.at(n,n-1) = (a*a*ht)/(hx*hx);
                p.at(n,n)   = -1.0 - 2.0*((a*a)*ht)/(hx*hx) - lambdaM*ht;
                p.at(n,n+1) = (a*a*ht)/(hx*hx);

                if (n>1) p.at(n,0) = 0.0;
                if (n==Xi[0]) p.at(n,0) = p.at(n,0) + (1/hx)*lambdaM*a*a*ht*k[0];
                if (n==Xi[1]) p.at(n,0) = p.at(n,0) + (1/hx)*lambdaM*a*a*ht*k[1];

                d1[n] = -psi.at(m+1, n);
                for (unsigned int j=0; j<L; j++)
                {
                    if (n==Xi[j])
                    {
                        for (unsigned i=0; n<L; i++)
                        {
                            d1[n] += (1/hx)*2*alpha2*ht*k[j]*k[i]*(u.at(m,Xi[i]) - z[i]);
                        }
                    }
                }

                for (unsigned int i=0; i<L; i++)
                {
                    if (n==Xi[i])
                    {
                        d1[n] += (1/hx)*2*alpha2*ht*k[i]*k[i]*(u.at(m,n) - z[i]);
                    }
                }
            }

            p.at(N,N-1) = (a*a*ht)/(hx*hx);
            p.at(N,N)   = -1.0 - (a*a*ht)/(hx*hx) - (lambdaL*a*a*ht)/hx - lambdaM*ht;
            d1[N]       = -psi.at(m+1,N);

            GaussianElimination(p, d1, rx);

            for (unsigned int n=0; n<=N; n++) psi.at(m,n) = rx[n];
        }
    }
}

double Problem2::mu(unsigned int i UNUSED_PARAM) const
{
    return 1.0;
}

void Problem2::calculateV(const DoubleVector &k)
{
    pk = &k;
    DoubleMatrix u;
    calculateU(u, ht, hx, M, N, lambdaM, lambdaL, lambdaR, a);
    V = u.row(M);
}


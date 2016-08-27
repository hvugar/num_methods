#include "problem3.h"

void Problem3::Main(int argc UNUSED_PARAM, char *argv[] UNUSED_PARAM)
{
    Problem3 p;

    DoubleVector v0(p.M);
    for (unsigned int m=0; m<=p.M; m++) v0[m] = 1.0;

//    puts("Analitic");
//    DoubleVector ga(p.M,0.0);
//    p.gradient(v0, ga);
//    IPrinter::printVector(ga);
//    ga.L2Normalize();
//    IPrinter::printVector(ga);

//    puts("Numerical");
//    double h = 0.001;
//    DoubleVector gn(p.M,0.0);
//    IGradient::Gradient(&p, h, v0, gn);
//    IPrinter::printVector(gn);
//    gn.L2Normalize();
//    IPrinter::printVector(gn);

    ConjugateGradient g;
    g.setFunction(&p);
    g.setGradient(&p);
    g.setPrinter(&p);
    g.setEpsilon1(0.0001);
    g.setEpsilon2(0.0001);
    g.setEpsilon3(0.0001);
    g.setR1MinimizeEpsilon(0.1, 0.0001);
    g.setNormalize(true);
    g.calculate(v0);

}

Problem3::Problem3()
{
    t0 = 0.0; t1 = 1.0;
    x0 = 0.0; x1 = 1.0;
    hx = 0.001;
    ht = 0.001;
    N = 100;
    M = 100;
    a = 1.0;

    lambdaM = 1.0;
    lambdaL = 1.0;
    lambdaR = 1.0;

    alpha1 = 1.0;
    alpha2 = 1.0;

    Te = 3.0;

    DoubleVector vs(M+1);
    for (unsigned int j=0; j<=M; j++) vs[j] = 4.0;
    calculateV(vs);
    IPrinter::printVector(V);
}

Problem3::~Problem3()
{}

double Problem3::fx(const DoubleVector &v)
{
    pv = &v;
    DoubleMatrix u;
    calculateU(u, ht, hx, M, N, lambdaM, lambdaL, lambdaR, a);
    pu = &u;

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

        double vm1 = v.at(m1);
        double vm2 = v.at(m2);

        double vs1 = vs(m1);
        double vs2 = vs(m2);

        norm = norm + ((vm1-vs1)*(vm1-vs1) + (vm2-vs2)*(vm2-vs2));
    }
    norm = 0.5*ht*norm;

    return alpha1*sum + alpha2*norm;
}

void Problem3::gradient(const DoubleVector &v, DoubleVector &g)
{
    pv = &v;
    DoubleMatrix u;
    calculateU(u, ht, hx, M, N, lambdaM, lambdaL, lambdaR, a);
    pu = &u;

    DoubleMatrix psi;
    calculateP(psi, ht, hx, M, N, lambdaM, lambdaL, lambdaR, a);
    pp = &psi;

    for (unsigned int m=0; m<=M; m++)
    {
        g[m] = -lambdaM*a*a*psi.at(m, 0) + 2.0*alpha2*(v[m]-vs(m));
    }
}

void Problem3::print(unsigned int i, const DoubleVector &v, const DoubleVector &g, double alpha, RnFunction *fn) const
{
    C_UNUSED(alpha);

    Problem3 *hc = dynamic_cast<Problem3*>(fn);
    printf("J[%d]: %.16f\n", i, hc->fx(v));
    IPrinter::printVector(v);
    IPrinter::printVector(g);
}

double Problem3::initial(unsigned int i UNUSED_PARAM) const
{
    return 1.0;
}

double Problem3::v(unsigned int j UNUSED_PARAM) const
{
    const DoubleVector &_v = *pv;
    return _v[j];
}

void Problem3::calculateU(DoubleMatrix &u, double ht, double hx, unsigned int M, unsigned int N, double lambdaM, double lambdaL, double lambdaR, double a)
{
    u.clear();
    u.resize(M+1, N+1);

    DoubleVector rb(N+1);
    DoubleVector rx(N+1);

    double a1 = -(a*a*ht)/(hx*hx);
    double b1 = 1.0 + (a*a*ht)/(hx*hx) + (lambdaL*a*a*ht)/hx+lambdaM*ht;
    double c1 = 1.0 + 2.0*(a*a*ht)/(hx*hx) + lambdaM*ht;

    for (unsigned int j=0; j<=M; j++)
    {
        if (j==0)
        {
            for (unsigned int i=0; i<=N; i++) u.at(0,i) = initial(i);
        }
        else
        {
            DoubleMatrix ra(N+1, N+1, 0.0);


            ra(0,0) = b1;
            ra(0,1) = a1;
            rb[0] = u.at(j-1,0) + ((lambdaL*a*a*ht)/(hx))*v(j) + lambdaM*ht*Te;

            for (unsigned int i=1; i<=N-1; i++)
            {
                ra(i,i-1) = a1;
                ra(i,i)   = c1;
                ra(i,i+1) = a1;
                rb[i] = u.at(j-1, i) + lambdaM*ht*Te;
            }

            ra(N,N-1) = a1;
            ra(N,N)   = b1;
            rb[N] = u.at(j-1,N) + ((lambdaR*a*a*ht)/(hx))*Te + lambdaM*ht*Te;

            GaussianElimination(ra, rb, rx);

            for (unsigned int i=0; i<=N; i++) u.at(j,i) = rx[i];
        }
    }
}

void Problem3::calculateP(DoubleMatrix &psi, double ht, double hx, unsigned int M, unsigned int N, double lambdaM, double lambdaL, double lambdaR UNUSED_PARAM, double a)
{
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

            for (unsigned int n=1; n<=N-1; n++)
            {
                p.at(n,n-1) = (a*a*ht)/(hx*hx);
                p.at(n,n)   = -1.0 - 2.0*((a*a)*ht)/(hx*hx) - lambdaM*ht;
                p.at(n,n+1) = (a*a*ht)/(hx*hx);

                d1[n] = -psi.at(m+1, n);
            }

            p.at(N,N-1) = (a*a*ht)/(hx*hx);
            p.at(N,N)   = -1.0 - (a*a*ht)/(hx*hx) - (lambdaL*a*a*ht)/hx - lambdaM*ht;
            d1[N]       = -psi.at(m+1,N);

            GaussianElimination(p, d1, rx);

            for (unsigned int n=0; n<=N; n++) psi.at(m,n) = rx[n];
        }
    }
}

double Problem3::mu(unsigned int i UNUSED_PARAM) const
{
    return 1.0;
}

void Problem3::calculateV(const DoubleVector &v)
{
    pv = &v;
    DoubleMatrix u;
    calculateU(u, ht, hx, M, N, lambdaM, lambdaL, lambdaR, a);
    V = u.row(M);
}

double Problem3::vs(unsigned int j) const
{
    return 4.0;
}


#include "problem2.h"

void Problem2::Main(int argc UNUSED_PARAM, char *argv[] UNUSED_PARAM)
{
    Problem2 p;
    DoubleMatrix u;
    p.calculateU(u, p.ht, p.hx, p.M, p.N, p.lambdaM, p.lambdaL, p.lambdaR, p.a);
    IPrinter::printMatrix(u);
}

Problem2::Problem2()
{
    t0 = 0.0; t1 = 1.0;
    x0 = 0.0; x1 = 1.0;
    hx = 0.01;
    ht = 0.001;
    N = 100;
    M = 1000;
    a = 1.0;

    lambdaM = 1.0;
    lambdaL = 1.0;
    lambdaR = 1.0;

    L = 2;

    k.resize(L);
    k[0] = 1.0;
    k[1] = 1.0;

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

    V.resize(N+1);
    DoubleMatrix u;
    calculateU(u, ht, hx, M, N, lambdaM, lambdaL, lambdaR, a);
    V = u.row(M);
    IPrinter::printVector(V);
}

Problem2::~Problem2()
{}

double Problem2::fx(const DoubleVector &x)
{
    k = x;
    DoubleMatrix u;
    calculate(u, ht, hx, M, N, lambdaM, lambdaL, lambdaR, a);
    pu = &u;

    double sum = 0.0;
    for (unsigned int i=0; i<=N-1; i++)
    {
        unsigned int i1 = i;
        unsigned int i2 = i + 1;
        double f1 = mu(i1) * (u(M, i1) - V[i1])*(u(M, i1) - V[i1]);
        double f2 = mu(i2) * (u(M, i2) - V[i2])*(u(M, i2) - V[i2]);
        sum += (f1 + f2);
    }
    sum *= 0.5*hx;

    double norm = 0.0;
    for (unsigned int j=0; j<=M-1; j++)
    {
        unsigned int j1 = j;
        unsigned int j2 = j + 1;
        sum += (vl(j1) + vl(j2));
    }
    norm *= 0.5*ht;

    return alpha1*sum + alpha2*norm;
}

void Problem2::gradient(const DoubleVector &x, DoubleVector &g)
{
}

double Problem2::initial(unsigned int i UNUSED_PARAM) const { return 1.0; }

double Problem2::vm(unsigned int j UNUSED_PARAM) const { return Te; }

double Problem2::vl(unsigned int j UNUSED_PARAM) const
{
    const DoubleMatrix &u =*pu;
    return k[0] * (u.at(j, Xi[0]) - z[0]) + k[1] * (u.at(j, Xi[1]) - z[1]);
}

double Problem2::vr(unsigned int j UNUSED_PARAM) const { return 3.0; }

void Problem2::calculateU(DoubleMatrix &m, double ht, double hx, unsigned int M, unsigned int N, double lambdaM, double lambdaL, double lambdaR, double a)
{
    m.clear();
    m.resize(M+1, N+1);

    DoubleVector a1(N+1);
    DoubleVector b1(N+1);
    DoubleVector c1(N+1);
    DoubleVector d1(N+1);
    DoubleVector rx(N+1);

    for (unsigned int j=0; j<=M; j++)
    {
        if (j==0)
        {
            for (unsigned int i=0; i<=N; i++) m.at(0,i) = initial(i);
        }
        else
        {
            DoubleMatrix u(N+1, N+1, 0.0);
            u(0,0) = 1.0 + (a*a*ht)/(hx*hx) + (lambdaL*a*a*ht)/hx + lambdaM*ht;
            u(0,1) = -(a*a*ht)/(hx*hx);
            d1[0] = m.at(j-1,0) + ((lambdaL*a*a*ht)/(hx))*(-(k[0]*z[0]-k[1]*z[1])) + lambdaM*ht*vm(j);

            u(0, Xi[0]) = -k[0] * ((lambdaL*a*a*ht)/(hx));
            u(0, Xi[1]) = -k[1] * ((lambdaL*a*a*ht)/(hx));

            for (unsigned int i=1; i<=N-1; i++)
            {
                u(i,i-1) = -(a*a*ht)/(hx*hx);
                u(i,i) = 1.0 + 2.0*((a*a)*ht)/(hx*hx) + lambdaM*ht;
                u(i,i+1) = -(a*a*ht)/(hx*hx);
                d1[i] = m.at(j-1, i) + lambdaM*ht*vm(j);
            }

            u(N,N-1) = -(a*a*ht)/(hx*hx);
            u(N,N)   = 1.0 + (a*a*ht)/(hx*hx) + (lambdaL*a*a*ht)/hx + lambdaM*ht;
            d1[N] = m.at(j-1,N) + ((lambdaR*a*a*ht)/(hx))*vr(j) + lambdaM*ht*vm(j);

            GaussianElimination(u, d1, rx);

//            a1[0] = 0.0;
//            b1[0] = 1.0 + (a*a*ht)/(hx*hx) + (lambdaL*a*a*ht)/hx + lambdaM*ht;
//            c1[0] = -(a*a*ht)/(hx*hx);
//            d1[0] = m.at(j-1,0) + ((lambdaL*a*a*ht)/(hx))*vl(j) + lambdaM*ht*vm(j);

//            for (unsigned int i=1; i<=N-1; i++)
//            {
//                a1[i] = -(a*a*ht)/(hx*hx);
//                b1[i] = 1.0 + 2.0*((a*a)*ht)/(hx*hx) + lambdaM*ht;
//                c1[i] = -(a*a*ht)/(hx*hx);
//                d1[i] = m.at(j-1, i) + lambdaM*ht*vm(j);
//            }

//            a1[N] = -(a*a*ht)/(hx*hx);
//            b1[N] = 1.0 + (a*a*ht)/(hx*hx) + lambdaR*(a*a*ht)/hx + lambdaM*ht;
//            c1[N] = 0.0;
//            d1[N] = m.at(j-1,N) + ((lambdaR*a*a*ht)/(hx))*vr(j) + lambdaM*ht*vm(j);

//            tomasAlgorithm(a1.data(), b1.data(), c1.data(), d1.data(), rx.data(), rx.size());

            for (unsigned int i=0; i<=N; i++) m.at(j,i) = rx[i];
        }
    }
}

void Problem2::calculateP(DoubleMatrix &psi, double ht, double hx, unsigned int M, unsigned int N, double lambdaM, double lambdaL, double lambdaR, double a)
{
    const DoubleMatrix &u1 = *pu;

    psi.clear();
    psi.resize(M+1, N+1);

    DoubleVector d1(N+1);
    DoubleVector rx(N+1);

    for (unsigned int j1=0; j1<=M; j1++)
    {
        unsigned int j = M-j1;

        if (j==M)
        {
            for (unsigned int i=0; i<=N; i++) psi.at(0,i) = -2.0*alpha1*mu(i)(u1(M,i) - V(i));
        }
        else
        {
            DoubleMatrix u(N+1, N+1, 0.0);
            u(0,0) = 1.0 + (a*a*ht)/(hx*hx) + (lambdaL*a*a*ht)/hx + lambdaM*ht;
            u(0,1) = -(a*a*ht)/(hx*hx);
            d1[0] = psi.at(j-1,0) + ((lambdaL*a*a*ht)/(hx))*(-(k[0]*z[0]-k[1]*z[1])) + lambdaM*ht*vm(j);

            u(0, Xi[0]) = -k[0] * ((lambdaL*a*a*ht)/(hx));
            u(0, Xi[1]) = -k[1] * ((lambdaL*a*a*ht)/(hx));

            for (unsigned int i=1; i<=N-1; i++)
            {
                u(i,i-1) = -(a*a*ht)/(hx*hx);
                u(i,i) = 1.0 + 2.0*((a*a)*ht)/(hx*hx) + lambdaM*ht;
                u(i,i+1) = -(a*a*ht)/(hx*hx);
                d1[i] = psi.at(j-1, i) + lambdaM*ht*vm(j);
            }

            u(N,N-1) = -(a*a*ht)/(hx*hx);
            u(N,N)   = 1.0 + (a*a*ht)/(hx*hx) + (lambdaL*a*a*ht)/hx + lambdaM*ht;
            d1[N] = psi.at(j-1,N) + ((lambdaR*a*a*ht)/(hx))*vr(j) + lambdaM*ht*vm(j);

            GaussianElimination(u, d1, rx);

//            a1[0] = 0.0;
//            b1[0] = 1.0 + (a*a*ht)/(hx*hx) + (lambdaL*a*a*ht)/hx + lambdaM*ht;
//            c1[0] = -(a*a*ht)/(hx*hx);
//            d1[0] = m.at(j-1,0) + ((lambdaL*a*a*ht)/(hx))*vl(j) + lambdaM*ht*vm(j);

//            for (unsigned int i=1; i<=N-1; i++)
//            {
//                a1[i] = -(a*a*ht)/(hx*hx);
//                b1[i] = 1.0 + 2.0*((a*a)*ht)/(hx*hx) + lambdaM*ht;
//                c1[i] = -(a*a*ht)/(hx*hx);
//                d1[i] = m.at(j-1, i) + lambdaM*ht*vm(j);
//            }

//            a1[N] = -(a*a*ht)/(hx*hx);
//            b1[N] = 1.0 + (a*a*ht)/(hx*hx) + lambdaR*(a*a*ht)/hx + lambdaM*ht;
//            c1[N] = 0.0;
//            d1[N] = m.at(j-1,N) + ((lambdaR*a*a*ht)/(hx))*vr(j) + lambdaM*ht*vm(j);

//            tomasAlgorithm(a1.data(), b1.data(), c1.data(), d1.data(), rx.data(), rx.size());

            for (unsigned int i=0; i<=N; i++) psi.at(j,i) = rx[i];
        }
    }
}

double Problem2::mu(unsigned int i UNUSED_PARAM) const
{
    return 1.0;
}


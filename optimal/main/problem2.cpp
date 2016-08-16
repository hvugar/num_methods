#include "problem2.h"

Problem2::Problem2()
{
    t0 = 0.0; t1 = 1.0;
    x0 = 0.0; x1 = 1.0;
    hx = ht = 0.001;
    N = M = 1000;
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

    z.resize(L);
    z[0] = 2.89;
    z[1] = 2.81;
}

Problem2::~Problem2()
{}

double Problem2::initial(unsigned int i UNUSED_PARAM) const { return 1.0; }

double Problem2::vm(unsigned int j UNUSED_PARAM) const { return 3.0; }

double Problem2::vl(unsigned int j UNUSED_PARAM) const
{
//    double v = 0.0;
//    for (unsigned int i=0; i<L; i++)
//    {
//        v = v + k[i] * z[i];
//    }
//    return -v;
    return 4.0;
}

double Problem2::vr(unsigned int j UNUSED_PARAM) const { return 3.0; }

void Problem2::calculate1(DoubleMatrix &m, double ht, double hx, unsigned int M, unsigned int N, double lambdaM, double lambdaL, double lambdaR, double a)
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
            DoubleMatrix u(N+1,N+1,0.0);
            u(0,0) = 1.0 + (a*a*ht)/(hx*hx) + (lambdaL*a*a*ht)/hx + lambdaM*ht;
            u(0,1) = -(a*a*ht)/(hx*hx);
            d1[0] = m.at(j-1,0) + ((lambdaL*a*a*ht)/(hx))*vl(j) + lambdaM*ht*vm(j);

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


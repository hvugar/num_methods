#include "newtonheatprocess.h"
#include <cmethods.h>

void NewtonHeatProcess::calculate(DoubleMatrix &m, double ht, double hx, unsigned int M, unsigned int N, double lambdaM, double lambdaL, double lambdaR, double a)
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
            a1[0] = 0.0;
            b1[0] = 1.0 + (a*a*ht)/(hx*hx) + (lambdaL*a*a*ht)/hx + lambdaM*ht;
            c1[0] = -(a*a*ht)/(hx*hx);
            d1[0] = m.at(j-1,0) + ((lambdaL*a*a*ht)/(hx))*vl(j) + lambdaM*ht*vm(j);

            for (unsigned int i=1; i<=N-1; i++)
            {
                a1[i] = -(a*a*ht)/(hx*hx);
                b1[i] = 1.0 + 2.0*((a*a)*ht)/(hx*hx) + lambdaM*ht;
                c1[i] = -(a*a*ht)/(hx*hx);
                d1[i] = m.at(j-1, i) + lambdaM*ht*vm(j);
            }

            a1[N] = -(a*a*ht)/(hx*hx);
            b1[N] = 1.0 + (a*a*ht)/(hx*hx) + lambdaR*(a*a*ht)/hx + lambdaM*ht;
            c1[N] = 0.0;
            d1[N] = m.at(j-1,N) + ((lambdaR*a*a*ht)/(hx))*vr(j) + lambdaM*ht*vm(j);

            tomasAlgorithm(a1.data(), b1.data(), c1.data(), d1.data(), rx.data(), rx.size());

            for (unsigned int i=0; i<=N; i++) m.at(j,i) = rx[i];
        }
    }
}

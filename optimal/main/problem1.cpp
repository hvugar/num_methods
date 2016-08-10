#include "problem1.h"

Problem1::Problem1()
{
    //t0 = 0.0;
    //t1 = 2.0;
    //x0 = 0.0;
    //x1 = 1.0;
    a = 1.0;
    lambda = -1.0;
    Te = 1.0;
    hx = 0.001;
    ht = 0.001;
    N = 1000;
    M = 10000;
}

double Problem1::initial(unsigned int i) const
{
    double x = i*ht;
    return x*x*x;
}

double Problem1::boundary(Boundary type, unsigned int j) const
{
    double t = j*ht;
    if (type == Left) return lambda * (t*t - u.at(0, j));
    if (type == Right) return lambda * (u.at(0, j));
    return 0.0;
}

double Problem1::f(unsigned int i, unsigned int j) const
{
    return 0.0;
}

double Problem1::v(unsigned int j) const
{
    return 2.0;
}

void Problem1::calculate()
{
    DoubleMatrix u(M+1, N+1);

    DoubleVector a1(N+1);
    DoubleVector b1(N+1);
    DoubleVector c1(N+1);
    DoubleVector d1(N+1);
    DoubleVector rx(N+1);

    for (uint32_t j=0; j<=M; j++)
    {
        if (j==0)
        {
            for (uint32_t i=0; i<=N; i++)
                u.at(j,i) = Te;
        }
        else
        {
            a1[0] = 0.0;
            b1[0] = 1.0 + (a*a*ht)/(hx*hx) + (a*a*ht*lambda)/hx - lambda*ht;
            c1[0] = -(a*a*ht)/(hx*hx);
            d1[0] = u.at(j-1,0) + v(j)*((lambda*(a*a)*ht)/(hx)) - lambda*ht*Te;

            for (uint32_t i=1; i<=N-1; i++)
            {
                a1[i] = -(a*a*ht)/(hx*hx);
                b1[i] = 1.0 + 2.0*((a*a)*ht)/(hx*hx) - ht*lambda;
                c1[i] = -(a*a*ht)/(hx*hx);
                d1[i] = u.at(j-1, i) - lambda*ht*Te;
            }

            a1[N] = -(a*a*ht)/(hx*hx);
            b1[N] = 1.0 + (a*a*ht)/(hx*hx) - (a*a*ht*lambda)/hx - lambda*ht;
            c1[N] = 0.0;
            d1[N] = u.at(j-1,N) - ((lambda*(a*a)*ht)/(hx))*Te - lambda*ht*Te;

            tomasAlgorithm(a1.data(), b1.data(), c1.data(), d1.data(), rx.data(), rx.size());

            for (uint32_t i=0; i<=N; i++) u.at(j,i) = rx[i];
        }
    }

    IPrinter::printVector(u[0]);
    IPrinter::printVector(u[1]);
    IPrinter::printVector(u[2]);
    IPrinter::printMatrix(u);
}


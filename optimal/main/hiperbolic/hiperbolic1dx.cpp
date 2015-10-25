#include "hiperbolic1dx.h"
#include <tomasmethod.h>
#include <stdlib.h>

Hiperbolic1DX::Hiperbolic1DX()
{
    t0 = 0.0;
    t1 = 1.0;
    x0 = 0.0;
    x1 = 1.0;
    M = 1000;
    N = 1000;
    hx = (x1-x0)/N;
    ht = (t1-t0)/M;
    a = 1.0;
    lamda = 0.25;
}

Hiperbolic1DX::~Hiperbolic1DX()
{

}

double Hiperbolic1DX::fx(const DoubleVector &e)
{
    return 0.0;
}

void Hiperbolic1DX::gradient(const DoubleVector &e, DoubleVector &g, double gradient_step)
{

}

void Hiperbolic1DX::project(DoubleVector &e, int i)
{

}

void Hiperbolic1DX::print(unsigned int i, const DoubleVector &e, const DoubleVector &g, double a, RnFunction *fn) const
{

}

void Hiperbolic1DX::calculateU(const DoubleVector &e, DoubleVector &ut)
{
    ut.clear();
    ut.resize(N+1);

    DoubleVector u0(N+1);
    DoubleVector u1(N+1);
    DoubleVector U1(N+1);

    DoubleVector a1;
    DoubleVector b1;
    DoubleVector c1;
    DoubleVector d1;
    DoubleVector x1;

    a1.resize(N-1);
    b1.resize(N-1);
    c1.resize(N-1);
    d1.resize(N-1);
    x1.resize(N-1);

    double alpha1 = -0.25;//-(lamda*a*a)*((ht*ht)/(hx*hx));
    double beta1  = 1.5;//1.0 + (2.0*lamda*a*a*(ht*ht))/(hx*hx);
    double alpha2 = 0.5;//(1.0-2.0*lamda)*a*a*((ht*ht)/(hx*hx));
    double alpha3 = 0.25;//+(lamda*a*a)*((ht*ht)/(hx*hx));

    for (unsigned int j=0; j<M; j++)
    {
        switch (j)
        {
        case 0:
        {
            for (unsigned int i=0; i<=N; i++)
            {
                u0[i] = fi1(i*hx);
            }
            Printer::printVector(u0, N/10);
        }
            break;
        case 1:
        {
            for (unsigned int i=0; i<=N; i++)
            {
                u1[i] = u0[i] + ht*fi2(i*hx);
            }
            Printer::printVector(u1, N/10);
        }
            break;
        default:
        {
            for (unsigned int i=1; i<=N-1; i++)
            {
                a1[i-1] = alpha1;
                b1[i-1] = beta1;
                c1[i-1] = alpha1;
                d1[i-1] = alpha2*(u1[i-1]-2.0*u1[i]+u1[i+1]) + 2.0*u1[i] + alpha3*(u0[i+1] - 2.0*u0[i] + u0[i-1]) + (ht*ht)*f(i*hx, j*ht);
            }

            a1[0]   = 0.0;
            c1[N-2] = 0.0;
            d1[0]   -= alpha1 * mu1((j+1)*ht);
            d1[N-2] -= alpha1 * mu2((j+1)*ht);
            TomasAlgorithm(a1, b1, c1, d1, x1);

            ut[0] = mu1((j+1)*ht);
            for (unsigned int i=1; i<=N-1; i++)
            {
                ut[i] = x1[i-1];
            }
            ut[N] = mu2((j+1)*ht);

            for (unsigned int i=0; i<=N; i++)
            {
                u0[i] = u1[i];
                u1[i] = ut[i];
                U1[i] = u(i*hx, j*ht);
            }

            Printer::printVector(ut, N/10);
            Printer::printVector(U1, N/10);
            exit(-1);
        }
            break;
        }
    }

    a1.clear();
    b1.clear();
    c1.clear();
    d1.clear();
    x1.clear();

    Printer::printVector(ut, N/10);
}

double Hiperbolic1DX::u(double x, double t)  { return x*x*x*x + t*t*t*t; }
double Hiperbolic1DX::fi1(double x) { return u(x, 0.0); }
double Hiperbolic1DX::fi2(double x) { return 0.0/*(u(x, +0.000001) - u(x, -0.000001))/(2.0*0.000001)*/; }
double Hiperbolic1DX::mu1(double t) { return u(0.0, t); }
double Hiperbolic1DX::mu2(double t) { return u(1.0, t); }
double Hiperbolic1DX::f(double x, double t) { return 12.0*(t*t-x*x*a*a); }

void Hiperbolic1DX::main()
{
    DoubleVector e;
    DoubleVector u;

    Hiperbolic1DX h;
    h.calculateU(e, u);

    //    unsigned int M = 1000;
    //    unsigned int N = 100;
    //    double hx = 1.0 / N;
    //    double ht = 1.0 / M;

    //    DoubleMatrix U;

    //    U.resize(M+1); for (unsigned int j=0; j<=M; j++) U[j].resize(N+1);

    //    for (unsigned int j=0; j<=M; j++)
    //    {
    //        for (unsigned int i=0; i<=N; i++)
    //        {
    //            U[j][i] = u(i*hx, j*ht);
    //        }
    //    }

    //    DoubleVector u0(N+1);
    //    for (unsigned int i=0; i<=N; i++)  u0[i] = fi1(i*hx);
    //    DoubleVector u1(N+1);
    //    for (unsigned int i=0; i<=N; i++)  u1[i] = fi1(i*hx) + ht * fi2(i*hx);

    //    //Printer::printMatrix(U, M/10, N/10);
    //    puts("------------------------------------------------");
    //    Printer::printVector(U[0], N/10);
    //    Printer::printVector(u0, N/10);
    //    puts("------------------------------------------------");
    //    Printer::printVector(U[1], N/10);
    //    Printer::printVector(u1, N/10);
    //    puts("------------------------------------------------");

    //    DoubleVector u2(N+1);
    //    double alpha = (ht*ht)/(hx*hx);
    //    double betta = 2.0*(1.0 - alpha);

    //    for (unsigned int j=2; j<=2; j++)
    //    {
    //        u2[0] = mu1(j*ht);
    //        for (unsigned int i=1; i<=N-1; i++)
    //        {
    //            u2[i] = alpha*u1[i+1] + betta*u1[i] + alpha*u1[i-1] + u0[i] + (ht*ht)*f(i*hx, j*ht);
    //        }
    //        u2[N] = mu2(j*ht);
    //        for (unsigned int i=0; i<=N; i++)
    //        {
    //            u0[i] = u1[i];
    //            u1[i] = u2[i];
    //        }
    //        Printer::printVector(U[j], N/10);
    //        Printer::printVector(u2, N/10);
    //        puts("------------------------------------------------");
    //    }
    //    puts("------------------------------------------------");
    //    Printer::printMatrix(U, M/10, N/10);
}


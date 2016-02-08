#include "discretehyperbolic.h"
#include <tomasmethod.h>
#include <gradient_cjt.h>

void DiscreteHyperbolic::main()
{
    DiscreteHyperbolic dh;

    DoubleVector f0((dh.M+1)*(dh.N+1));
    for (unsigned int j=0; j<=dh.M; j++)
    {
        for (unsigned int i=0; i<=dh.N; i++)
        {
            //double x = i*dh.hx;
            //double t = j*dh.ht;
            f0[j*(dh.N+1)+i] = 0.0;//6.0*t - 6.0*x*dh.a;
        }
    }

    //DoubleVector g(f0.size());
    //dh.gradient(f0, g);

    ConjugateGradient g2;
    g2.setFunction(&dh);
    //    g2.setEpsilon1(0.000000001);
    //    g2.setEpsilon2(0.000000001);
    g2.setEpsilon1(0.01);
    g2.setEpsilon2(0.01);
    //g2.setGradientStep(0.001);
    g2.setR1MinimizeEpsilon(0.1, 0.00001);
    g2.setPrinter(&dh);
    g2.calculate(f0);

    IPrinter::printAsMatrix(f0, dh.M, dh.N);
    DoubleVector u;
    dh.pf = &f0;
    dh.calculateU(u, dh.hx, dh.ht, dh.M, dh.N);
    puts("-----");
    IPrinter::printVector(dh.U);
    IPrinter::printVector(u);
}

DiscreteHyperbolic::DiscreteHyperbolic()
{
    t0 = x0 = 0.0;
    t1 = x1 = 1.0;
    M = 100;
    N = 100;
    ht = (t1 - t0) / M;
    hx = (x1 - x0) / N;
    a = 1.0;
    lamda = 0.25;
    U.resize(N+1);
    for (unsigned int i=0; i<=N; i++)
    {
        double x = i*hx;
        U[i] = x*x*x + 1.0;
    }
}

double DiscreteHyperbolic::fx(const DoubleVector& f0)
{
    pf = &f0;

    double sum = 0.0;

    DoubleVector u;
    IHyperbolicEquation::calculateU(u, hx, ht, M, N, a);
    for (unsigned int i=0; i<=N; i++)
    {
        double b = 1.0;
        if (i==0 || i==N) b = 0.5;
        sum += b*(u[i]-U[i])*(u[i]-U[i]);
    }
    sum = hx*sum;

    double norm = 0.0;
    for (unsigned int j=0; j<=M; j++)
    {
        for (unsigned int i=0; i<=N; i++)
        {
            double alpha = 1.0;
            if (i==0 || i==N || j==0 || j==M) alpha = 0.5;
            if (i==0 && j==0) alpha = 0.25;
            if (i==0 && j==M) alpha = 0.25;
            if (i==N && j==0) alpha = 0.25;
            if (i==N && j==M) alpha = 0.25;
            double f1 = (f0[j*(N+1)+i] - F(i, j));
            norm += alpha*f1*f1;
        }
    }
    norm = hx*ht*norm;

    return sum+norm;
}

void DiscreteHyperbolic::gradient(const DoubleVector& f0, DoubleVector& g)
{
    pf = &f0;
    DoubleMatrix u;
    IHyperbolicEquation::calculateU(u, hx, ht, M, N, a, lamda);
    DoubleMatrix psi;
    calculateP(f0, u, psi, g);
}

void DiscreteHyperbolic::print(unsigned int i, const DoubleVector &x, const DoubleVector &g , double alpha, RnFunction *fn) const
{
    C_UNUSED(alpha);
    C_UNUSED(g);
    printf("J[%d]: %.12f\n", i, fn->fx(x));
}

double DiscreteHyperbolic::fi1(unsigned int i) const
{
    double x = i*hx;
    return x*x*x;
}

double DiscreteHyperbolic::fi2(unsigned int i) const
{
    return 0.0;
}

double DiscreteHyperbolic::m1(unsigned int j) const
{
    double t = j*ht;
    return t*t*t;
}

double DiscreteHyperbolic::m2(unsigned int j) const
{
    double t = j*ht;
    return t*t*t+1.0;
}

double DiscreteHyperbolic::f(unsigned int i, unsigned int j) const
{
    return (*pf)[j*(N+1)+i];
    //    double x = i*hx;
    //    double t = j*ht;
    //    return 6.0*t - 6.0*x*a;
}

double DiscreteHyperbolic::F(unsigned int i, unsigned int j) const
{
    double x = i*hx;
    double t = j*ht;
    return 6.0*t - 6.0*x*a;
}

void DiscreteHyperbolic::calculateP(const DoubleVector& f0, const DoubleMatrix &u, DoubleMatrix &psi, DoubleVector &g)
{
    for (unsigned int j=0; j<psi.size(); j++) psi[j].clear();
    psi.clear();

    psi.resize(M+1);
    for (unsigned int j=0; j<psi.size(); j++) psi[j].resize(N+1);

    double A1 = -(lamda*a*a*ht*ht)/(hx*hx);
    double B0 = 1.0 + (2.0*lamda*a*a*ht*ht)/(hx*hx);
    double A2 = -(lamda*a*a*ht*ht)/(hx*hx);

    double C1 = -(1.0-2.0*lamda)*(a*a*ht*ht)/(hx*hx);
    double D0 = 2.0*(((1.0-2.0*lamda)*(a*a*ht*ht) - hx*hx)/(hx*hx));
    double C2 = -(1.0-2.0*lamda)*(a*a*ht*ht)/(hx*hx);

    double E1 = -(lamda*a*a*ht*ht)/(hx*hx);
    double F0 = 1.0 + (2.0*lamda*a*a*ht*ht)/(hx*hx);
    double E2 = -(lamda*a*a*ht*ht)/(hx*hx);
    double G0 = -ht*ht;

    A1 /= G0;
    B0 /= G0;
    A2 /= G0;
    C1 /= G0;
    D0 /= G0;
    C2 /= G0;
    E1 /= G0;
    F0 /= G0;
    E2 /= G0;
    G0 = 1.0;

    DoubleVector da(N-1);
    DoubleVector db(N-1);
    DoubleVector dc(N-1);
    DoubleVector rd(N-1);
    DoubleVector rx(N-1);

    for (unsigned int j1=0; j1<=M; j1++)
    {
        unsigned int j = M-j1;

        if (j==M)
        {
            for (unsigned int i=1; i<=N-1; i++)
            {
                da[i-1] = A1;
                db[i-1] = B0;
                dc[i-1] = A2;
                rd[i-1] = -2.0*hx*1.0*(u[M][i]-U[i]);
            }
            da[0]=0.0;
            dc[N-2]=0.0;
            TomasAlgorithm(da, db, dc, rd, rx);
            for (unsigned int i=1; i<=N-1; i++)
            {
                psi[M][i] = rx[i-1];
            }
            psi[M][0]   = -(A2*psi[M][1]   + 2.0*hx*0.5*(u[M][0]-U[0]));
            psi[M][N]   = -(A1*psi[M][N-1] + 2.0*hx*0.5*(u[M][N]-U[N]));
        }
        else if (j==(M-1))
        {
            for (unsigned int i=1; i<=N-1; i++)
            {
                da[i-1] = A1;
                db[i-1] = B0;
                dc[i-1] = A2;
                rd[i-1] = 0.0;

                if (i==1)
                {
                    rd[i-1] = -(D0*psi[M][1] + C2*psi[M][2]);
                }
                else if (i==N-1)
                {
                    rd[i-1] = -(C1*psi[M][N-2] + D0*psi[M][N-1]);
                }
                else
                {
                    rd[i-1] = -(C1*psi[M][i-1] + D0*psi[M][i] + C2*psi[M][i+1]);
                }
            }
            da[0]=0.0;
            dc[N-2]=0.0;
            TomasAlgorithm(da, db, dc, rd, rx);
            for (unsigned int i=1; i<=N-1; i++)
            {
                psi[M-1][i] = rx[i-1];
            }
            psi[M-1][0] = -(A2*psi[M-1][1]  +C2*psi[M][1]);
            psi[M-1][N] = -(A1*psi[M-1][N-1]+C1*psi[M][N-1]);
        }
        else
        {
            if (j>=2)
            {
                for (unsigned int i=1; i<=N-1; i++)
                {
                    da[i-1] = A1;
                    db[i-1] = B0;
                    dc[i-1] = A2;
                    rd[i-1] = 0.0;

                    if (i==1)
                    {
                        rd[i-1] = -(D0*psi[j+1][1] + C2*psi[j+1][2] + F0*psi[j+2][1] + E2*psi[j+2][2]);
                    }
                    else if (i==N-1)
                    {
                        rd[i-1] = -(C1*psi[j+1][N-2] + D0*psi[j+1][N-1] + E1*psi[j+2][N-2] + F0*psi[j+2][N-1]);
                    }
                    else
                    {
                        rd[i-1] = -(C1*psi[j+1][i-1] + D0*psi[j+1][i] + C2*psi[j+1][i+1] + E1*psi[j+2][i-1] + F0*psi[j+2][i] + E2*psi[j+2][i+1]);
                    }
                }

                da[0]=0.0;
                dc[N-2]=0.0;

                TomasAlgorithm(da, db, dc, rd, rx);

                for (unsigned int i=1; i<=N-1; i++)
                {
                    psi[j][i] = rx[i-1];
                }

                psi[j][0] = -(A2*psi[j][1]   + C2*psi[j+1][1]   + E2*psi[j+2][1]);
                psi[j][N] = -(A1*psi[j][N-1] + C1*psi[j+1][N-1] + E2*psi[j+2][N-1]);
            }
            if (j==1)
            {
                psi[1][0] = -(C2*psi[2][1] + E2*psi[3][1]);
                for (unsigned int i=1; i<=N-1; i++)
                {
                    if (i==1)
                    {
                        psi[1][1] = -(D0*psi[2][1] + C2*psi[2][2] + F0*psi[3][1] + E2*psi[3][2]);
                    }
                    else if (i==N-1)
                    {
                        psi[1][N-1] = -(C1*psi[2][N-2] + D0*psi[2][N-1] + E1*psi[3][N-2] + F0*psi[3][N-1]);
                    }
                    else
                    {
                        psi[1][i] = -(C1*psi[2][i-1] + D0*psi[2][i] + C2*psi[2][i+1] + E1*psi[3][i-1] + F0*psi[3][i] + E2*psi[3][i+1]);
                    }
                }
                psi[1][N] = -(C1*psi[2][N-1] + E1*psi[3][N-1]);
            }
            if (j==0)
            {
                psi[0][0] = psi[1][0] - E2*psi[2][1];
                for (unsigned int i=1; i<=N-1; i++)
                {
                    if (i==1)
                    {
                        psi[0][1] = psi[1][1] - F0*psi[2][1] - E2*psi[2][2];
                    }
                    else if (i==N-1)
                    {
                        psi[0][N-1] = psi[1][N-1] - E1*psi[2][N-2] - F0*psi[2][N-1];
                    }
                    else
                    {
                        psi[0][i] = psi[1][i] - E1*psi[2][i-1] - F0*psi[2][i] - E2*psi[2][i+1];
                    }
                }
                psi[0][N] = psi[1][N]-E1*psi[2][N-1];
            }
        }

    }

    rx.clear();
    rd.clear();
    dc.clear();
    db.clear();
    da.clear();

    for (unsigned int j=0; j<=M; j++)
    {
        for (unsigned int i=0; i<=N; i++)
        {
            double alpha = 1.0;
            if (i==0 || i==N || j==0 || j==M) alpha = 0.5;
            if (i==0 && j==0) alpha = 0.25;
            if (i==0 && j==M) alpha = 0.25;
            if (i==N && j==0) alpha = 0.25;
            if (i==N && j==M) alpha = 0.25;

            unsigned int k = j*(N+1)+i;
            g[k] = 2.0*alpha*(f0[k]-F(i, j));
            //            g[k] = 2.0*alpha*hx*ht*(f0[k]-F(i, j));
            if (j>=2 && 1<=i && i<=(N-1))
                g[k] += psi[j][i];
        }
    }
}

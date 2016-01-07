#include "discretehyperbolic1.h"
#include <tomasmethod.h>
#include <gradient_cjt.h>

void Discretehyperbolic1::main()
{
    Discretehyperbolic1 dh;

    DoubleVector f0((dh.M+1)*(dh.N+1));
    for (unsigned int j=0; j<=dh.M; j++)
    {
        for (unsigned int i=0; i<=dh.N; i++)
        {
            double x = i*dh.hx;
            double t = j*dh.ht;
            f0[j*(dh.N+1)+i] = 6.0;//6.0*t - 6.0*x*dh.a;
        }
    }

    //DoubleVector g(f0.size());
    //dh.gradient(f0, g);

    ConjugateGradient g2;
    g2.setFunction(&dh);
    g2.setEpsilon1(0.0001);
    g2.setEpsilon2(0.0001);
    g2.setGradientStep(0.001);
    g2.setR1MinimizeEpsilon(0.1, 0.000001);
    g2.setPrinter(&dh);
    g2.calculate(f0);

    Printer::printAsMatrix(f0, dh.M, dh.N);


//    DoubleMatrix um;
//    dh.calculateU(um, dh.hx, dh.ht, dh.M, dh.N, dh.a, dh.lamda);
//    DoubleVector uv;
//    dh.calculateU(uv, dh.hx, dh.ht, dh.M, dh.N, dh.a, dh.lamda);
//    Printer::printMatrix(um);
//    Printer::printVector(uv);
//    printf("%.16f\n", dh.fx(f0));
//    //Printer::printAsMatrix(f0, dh.M, dh.N);
//    DoubleMatrix psi;
//    dh.calculateP(uv, psi);
//    //Printer::printMatrix(psi, 100, 10);
}

Discretehyperbolic1::Discretehyperbolic1()
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

double Discretehyperbolic1::fx(const DoubleVector& f0)
{
    pf = &f0;

    double sum = 0.0;

    DoubleVector u;
    calculateU(u, hx, ht, M, N, a);
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
            double b = 1.0;
            if (i==0 || i==N || j==0 || j==M) b = 0.5;
            if ((i==0 && j==0) || (i==0 && j==M) || (i==N && j==0) || (i==N && j==M)) b = 0.25;
            double f1 = (f0[j*(N+1)+i] - F(i, j));
            norm += b*f1*f1;
        }
    }
    norm = hx*ht*norm;

    return sum+norm;
}

void Discretehyperbolic1::gradient(const DoubleVector& f0, DoubleVector& g, double)
{
    pf = &f0;
    double G0 = -ht*ht;

    DoubleVector u;
    calculateU(u, hx, ht, M, N, a, lamda);
    //Printer::printVector(u);
    DoubleMatrix psi;
    calculateP(u, psi);
    //Printer::printMatrix(psi);

    for (unsigned int j=0; j<=M; j++)
    {
        for (unsigned int i=0; i<=N; i++)
        {
            unsigned int k = j*(N+1)+i;
            if (j==0)
            {
                if (i==0) g[k] = 2.0*0.25*hx*ht*(f0[k]-F(i, j));
                if (i==N) g[k] = 2.0*0.25*hx*ht*(f0[k]-F(i, j));
                if (0<i && i<N) g[k] = 2.0*0.5*hx*ht*(f0[k]-F(i, j)) + psi[j][i]*G0;
            }
            if (j==M)
            {
                if (i==0) g[k] = 2.0*0.25*hx*ht*(f0[k]-F(i, j));
                if (i==N) g[k] = 2.0*0.25*hx*ht*(f0[k]-F(i, j));
                if (0<i && i<N) g[k] = 2.0*0.5*hx*ht*(f0[k]-F(i, j)) + psi[j][i]*G0;
            }
            if (0<j && j<M)
            {
                if (i==0) g[k] = 2.0*0.5*hx*ht*(f0[k]-F(i, j));
                if (i==N) g[k] = 2.0*0.5*hx*ht*(f0[k]-F(i, j));
                if (0<i && i<N) g[k] = 2.0*1.0*hx*ht*(f0[k]-F(i, j)) + psi[j][i]*G0;
            }
        }
    }
    puts("---");
    Printer::printAsMatrix(g, M, N);
}

void Discretehyperbolic1::print(unsigned int iteration, const DoubleVector &x, const DoubleVector &gradient, double alpha, RnFunction *fn) const
{
    printf("J[%d]: %.12f\n", iteration, fn->fx(x));
}

double Discretehyperbolic1::fi1(unsigned int i) const
{
    double x = i*hx;
    return x*x*x;
}

double Discretehyperbolic1::fi2(unsigned int i) const
{
    return 0.0;
}

double Discretehyperbolic1::m1(unsigned int j) const
{
    double t = j*ht;
    return t*t*t;
}

double Discretehyperbolic1::m2(unsigned int j) const
{
    double t = j*ht;
    return t*t*t+1.0;
}

double Discretehyperbolic1::f(unsigned int i, unsigned int j) const
{
    return (*pf)[j*(N+1)+i];
//    double x = i*hx;
//    double t = j*ht;
//    return 6.0*t - 6.0*x*a;
}

double Discretehyperbolic1::F(unsigned int i, unsigned int j) const
{
    double x = i*hx;
    double t = j*ht;
    return 6.0*t - 6.0*x*a;
}

void Discretehyperbolic1::calculateP(const DoubleVector &u, DoubleMatrix &psi)
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

//        double G0 = -ht*ht;

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
                da[i-1] = A2;
                db[i-1] = B0;
                dc[i-1] = A1;
                rd[i-1] = -2.0*hx*1.0*(u[i]-U[i]);
            }
            da[0]=0.0;
            dc[N-2]=0.0;
            TomasAlgorithm(da, db, dc, rd, rx);
            for (unsigned int i=1; i<=N-1; i++)
            {
                psi[M][i] = rd[i-1];
            }
            psi[M][0]   = -(A1*psi[M][1]   + 2.0*hx*0.5*(u[0]-U[0]));
            psi[M][N]   = -(A2*psi[M][N-1] + 2.0*hx*0.5*(u[N]-U[N]));
            //Printer::printVector(psi[M]);
        }
        else if (j==(M-1))
        {
            for (unsigned int i=1; i<=N-1; i++)
            {
                da[i-1] = A2;
                db[i-1] = B0;
                dc[i-1] = A1;
                rd[i-1] = 0.0;

                if (i==1)
                {
                    rd[i-1] = -(D0*psi[M][1] + C1*psi[M][2]);
                }
                else if (i==N-1)
                {
                    rd[i-1] = -(C2*psi[M][N-2] + D0*psi[M][N-1]);
                }
                else
                {
                    rd[i-1] = -(C2*psi[M][i-1] + D0*psi[M][i] + C1*psi[M][i+1]);
                }
            }
            da[0]=0.0;
            dc[N-2]=0.0;
            TomasAlgorithm(da, db, dc, rd, rx);
            for (unsigned int i=1; i<=N-1; i++)
            {
                psi[M-1][i] = rd[i-1];
            }
            psi[M-1][0] = -(A1*psi[M-1][1]  +C1*psi[M][1]);
            psi[M-1][N] = -(A2*psi[M-1][N-1]+C2*psi[M][N-1]);
            //Printer::printVector(psi[j]);
        }
        else
        {
            if (j>=2)
            {
                for (unsigned int i=1; i<=N-1; i++)
                {
                    da[i-1] = A2;
                    db[i-1] = B0;
                    dc[i-1] = A1;
                    rd[i-1] = 0.0;

                    if (i==1)
                    {
                        rd[i-1] = -(D0*psi[j+1][1] + C1*psi[j+1][2] + F0*psi[j+2][1] + E1*psi[j+2][2]);
                    }
                    else if (i==N-1)
                    {
                        rd[i-1] = -(C2*psi[j+1][N-2] + D0*psi[j+1][N-1] + E2*psi[j+2][N-2] + F0*psi[j+2][N-1]);
                    }
                    else
                    {
                        rd[i-1] = -(C2*psi[j+1][i-1] + D0*psi[j+1][i] + C1*psi[j+1][i+1] + E2*psi[j+2][i-1] + F0*psi[j+2][i] + E1*psi[j+2][i+1]);
                    }

                    ////                if (j>=M)
                    ////                {
                    ////                    rd[i-1] -= 2.0*hx*ht*(u[j][i]-U);
                    ////                }

                    //                if (j==1)
                    //                {
                    //                    rd[i-1] += -psi[j][i];
                    //                }
                    //                if (j==0)
                    //                {
                    //                    rd[i-1] += (-psi[j][i] + psi[j+1][i]);
                    //                }
                }

                da[0]=0.0;
                dc[N-2]=0.0;

                TomasAlgorithm(da, db, dc, rd, rx);

                for (unsigned int i=1; i<=N-1; i++)
                {
                    psi[j][i] = rd[i-1];
                }

                psi[j][0] = -(A1*psi[j][1]   + C1*psi[j+1][1]   + E1*psi[j+2][1]);
                psi[j][N] = -(A2*psi[j][N-1] + C2*psi[j+1][N-1] + E2*psi[j+2][N-1]);

                ////            if (j>=M)
                ////            {
                ////                psi[j][0] -= 2.0*hx*ht*0.5*(u[j][0]-U);
                ////                psi[j][N] -= 2.0*hx*ht*0.5*(u[j][N]-U);
                ////            }

                //            if (j==1)
                //            {
                //                 psi[j][0] = -0.5*(A1*psi[j][1]   + C1*psi[j+1][1]   + E1*psi[j+2][1]);
                //                 psi[j][N] = -0.5*(A1*psi[j][N-1] + C1*psi[j+1][N-1] + E1*psi[j+2][N-1]);
                //            }
                //            if (j==0)
                //            {
                //                psi[j][0] = -0.5*(A1*psi[j][1]   + C1*psi[j+1][1]   + E1*psi[j+2][1] - psi[j+1][0]);
                //                psi[j][N] = -0.5*(A1*psi[j][N-1] + C1*psi[j+1][N-1] + E1*psi[j+2][N-1] - psi[j+1][N]);
                //            }
            }
            if (j==1)
            {
                psi[1][0] = -(C1*psi[2][1]   + E1*psi[3][1]);
                for (unsigned int i=1; i<=N-1; i++)
                {
                    if (i==1)
                    {
                        psi[1][1] = -(D0*psi[j+1][1] + C1*psi[j+1][2] + F0*psi[j+2][1] + E1*psi[j+2][2]);
                    }
                    else if (i==N-1)
                    {
                        psi[1][N-1] = -(C2*psi[j+1][N-2] + D0*psi[j+1][N-1] + E2*psi[j+2][N-2] + F0*psi[j+2][N-1]);
                    }
                    else
                    {
                        psi[1][i] = -(C2*psi[j+1][i-1] + D0*psi[j+1][i] + C1*psi[j+1][i+1] + E2*psi[j+2][i-1] + F0*psi[j+2][i] + E1*psi[j+2][i+1]);
                    }
                }
                psi[1][N] = -(C2*psi[2][N-1] + E2*psi[3][N-1]);
            }
            if (j==0)
            {
                psi[0][0] = -(E1*psi[3][1]);
                for (unsigned int i=1; i<=N-1; i++)
                {
                    if (i==1)
                    {
                        psi[0][1] = -(F0*psi[j+2][1] + E1*psi[j+2][2]);
                    }
                    else if (i==N-1)
                    {
                        psi[0][N-1] = -(E2*psi[j+2][N-2] + F0*psi[j+2][N-1]);
                    }
                    else
                    {
                        psi[0][i] = -(E2*psi[j+2][i-1] + F0*psi[j+2][i] + E1*psi[j+2][i+1]);
                    }
                }
                psi[0][N] = -(E2*psi[3][N-1]);
            }
            //Printer::printVector(psi[j]);
        }

    }

    rx.clear();
    rd.clear();
    dc.clear();
    db.clear();
    da.clear();

    //    Printer::printMatrix(u, 10, 10, "u");
    //    puts("------------------------------------");
    //    Printer::printMatrix(psi, 10, 10, "psi");
    //    puts("------------------------------------");
    //    Printer::printVector(psi[M+D]);
    //    Printer::printVector(psi[M+D-1]);
    //    Printer::printVector(psi[M+D-2]);
    //    Printer::printVector(psi[M+D-3]);
    //    Printer::printVector(psi[M+D-4]);
    //    Printer::printVector(psi[M+D-5]);
}

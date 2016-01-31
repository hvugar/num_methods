#include "discretehyperbolic1.h"
#include <tomasmethod.h>
#include <gradient_cjt.h>

void DiscreteHyperbolic1::main()
{
    DiscreteHyperbolic1 dh;
    dh.fx(0.92);

//    double a = 0.8;
//    double b = 1.0;
//    double t = 0.0;
//    goldenSectionSearch(a, b, t, &dh, 0.00001);
//    printf("%f\n", t);
}

double DiscreteHyperbolic1::fx(double t)
{
    t0 = 0.0;
    x0 = 0.0;
    x1 = 1.0;
    t1 = t;

    hx = 0.01;
    N = 100;
    ht = 0.01;
    M = (unsigned int) round((t1-t0)/ht);
    L = 1.0;
    X1 = 40;

    D = 10;
    a = 1.0;
    lamda = 0.25;
    U = 0.0;

    //fprintf(file, "---\n");

    DoubleVector v((2+L)*(M+D+1));
    for (unsigned int j=0; j<=M+D; j++)
    {
        //double t = j*ht;
        v[0*(M+D+1)+j] = 1.0;//t*t;
        v[1*(M+D+1)+j] = 1.0;//t*t+1.0;
        v[2*(M+D+1)+j] = 1.0;
    }

    ConjugateGradient g;
    g.setGradient(this);
    g.setFunction(this);
    g.setEpsilon1(0.000001);
    g.setEpsilon2(0.000001);
    g.setEpsilon3(0.000001);
    g.setR1MinimizeEpsilon(1.0, 0.01);
    g.setPrinter(this);
    g.calculate(v);

    DoubleVector gr1(v.size());
    IGradient::Gradient(this, 0.001, v, gr1);

    DoubleVector gr2(v.size());
    gradient(v, gr2);
    gr2.L2Normalize();

    FILE* file = fopen("20160131.txt", "a");
    fprintf(file, "------------------------------------------------------------------------------------------------------------------------\n");
    fprintf(file, "T:%f hx:%f ht:%f M:%d N:%d x:%f X:%d J[v]:%.20f\n", t, hx, ht, M, N, x1, X1, fx(v));
    unsigned int div = (M+D+1);
    fprintf(file, "Controls. Count:%d\n", div);
    IPrinter::printVector(v, "v1: ", div, 0*(M+D+1), 0*(M+D+1)+(M+D), file);
    IPrinter::printVector(v, "v2: ", div, 1*(M+D+1), 1*(M+D+1)+(M+D), file);
    IPrinter::printVector(v, "v3: ", div, 2*(M+D+1), 2*(M+D+1)+(M+D), file);
    fprintf(file, "Numerical gradients. Count:%d\n", div);
    IPrinter::printVector(gr1, "gr1:", (M+D+1), 0*(M+D+1), 0*(M+D+1)+(M+D), file);
    IPrinter::printVector(gr1, "gr2:", (M+D+1), 1*(M+D+1), 1*(M+D+1)+(M+D), file);
    IPrinter::printVector(gr1, "gr3:", (M+D+1), 2*(M+D+1), 2*(M+D+1)+(M+D), file);
    fprintf(file, "Analytic gradient. Count:%d\n", div);
    IPrinter::printVector(gr2, "gr1:", (M+D+1), 0*(M+D+1), 0*(M+D+1)+(M+D), file);
    IPrinter::printVector(gr2, "gr2:", (M+D+1), 1*(M+D+1), 1*(M+D+1)+(M+D), file);
    IPrinter::printVector(gr2, "gr3:", (M+D+1), 2*(M+D+1), 2*(M+D+1)+(M+D), file);
    fputs("Amplitudes:\n", file);
    DoubleMatrix u;
    pv = &v;
    IHyperbolicEquation::calculateU(u, hx, ht, M+D, N);
    for (unsigned int j=M; j<=M+D; j++)
    {
        char buffer[20];
        int n = sprintf(buffer, "u[%d]: ", j);
        buffer[n] = 0;
        IPrinter::printVector(u[j], buffer, u[j].size(), 0, 0, file);
    }
    fputs("------------------------------------------------------------------------------------------------------------------------\n", file);
    fclose(file);

    return fx(v);
}

double DiscreteHyperbolic1::fx(const DoubleVector& v)
{
    pv = &v;

    double sum = 0.0;

    DoubleMatrix u;
    IHyperbolicEquation::calculateU(u, hx, ht, M+D, N, a);

    for (unsigned int j=M; j<=M+D; j++)
    {
        for (unsigned int i=0; i<=N; i++)
        {
            double alpha = 1.0;
            if (i==0 || i==N || j==M+D) alpha = 0.5;
            if (i==0 && j==M+D) alpha = 0.25;
            if (i==N && j==M+D) alpha = 0.25;
            sum += alpha*(u[j][i]-U)*(u[j][i]-U);
        }
    }
    sum = hx*ht*sum;

    return sum;
}

void DiscreteHyperbolic1::gradient(const DoubleVector &v, DoubleVector &g)
{
    pv = &v;
    DoubleMatrix u;
    calculateU(u, hx, ht, M+D, N, a, lamda);
    DoubleMatrix psi;
    calculateP(v, u, psi, g);

    for (unsigned int j=2; j<=M+D; j++)
    {
        g[0*(M+D+1)+j] = -psi[j][0];
        g[1*(M+D+1)+j] = -psi[j][N];

        g[2*(M+D+1)+j] = 0.0;
        for (unsigned int i=X1; i<=N; i++)
        {
            g[2*(M+D+1)+j] += psi[j][i];
        }
        g[2*(M+D+1)+j] *= -(ht*ht);
    }

    g[0*(M+D+1)+0] = 0.0;
    g[0*(M+D+1)+1] = 0.0;

    g[1*(M+D+1)+0] = 0.0;
    g[1*(M+D+1)+1] = 0.0;

    g[2*(M+D+1)+0] = 0.0;
    g[2*(M+D+1)+1] = 0.0;
}

void DiscreteHyperbolic1::calculateP(const DoubleVector& f0, const DoubleMatrix &u, DoubleMatrix &psi, DoubleVector &g)
{
    for (unsigned int j=0; j<psi.size(); j++) psi[j].clear();
    psi.clear();

    psi.resize(M+D+1);
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

    for (unsigned int j1=0; j1<=M+D; j1++)
    {
        unsigned int j = M+D-j1;

        if (j==M+D)
        {
            for (unsigned int i=1; i<=N-1; i++)
            {
                //double alpha = 1.0;
                //if (j==M+D) alpha = 0.50;

                da[i-1] = A1;
                db[i-1] = B0;
                dc[i-1] = A2;
                rd[i-1] = -2.0*hx*ht*0.50*(u[M+D][i]-U);
            }
            da[0]=0.0;
            dc[N-2]=0.0;
            TomasAlgorithm(da, db, dc, rd, rx);
            for (unsigned int i=1; i<=N-1; i++)
            {
                psi[j][i] = rx[i-1];
            }

            //double alpha = 0.50;
            //if (j==M+D) alpha = 0.25;

            psi[M+D][0] = -(A2*psi[M+D][1]   + 2.0*hx*ht*0.25*(u[M+D][0]-U));
            psi[M+D][N] = -(A1*psi[M+D][N-1] + 2.0*hx*ht*0.25*(u[M+D][N]-U));
        }
        else if (j==(M+D-1))
        {
            for (unsigned int i=1; i<=N-1; i++)
            {
                da[i-1] = A1;
                db[i-1] = B0;
                dc[i-1] = A2;
                rd[i-1] = 0.0;

                if (i==1)
                {
                    rd[i-1] = -(D0*psi[M+D][1] + C2*psi[M+D][2]);
                }
                else if (i==N-1)
                {
                    rd[i-1] = -(C1*psi[M+D][N-2] + D0*psi[M+D][N-1]);
                }
                else
                {
                    rd[i-1] = -(C1*psi[M+D][i-1] + D0*psi[M+D][i] + C2*psi[M+D][i+1]);
                }

                rd[i-1] -= 2.0*hx*ht*1.0*(u[j][i]-U);
            }
            da[0]=0.0;
            dc[N-2]=0.0;
            TomasAlgorithm(da, db, dc, rd, rx);
            for (unsigned int i=1; i<=N-1; i++)
            {
                psi[M+D-1][i] = rx[i-1];
            }
            psi[M+D-1][0] = -(A2*psi[M+D-1][1]  +C2*psi[M+D][1]  + 2.0*hx*ht*0.5*(u[j][0]-U));
            psi[M+D-1][N] = -(A1*psi[M+D-1][N-1]+C1*psi[M+D][N-1]+ 2.0*hx*ht*0.5*(u[j][N]-U));
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

                    if (j>=M)
                    {
                        rd[i-1] -= 2.0*hx*ht*1.0*(u[j][i]-U);
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
                psi[j][N] = -(A1*psi[j][N-1] + C1*psi[j+1][N-1] + E1*psi[j+2][N-1]);

                if (j>=M)
                {
                    psi[j][0] -= 2.0*hx*ht*0.5*(u[j][0]-U);
                    psi[j][N] -= 2.0*hx*ht*0.5*(u[j][N]-U);
                }
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
}

void DiscreteHyperbolic1::print(unsigned int iteration, const DoubleVector &v, const DoubleVector &g, double alpha, RnFunction *fn) const
{
    printf("J[%d]: %.16f\n", iteration, fn->fx(v));
}

double DiscreteHyperbolic1::fi1(unsigned int i) const
{
    return 2.0;
}

double DiscreteHyperbolic1::fi2(unsigned int i) const
{
    return 0.0;
}

double DiscreteHyperbolic1::m1(unsigned int j) const
{
    const DoubleVector &v = *pv;
    double v1 = v[0*(M+D+1)+j];
    return v1;
}

double DiscreteHyperbolic1::m2(unsigned int j) const
{
    const DoubleVector &v = *pv;
    double v2 = v[1*(M+D+1)+j];
    return v2;
}

double DiscreteHyperbolic1::f(unsigned int i, unsigned int j) const
{
    const DoubleVector &v = *pv;
    double v3 = v[2*(M+D+1)+j];

    double sum = 0.0;
    if (i>=X1)
    {
        sum = v3;
    }
    return sum;
}

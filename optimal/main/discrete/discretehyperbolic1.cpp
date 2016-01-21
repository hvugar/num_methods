#include "discretehyperbolic1.h"
#include <tomasmethod.h>
#include <gradient_cjt.h>


FILE *file;
void DiscreteHyperbolic1::main()
{
    file = fopen("20160120.txt", "w");
    //file = stdout;
    DiscreteHyperbolic1 dh;
//    dh.fx(0.996462);

    double a = 0.8;
    double b = 1.0;
    double t = 0.0;
    goldenSectionSearch(a, b, t, &dh, 0.00001);
    printf("%f\n", t);

    fclose(file);
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

//    v[0] = 0.0;
//    v[1] = 0.0;//dh.ht*dh.ht;
//    v[(M+D+1)] = 0.0;//1.0;
//    v[(M+D+1)] = 0.0;//dh.ht*dh.ht+1.0;
    printf("t: %f %d\n", t, M+D);
    fprintf(file, "%f %d\n", t, M+D);
    fflush(file);

    ConjugateGradient g;
    g.setGradient(this);
    g.setFunction(this);
    g.setEpsilon1(0.000001);
    g.setEpsilon2(0.000001);
    g.setEpsilon3(0.000001);
    g.setR1MinimizeEpsilon(1.0, 0.01);
    g.setPrinter(this);
    g.calculate(v);

    DoubleVector gr(v.size());
    for (unsigned int i=0; i<v.size(); i++)
    {
        DoubleVector v1 = v;
        DoubleVector v2 = v;
        double h = 0.01;
        v1[i] = v[i] + h;
        v2[i] = v[i] - h;
        gr[i] = (fx(v1)-fx(v2))/(2.0*h);
    }

    fprintf(file, "----\n");

    IPrinter::printVector(v, "v1:", M+D+1, 0*(M+D+1), 0*(M+D+1)+(M+D), file);
    IPrinter::printVector(v, "v2:", M+D+1, 1*(M+D+1), 1*(M+D+1)+(M+D), file);
    IPrinter::printVector(v, "v3:", M+D+1, 2*(M+D+1), 2*(M+D+1)+(M+D), file);

    IPrinter::printVector(gr, "g1:", M+D+1, 0*(M+D+1), 0*(M+D+1)+(M+D), file);
    IPrinter::printVector(gr, "g2:", M+D+1, 1*(M+D+1), 1*(M+D+1)+(M+D), file);
    IPrinter::printVector(gr, "g3:", M+D+1, 2*(M+D+1), 2*(M+D+1)+(M+D), file);

//    Printer::printVector(v, "v1:", D+1, 0, M+D, file);
//    Printer::printVector(v, "v2:", D+1, M+D+1, (2*(M+D)+1), file);

//    Printer::printVector(v, "v1:", D+1, 0*(M+D+1)+M, 0*(M+D+1)+(M+D), file);
//    Printer::printVector(v, "v2:", D+1, 1*(M+D+1)+M, 1*(M+D+1)+(M+D), file);
//    Printer::printVector(v, "v3:", D+1, 2*(M+D+1)+M, 2*(M+D+1)+(M+D), file);

//    DoubleMatrix u;
//    pv = &v;
//    calculateU(u, hx, ht, M+D, N);
//    puts("-----");
//    for (unsigned int j=M; j<=M+D; j++)
//    {
//        char buffer[20];
//        int n = sprintf(buffer, "u[%d]: ", j);
//        buffer[n] = 0;
//        Printer::printVector(u[j], buffer, u[j].size(), 0, 0, file);
//    }

    return fx(v);
}

double DiscreteHyperbolic1::fx(const DoubleVector& v)
{
    pv = &v;

    double sum = 0.0;

    DoubleMatrix u;
    calculateU(u, hx, ht, M+D, N, a);

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
    for (unsigned int i=0; i<v.size(); i++)
    {
        DoubleVector v1 = v;
        DoubleVector v2 = v;
        double h = 0.01;
        v1[i] = v[i] + h;
        v2[i] = v[i] - h;
        g[i] = (fx(v1)-fx(v2))/(2.0*h);
    }

//    Printer::printVector(v, "v1:", 10, 0*(M+D+1), 0*(M+D+1)+(M+D), file);
//    Printer::printVector(v, "v2:", 10, 1*(M+D+1), 1*(M+D+1)+(M+D), file);
//    Printer::printVector(v, "v3:", 10, 2*(M+D+1), 2*(M+D+1)+(M+D), file);

//    Printer::printVector(g, "g1:", 10, 0*(M+D+1), 0*(M+D+1)+(M+D), file);
//    Printer::printVector(g, "g2:", 10, 1*(M+D+1), 1*(M+D+1)+(M+D), file);
//    Printer::printVector(g, "g3:", 10, 2*(M+D+1), 2*(M+D+1)+(M+D), file);
}
/*
void DiscreteHyperbolic1::gradient(const DoubleVector &v, DoubleVector &g)
{
    pv = &v;
    DoubleMatrix u;
    calculateU(u, hx, ht, M+D, N, a, lamda);
    DoubleMatrix psi;
    calculateP(v, u, psi, g);

    //double sgm = 3.0*hx;
    //double a = 1.0/(sgm*sqrt(2.0*M_PI));
    //double b = 2.0*sgm*sgm;

    for (unsigned int j=2; j<=M+D; j++)
    {
        g[0*(M+D+1)+j] = -psi[j][0];
        g[1*(M+D+1)+j] = -psi[j][N];
        g[2*(M+D+1)+j] = (ht*ht)*(1/hx)*psi[j][X1];

//        g[2*(M+D+1)+j] = 0.0;
//        for (unsigned int i=0; i<=N; i++)
//        {
//            double x = i*hx;
//            double c = a * exp(-((x-0.4)*(x-0.4))/b);

//            g[2*(M+D+1)+j] += psi[j][i] * c / hx;
//        }
        //printf("%.12f %.12f %.12f %.12f\n", g[2*(M+D+1)+j], g[1*(M+D+1)+j], g[1*(M+D+1)+j], psi[j][X1]);
    }

    g[0*(M+D+1)+0] = 0.0;
    g[0*(M+D+1)+1] = 0.0;

    g[1*(M+D+1)+0] = 0.0;
    g[1*(M+D+1)+1] = 0.0;

    g[2*(M+D+1)+0] = 0.0;
    g[2*(M+D+1)+1] = 0.0;

    Printer::printVector(v, "v1:", D+1, 0*(M+D+1)+M, 0*(M+D+1)+(M+D), file);
    Printer::printVector(v, "v2:", D+1, 1*(M+D+1)+M, 1*(M+D+1)+(M+D), file);
    Printer::printVector(v, "v3:", D+1, 2*(M+D+1)+M, 2*(M+D+1)+(M+D), file);

    Printer::printVector(g, "g1:", D+1, 0*(M+D+1)+M, 0*(M+D+1)+(M+D), file);
    Printer::printVector(g, "g2:", D+1, 1*(M+D+1)+M, 1*(M+D+1)+(M+D), file);
    Printer::printVector(g, "g3:", D+1, 2*(M+D+1)+M, 2*(M+D+1)+(M+D), file);

    //fprintf(file, "psi1: "); for (unsigned int j=M; j<=M+D; j++) fprintf(file, psi[j][0]>0?"+%.10f ":"%.10f ", psi[j][0]); fprintf(file, "\n");
    //fprintf(file, "psi2: "); for (unsigned int j=M; j<=M+D; j++) fprintf(file, psi[j][N]>0?"+%.10f ":"%.10f ", psi[j][N]); fprintf(file, "\n");
    //fprintf(file, "psi3: "); for (unsigned int j=M; j<=M+D; j++) fprintf(file, psi[j][X1]>0?"+%.10f ":"%.10f ", psi[j][X1]); fprintf(file, "\n");

    //Printer::printMatrix(psi,10,11,NULL,file);
    for (unsigned int j=M; j<=M+D; j++)
    {
        char buffer[20];
        int n = sprintf(buffer, "psi[%d]: ", j);
        buffer[n] = 0;
        Printer::printVector(psi[j], buffer, 10, 0, 0, file);
    }


}
*/
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
    //double G0 = -ht*ht;

//    A1 /= G0;
//    B0 /= G0;
//    A2 /= G0;
//    C1 /= G0;
//    D0 /= G0;
//    C2 /= G0;
//    E1 /= G0;
//    F0 /= G0;
//    E2 /= G0;
//    G0 = 1.0;

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

//    DoubleMatrix u;
//    calculateU(u, hx, ht, M+D, N, a);
//    for (unsigned int j=M; j<=M+D; j++)
//    {
//        char buffer[20];
//        int n = sprintf(buffer, "u[%d]:\t", j);
//        buffer[n] = 0;
//        Printer::printVector(u[j], buffer, 10, 0, 0, file);
//    }
//    Printer::printVector(g, "g1:\t", 11, 0, M+D, file);
//    Printer::printVector(g, "g2:\t", 11, M+D+1, (2*(M+D)+1), file);
}

double DiscreteHyperbolic1::fi1(unsigned int i) const
{
    double x = i*hx;
    return x*x;
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
    if (i==X1)
    {
        sum = (1.0/hx)*v3;
    }
    return sum;

//    double x = i*hx;
//    double sum = 0.0;

//    double sgm = 6.0*hx;
//    double a = 1.0/(sgm*sqrt(2.0*M_PI));
//    double b = 2.0*sgm*sgm;
//    double g = a * exp(-((x-0.4)*(x-0.4))/b);
//    sum += v3 * g / hx;
//    return sum;


    //return 2.0-2.0*a*a;
}

#include "discretehyperbolic.h"
#include <tomasmethod.h>
#include <gradient_cjt.h>

void DiscreteHyperbolic::main()
{
    DiscreteHyperbolic hc;
    hc.calculateSettings();
    hc.e[0] = 0.4;
    hc.count = 0;
    hc.fx(1.0);
    printf("Function call count: %u\n", hc.count);
}

DiscreteHyperbolic::DiscreteHyperbolic()
{

}

DiscreteHyperbolic::~DiscreteHyperbolic()
{

}

void DiscreteHyperbolic::calculateSettings()
{
    t0 = 0.0;

    x0 = 0.0;
    x1 = 1.0;

    N = 100;
    hx = (x1-x0)/N;

    ht = 0.01;
    M = (unsigned int) round((t1-t0)/ht);
    D = 10;

    a = 1.0;
    lamda = 0.25;
    R = 1.0;
    U = 0.0;

    L = 1;
    e.resize(L);
}

double DiscreteHyperbolic::fx(const DoubleVector& x)
{
    return 0.0;
}

void DiscreteHyperbolic::gradient(const DoubleVector& x, DoubleVector& g, double)
{

}

double DiscreteHyperbolic::fx(double t)
{
    t1 = t;
    calculateSettings();

    DoubleVector v((M+D+1)*(L+2));
    for (unsigned int j=0; j<=(M+D); j++)
    {
        //double t = j*ht;
        v[0*(M+D+1)+j] = 0.0;
        v[1*(M+D+1)+j] = 0.0;
        v[2*(M+D+1)+j] = 0.0;
    }

    DoubleVector g;
    DoubleMatrix u1;
    calculateU(v, u1);
    calculateP(u1, g);
    return 0.0;

    double min_step = 2.0;
    double gold_eps = 0.0001;

    ConjugateGradient cg;
    cg.setFunction(this);
    cg.setEpsilon1(0.000001);
    cg.setEpsilon2(0.000001);
    //g.setGradientStep(0.000001);
    //    g.setR1MinimizeEpsilon(20.0, 0.01);
    cg.setR1MinimizeEpsilon(min_step, gold_eps);
    cg.setPrinter(this);
    cg.setProjection(this);
    cg.setNormalize(true);
    cg.showEndMessage(false);
    cg.calculate(v);

    double rf = fx(v);

    //    DoubleVector v1(D+1); for (unsigned j=M; j<=M+D; j++) v1[j] = v[0*(M+D+1)+j];
    //    DoubleVector v2(D+1); for (unsigned j=M; j<=M+D; j++) v2[j] = v[1*(M+D+1)+j];
    //    DoubleVector v3(D+1); for (unsigned j=M; j<=M+D; j++) v3[j] = v[2*(M+D+1)+j];
    //    Printer::printVector(v1, 10, "v1:\t");
    //    Printer::printVector(v2, 10, "v2:\t");
    //    Printer::printVector(v3, 10, "v2:\t");

    DoubleMatrix u;
    DoubleVector gr(3*(M+D+1));
    calculateU(v, u);
    pu = &u;
    calculateP(u, gr);

    printf("Norm: %.12f %d\n", gr.L2Norm(), cg.count());

    FILE* file = fopen("20151220.txt", "a");
    fprintf(file, "------------------------------------------------------------\n");
    fprintf(file, "e1: %f T: %.8f Functional: %.16f hx: %f ht: %f step: %f gold_epsilon: %f N: %d M: %d\n", e[0], t, rf, hx, ht, min_step, gold_eps, N, M);

    //    Printer::printVector(v, "v1:\t", M+D+1, 0, M+D, file);
    //    Printer::printVector(v, "v2:\t", M+D+1, (M+D+1), (M+D+1)+M+D+1, file);
    //    Printer::printVector(v, "v3:\t", M+D+1, 2*(M+D+1), 2*(M+D+1)+M+D+1, file);
    //    Printer::printVector(g, "g1:\t", M+D+1, 0, M+D+1, file);
    //    Printer::printVector(g, "g2:\t", M+D+1, (M+D+1), (M+D+1)+M+D+1, file);
    //    Printer::printVector(g, "g3:\t", M+D+1, 2*(M+D+1), 2*(M+D+1)+M+D+1, file);

    fprintf(file, "v1:\t");
    for (unsigned int j=0; j<=M+D; j++)
    {
        double v1 = v[0*(M+D+1)+j];
        if (v1<0)
            fprintf(file, "%14.8f ", v1);
        else
            fprintf(file, "%+14.8f ", v1);
    }
    fprintf(file, "\n");

    fprintf(file, "v2:\t");
    for (unsigned int j=0; j<=M+D; j++)
    {
        double v2 = v[1*(M+D+1)+j];
        if (v2<0)
            fprintf(file, "%14.8f ", v2);
        else
            fprintf(file, "%+14.8f ", v2);
    }
    fprintf(file, "\n");

    fprintf(file, "v3:\t");
    for (unsigned int j=0; j<=M+D; j++)
    {
        double v3 = v[2*(M+D+1)+j];
        if (v3<0)
            fprintf(file, "%14.8f ", v3);
        else
            fprintf(file, "%+14.8f ", v3);
    }
    fprintf(file, "\n");

    fprintf(file, "g1:\t");
    for (unsigned int j=0; j<=M+D; j++)
    {
        double g1 = gr[0*(M+D+1)+j];
        if (g1<0)
            fprintf(file, "%14.8f ", g1);
        else
            fprintf(file, "%+14.8f ", g1);
    }
    fprintf(file, "\n");

    fprintf(file, "g2:\t");
    for (unsigned int j=0; j<=M+D; j++)
    {
        double g2 = gr[1*(M+D+1)+j];
        if (g2<0)
            fprintf(file, "%14.8f ", g2);
        else
            fprintf(file, "%+14.8f ", g2);
    }
    fprintf(file, "\n");

    fprintf(file, "g3:\t");
    for (unsigned int j=0; j<=M+D; j++)
    {
        double g3 = gr[2*(M+D+1)+j];
        if (g3<0)
            fprintf(file, "%14.8f ", g3);
        else
            fprintf(file, "%+14.8f ", g3);
    }
    fprintf(file, "\n");

    for (unsigned int j=M; j<=M+D; j++)
    {
        //printf("u[%d]:\t", j);
        //Printer::printVector(u[j]);

        fprintf(file, "u[%d]:\t", j);
        for (unsigned int i=0; i<=N; i++)
        {
            double uji = u[j][i];
            if (uji<0)
                fprintf(file, "%10.8f ", uji);
            else
                fprintf(file, "+%10.8f ", uji);
        }
        fprintf(file, "\n");
    }
    fclose(file);
    //Printer::printVector(u[M], 10, "u[M]:\t");
    //Printer::printVector(u[M+D], 10, "u[M+D]:");

    printf("e1: %f T: %.8f Integral: %.16f M: %d\n", e[0], t, rf, M);

    return rf;
}

void DiscreteHyperbolic::project(DoubleVector &x, int index)
{

}

void DiscreteHyperbolic::print(unsigned int iteration, const DoubleVector& x, const DoubleVector &gradient, double alpha, RnFunction* fn) const
{

}

void DiscreteHyperbolic::calculateU(const DoubleVector& v, DoubleMatrix& u)
{
    pv = &v;

    u.clear();
    u.resize(M+D+1);
    for (unsigned int j=0; j<=M+D; j++) u[j].resize(N+1);

    DoubleVector u0(N+1);
    DoubleVector u1(N+1);

    DoubleVector da;
    DoubleVector db;
    DoubleVector dc;
    DoubleVector rd;
    DoubleVector rx;

    da.resize(N-1);
    db.resize(N-1);
    dc.resize(N-1);
    rd.resize(N-1);
    rx.resize(N-1);

    double alpha1 = -(lamda*a*a)*((ht*ht)/(hx*hx));
    double beta1  = 1.0 + (2.0*lamda*a*a*(ht*ht))/(hx*hx);
    double alpha2 = (1.0-2.0*lamda)*a*a*((ht*ht)/(hx*hx));
    double alpha3 = +(lamda*a*a)*((ht*ht)/(hx*hx));

    for (unsigned int j=0; j<=M+D-1; j++)
    {
        if (j==0)
        {
            for (unsigned int i=0; i<=N; i++)
            {
                u0[i] = fi1(i);
                u1[i] = u0[i] + ht*fi2(i);
                u[0][i] = u0[i];
                u[1][i] = u1[i];
            }
        }
        else
        {
            for (unsigned int i=1; i<=N-1; i++)
            {
                da[i-1] = alpha1;
                db[i-1] = beta1;
                dc[i-1] = alpha1;
                rd[i-1] = alpha2*(u1[i-1]-2.0*u1[i]+u1[i+1]) + 2.0*u1[i] + alpha3*(u0[i+1] - 2.0*u0[i] + u0[i-1]) - u0[i] + (ht*ht)*f(i, j);
            }

            da[0]   = 0.0;
            dc[N-2] = 0.0;
            rd[0]   -= alpha1 * m1(j+1);
            rd[N-2] -= alpha1 * m2(j+1);
            TomasAlgorithm(da, db, dc, rd, rx);

            u[j+1][0] = m1(j+1);
            for (unsigned int i=1; i<=N-1; i++)
            {
                u[j+1][i] = rx[i-1];
            }
            u[j+1][N] = m2(j+1);

            for (unsigned int i=0; i<=N; i++)
            {
                u0[i] = u1[i];
                u1[i] = u[j+1][i];
            }
        }
    }

    da.clear();
    db.clear();
    dc.clear();
    rd.clear();
    rx.clear();

    u1.clear();
    u0.clear();
}

void DiscreteHyperbolic::calculateP(const DoubleMatrix &u, DoubleVector &g)
{
    DoubleMatrix psi;
    psi.clear();
    psi.resize(M+D+1);
    for (unsigned int j=0; j<=M+D; j++) psi[j].resize(N+1);

    double A1 = -(lamda*a*a*ht*ht)/(hx*hx);
    double B1 = 1.0 + (2.0*lamda*a*a*ht*ht)/(hx*hx);
    double C1 = -(1.0-2.0*lamda)*(a*a*ht*ht)/(hx*hx);
    double D1 = 2.0*(((1.0-2.0*lamda)*(a*a*ht*ht) - hx*hx)/(hx*hx));
    double E1 = 1.0 + (2.0*lamda*a*a*ht*ht)/(hx*hx);
    double F1 = -ht*ht;

    DoubleVector da(N-1);
    DoubleVector db(N-1);
    DoubleVector dc(N-1);
    DoubleVector rd(N-1);
    DoubleVector rx(N-1);

    for (unsigned int j1=0; j1<=M+D; j1++)
    {
        unsigned int j = M+D-j1;

        if (j==(M+D))
        {
            for (unsigned int i=1; i<=N-1; i++)
            {
                da[i-1] = A1;
                db[i-1] = B1;
                dc[i-1] = A1;
                rd[i-1] = -2.0*hx*ht*0.5*(u[j][i]-U);
            }
            da[0]=0.0;
            dc[N-2]=0.0;
            TomasAlgorithm(da, db, dc, rd, rx);
            for (unsigned int i=1; i<=N-1; i++)
            {
                psi[j][i] = rd[i-1];
            }
            psi[j][0]   = -(A1*psi[j][1]  +2.0*hx*ht*0.25*(u[j][0]-U));
            psi[j][N]   = -(A1*psi[j][N-1]+2.0*hx*ht*0.25*(u[j][N]-U));
            Printer::printVector(psi[j]);
        }
        else if (j==(M+D-1))
        {
            for (unsigned int i=1; i<=N-1; i++)
            {
                da[i-1] = A1;
                db[i-1] = B1;
                dc[i-1] = A1;
                rd[i-1] = 0.0;

                if (i==1)
                {
                    rd[i-1] = -(D1*psi[j+1][i] + C1*psi[j+1][i+1]+2.0*hx*ht*(u[j][i]-U));
                }
                else if (i==N-1)
                {
                    rd[i-1] = -(C1*psi[j+1][i-1] + D1*psi[j+1][i]+2.0*hx*ht*(u[j][i]-U));
                }
                else
                {
                    rd[i-1] = -(C1*psi[j+1][i-1] + D1*psi[j+1][i] + C1*psi[j+1][i+1]+2.0*hx*ht*(u[j][i]-U));
                }
            }
            da[0]=0.0;
            dc[N-2]=0.0;
            TomasAlgorithm(da, db, dc, rd, rx);
            for (unsigned int i=1; i<=N-1; i++)
            {
                psi[j][i] = rd[i-1];
            }
            psi[j][0] = -(A1*psi[j][1]  +C1*psi[j+1][1]  +2.0*hx*ht*0.5*(u[j][0]-U));
            psi[j][N] = -(A1*psi[j][N-1]+C1*psi[j+1][N-1]+2.0*hx*ht*0.5*(u[j][N]-U));
            Printer::printVector(psi[j]);
        }
        else
        {
            for (unsigned int i=1; i<=N-1; i++)
            {
                da[i-1] = A1;
                db[i-1] = B1;
                dc[i-1] = A1;
                rd[i-1] = 0.0;

                if (i==1)
                {
                    rd[i-1] = -(D1*psi[j+1][i] + C1*psi[j+1][i+1] + E1*psi[j+2][i] + A1*psi[j+2][i+1]);
                }
                else if (i==N-1)
                {
                    rd[i-1] = -(C1*psi[j+1][i-1] + D1*psi[j+1][i] + A1*psi[j+2][i-1] + E1*psi[j+2][i]);
                }
                else
                {
                    rd[i-1] = -(C1*psi[j+1][i-1] + D1*psi[j+1][i] + C1*psi[j+1][i+1] + A1*psi[j+2][i-1] + E1*psi[j+2][i] + A1*psi[j+2][i+1]);
                }

//                if (j>=M)
//                {
//                    rd[i-1] -= 2.0*hx*ht*(u[j][i]-U);
//                }

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

            psi[j][0] = -(A1*psi[j][1]   + C1*psi[j+1][1]   + A1*psi[j+2][1]);
            psi[j][N] = -(A1*psi[j][N-1] + C1*psi[j+1][N-1] + A1*psi[j+2][N-1]);

//            if (j>=M)
//            {
//                psi[j][0] -= 2.0*hx*ht*0.5*(u[j][0]-U);
//                psi[j][N] -= 2.0*hx*ht*0.5*(u[j][N]-U);
//            }

//            if (j==1)
//            {
//                psi[j][0] = -0.5*(A1*psi[j][1]   + C1*psi[j+1][1]   + A1*psi[j+2][1]);
//                psi[j][N] = -0.5*(A1*psi[j][N-1] + C1*psi[j+1][N-1] + A1*psi[j+2][N-1]);
//            }
//            if (j==0)
//            {
//                psi[j][0] = -0.5*(A1*psi[j][1]   + C1*psi[j+1][1]   + A1*psi[j+2][1] - psi[j+1][0]);
//                psi[j][N] = -0.5*(A1*psi[j][N-1] + C1*psi[j+1][N-1] + A1*psi[j+2][N-1] - psi[j+1][N]);
//            }
            Printer::printVector(psi[j]);
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

double DiscreteHyperbolic::fi1(unsigned int i) const
{
    double x = i*hx;
    return x*x;
}

double DiscreteHyperbolic::fi2(unsigned int i) const
{
    return 0.0;
}

double DiscreteHyperbolic::m1(unsigned int j) const
{
    const DoubleVector &v = *pv;
    double v1 = v[j];
    return v1;
}

double DiscreteHyperbolic::m2(unsigned int j) const
{
    const DoubleVector &v = *pv;
    double v2 = v[M+D+1 + j];
    return v2;
}

double DiscreteHyperbolic::f(unsigned int i, unsigned int j) const
{
    double x = i*hx;
    //double t = j*ht;
    double sum = 0.0;

    const DoubleVector &v = *pv;
    double v3 = v[2*(M+D+1)+j];

    static double sgm = 6.0*hx;
    static double a = 1.0/(sgm*sqrt(2.0*M_PI));
    static double b = 2.0*sgm*sgm;
    double g = a * exp(-((x-e[0])*(x-e[0]))/b);
    sum += v3 * g * hx;

    return sum;
}



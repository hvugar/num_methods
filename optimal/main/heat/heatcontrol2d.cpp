#include "heatcontrol2d.h"
#include <tomasmethod.h>
#include <gradient_cjt.h>
#include <gradient_sd.h>

HeatControl2D::HeatControl2D() : RnFunction()
{
    N1 = 10;
    N2 = 10;
    M  = 1000;

    C =(M+1)*(N2+1)*(N1+1);

    t0 = 0.0;
    t1 = 1.0;
    x10 = 0.0;
    x11 = 1.0;
    x20 = 0.0;
    x21 = 1.0;

    a1 = a2 = 1.0;

    ht = (t1-t0)/M;
    h1 = (x11-x10)/N1;
    h2 = (x21-x20)/N2;

    U.resize(N2+1);
    for (unsigned int j=0; j<=N2; j++)
    {
        U[j].resize(N1+1);
        for (unsigned int i=0; i<=N1; i++) U[j][i] = u(i*h1, j*h2, t1);
    }
}

HeatControl2D::~HeatControl2D()
{
}

double HeatControl2D::fx(const DoubleVector &f)
{
    calculateU(f);

    double sum = 0.0;

    for (unsigned int j=0; j<N2; j++)
    {
        for (unsigned int i=0; i<N1; i++)
        {
            int j1 = j;
            int j2 = j+1;
            int i1 = i;
            int i2 = i+1;

            double f1 = uT[j1][i1] - U[j1][i1];
            double f2 = uT[j1][i2] - U[j1][i2];
            double f3 = uT[j2][i1] - U[j2][i1];
            double f4 = uT[j2][i2] - U[j2][i2];

            sum = sum + (0.25*(h1*h2))*(f1*f1 + f2*f2 + f3*f3 + f4*f4);
        }
    }

    double norm = 0.0;

    for (unsigned int i=0; i<f.size(); i++)
    {
        double x1 = ((i%((N1+1)*(N2+1)))%(N2+1)) * h1;
        double x2 = ((i%((N1+1)*(N2+1)))%(N1+1)) * h2;
        double t = (i/((N1+1)*(N2+1))) * ht;
        norm += (f[i] - fxt(x1, x2, t)) * (f[i] - fxt(x1, x2, t));
    }

    return sum + norm;
}

void HeatControl2D::gradient(double step, const DoubleVector &f, DoubleVector &g)
{
    calculateU(f);
    calculateP(f, g);
}

void HeatControl2D::calculateU(const DoubleVector &f)
{
    DoubleMatrix u0;
    DoubleMatrix u1;

    u0.resize(N2+1); for (unsigned int j=0; j<=N2; j++) u0[j].resize(N1+1);
    u1.resize(N2+1); for (unsigned int j=0; j<=N2; j++) u1[j].resize(N1+1);

    double alpha1 = 0;
    double beta1  = 0;
    double alpha2 = 0;
    double beta2  = 0;

    DoubleVector a;
    DoubleVector b;
    DoubleVector c;
    DoubleVector d;
    DoubleVector x;

    for (unsigned int k=0; k<=M; k++)
    {
        if (k==0)
        {
            for (unsigned int j=0; j<=N2; j++)
            {
                for (unsigned i=0; i<=N1; i++)
                {
                    u0[j][i] = fi(i*h1, j*h2);
                }
            }
        }
        else
        {
            // Approximation to x1 direction
            alpha1 = -(a2*ht)/(2.0*h2*h2);
            beta1  = 1.0 + (a2*ht)/(h2*h2);
            alpha2 = (a1*ht)/(2.0*h1*h1);
            beta2  = 1.0 - (a1*ht)/(h1*h1);

            a.resize(N1-1);
            b.resize(N1-1);
            c.resize(N1-1);
            d.resize(N1-1);
            x.resize(N1-1);

            for (unsigned int i=1; i<N1; i++)
            {
                for (unsigned int j=1; j<N2; j++)
                {
                    a[j-1] = alpha1;
                    b[j-1] = beta1;
                    c[j-1] = alpha1;
                    d[j-1] = alpha2*u0[j][i-1] + beta2*u0[j][i] + alpha2*u0[j][i+1] + (ht/2.0) * f[k*(N1+1)*(N2+1)+i*(N2+1)+j]/*fxt(i*h1, j*h2, (k-1.0)*ht)*/;
                }

                a[0]     = 0.0;
                c[N2-2]  = 0.0;
                d[0]    -= alpha1 * m3(h1*i, ht*(k));
                d[N2-2] -= alpha1 * m4(h1*i, ht*(k));

                TomasAlgorithm(a, b, c, d, x);

                u1[0][i]  = m3(h1*i, ht*(k-0.5));
                for (unsigned int j=1; j<N2; j++)
                {
                    u1[j][i] = x[j-1];
                }
                u1[N2][i] = m4(h1*i, ht*(k-0.5));
            }

            for (unsigned int j=0; j<=N2; j++)
            {
                u1[j][0]  = m1(h2*j, ht*(k-0.5));
                u1[j][N2] = m2(h2*j, ht*(k-0.5));
            }

            a.clear();
            b.clear();
            c.clear();
            d.clear();
            x.clear();

            // Approximation to x2 direction
            alpha1 = -(a1*ht)/(2.0*h1*h1);
            beta1  = 1.0 + (a1*ht)/(h1*h1);
            alpha2 = (a2*ht)/(2.0*h2*h2);
            beta2  = 1.0 - (a2*ht)/(h2*h2);

            a.resize(N2-1);
            b.resize(N2-1);
            c.resize(N2-1);
            d.resize(N2-1);
            x.resize(N2-1);

            for (unsigned int j=1; j<N2; j++)
            {
                for (unsigned int i=1; i<N1; i++)
                {
                    a[i-1] = alpha1;
                    b[i-1] = beta1;
                    c[i-1] = alpha1;
                    d[i-1] = alpha2*u1[j-1][i] + beta2*u1[j][i] + alpha2*u1[j+1][i] + (ht/2.0) * f[k*(N1+1)*(N2+1)+j*(N1+1)+i]/*fxt(i*h1, j*h2, (k-1)*ht)*/;
                }
                a[0]     = 0.0;
                c[N1-2]  = 0.0;
                d[0]    -= alpha1 * m1(h2*j, ht*(k));
                d[N1-2] -= alpha1 * m2(h2*j, ht*(k));
                TomasAlgorithm(a, b, c, d, x);

                u0[j][0]  = m1(h2*j, ht*(k));
                for (unsigned int i=1; i<N1; i++)
                {
                    u0[j][i] = x[i-1];
                }
                u0[j][N1] = m2(h2*j, ht*(k));
            }

            for (unsigned int i=0; i<=N1; i++)
            {
                u0[0][i]  = m3(h1*i, ht*(k));
                u0[N1][i] = m4(h1*i, ht*(k));
            }

            a.clear();
            b.clear();
            c.clear();
            d.clear();
            x.clear();
        }
    }

    uT = u0;

//    Printer::printMatrix(u0, N2/10, N1/10);
//    puts("---");
}

void HeatControl2D::calculateP(const DoubleVector &f, DoubleVector& g)
{
    DoubleMatrix psi0;
    DoubleMatrix psi1;

    psi0.resize(N2+1); for (unsigned int j=0; j<=N2; j++) psi0[j].resize(N1+1);
    psi1.resize(N2+1); for (unsigned int j=0; j<=N2; j++) psi1[j].resize(N1+1);

    double alpha1 = 0;
    double beta1  = 0;
    double alpha2 = 0;
    double beta2  = 0;

    DoubleVector a;
    DoubleVector b;
    DoubleVector c;
    DoubleVector d;
    DoubleVector x;

    ht = -ht;
//    a1 = -a1;
//    a2 = -a2;

    for (unsigned int k1=0; k1<=M; k1++)
    {
        unsigned int k = M-k1;

        if (k==M)
        {
            for (unsigned int j=0; j<=N2; j++)
            {
                for (unsigned i=0; i<=N1; i++)
                {
                    psi0[j][i] = -2.0*(uT[j][i] - U[j][i]);/*u(i*h1, j*h2, 1.0)*/;
                }
            }
        }
        else
        {
            // Approximation to x1 direction
            alpha1 = (a1*ht)/(2.0*h1*h1);
            beta1  = 1.0 - (a1*ht)/(h1*h1);
            alpha2 = -(a2*ht)/(2.0*h2*h2);
            beta2  = 1.0 + (a2*ht)/(h2*h2);

            a.resize(N1-1);
            b.resize(N1-1);
            c.resize(N1-1);
            d.resize(N1-1);
            x.resize(N1-1);

            for (unsigned int j=1; j<N2; j++)
            {
                for (unsigned int i=1; i<N1; i++)
                {
                    a[i-1] = alpha1;
                    b[i-1] = beta1;
                    c[i-1] = alpha1;
                    d[i-1] = alpha2*psi0[j-1][i] + beta2*psi0[j][i] + alpha2*psi0[j+1][i];
                }

                a[0]     = 0.0;
                c[N1-2]  = 0.0;
                d[0]    -= alpha1 * pm1(h2*j, ht*(k+0.5));
                d[N1-2] -= alpha1 * pm2(h2*j, ht*(k+0.5));

                TomasAlgorithm(a, b, c, d, x);

                psi1[j][0]  = pm1(h2*j, ht*(k+0.5));
                for (unsigned int i=1; i<N1; i++)
                {
                    psi1[j][i] = x[i-1];
                }
                psi1[j][N1] = pm2(h2*j, ht*(k+0.5));
            }

            for (unsigned int i=0; i<=N1; i++)
            {
                psi1[0][i]  = pm3(h1*i, ht*(k+0.5));
                psi1[N2][i] = pm4(h1*i, ht*(k+0.5));
            }

            a.clear();
            b.clear();
            c.clear();
            d.clear();
            x.clear();

            // Approximation to x2 direction
            alpha1 = +(a2*ht)/(2.0*h2*h2);
            beta1  = 1.0 - (a2*ht)/(h2*h2);
            alpha2 = -(a1*ht)/(2.0*h1*h1);
            beta2  = 1.0 + (a1*ht)/(h1*h1);

            a.resize(N2-1);
            b.resize(N2-1);
            c.resize(N2-1);
            d.resize(N2-1);
            x.resize(N2-1);

            for (unsigned int i=1; i<N1; i++)
            {
                for (unsigned int j=1; j<N2; j++)
                {
                    a[j-1] = alpha1;
                    b[j-1] = beta1;
                    c[j-1] = alpha1;
                    d[j-1] = alpha2*psi1[j][i-1] + beta2*psi1[j][i] + alpha2*psi1[j][i+1];
                }
                a[0]     = 0.0;
                c[N2-2]  = 0.0;
                d[0]    -= alpha1 * pm3(h1*i, ht*(k));
                d[N2-2] -= alpha1 * pm4(h1*i, ht*(k));
                TomasAlgorithm(a, b, c, d, x);

                psi0[0][i]  = pm3(h1*i, ht*(k));
                for (unsigned int j=1; j<N2; j++)
                {
                    psi0[j][i] = x[j-1];
                }
                psi0[N2][i] = pm4(h1*i, ht*(k));
            }

            for (unsigned int j=0; j<=N2; j++)
            {
                psi0[j][0]  = pm1(h2*j, ht*(k));
                psi0[j][N1] = pm2(h2*j, ht*(k));
            }

            a.clear();
            b.clear();
            c.clear();
            d.clear();
            x.clear();
        }

        for (unsigned int j=0; j<=N2; j++)
        {
            for (unsigned i=0; i<=N1; i++)
            {
                int index = k*(N1+1)*(N2+1)+j*(N1+1)+i;
                g[index] = -psi0[j][i] + 2*(f[index] - fxt(i*h1, j*h2, k*ht));
            }
        }
    }

//    a1 = -a1;
//    a2 = -a2;
    ht = -ht;

//    Printer::printMatrix(psi0, N2/10, N1/10);
//    printf("-----");
}

void HeatControl2D::calculateU1(const DoubleVector &f)
{
    DoubleMatrix u0;
    DoubleMatrix u1;

    u0.resize(N2+1); for (unsigned int j=0; j<=N2; j++) u0[j].resize(N1+1);
    u1.resize(N2+1); for (unsigned int j=0; j<=N2; j++) u1[j].resize(N1+1);

    double alpha1 = 0;
    double beta1  = 0;
    double alpha2 = 0;
    double beta2  = 0;

    DoubleVector a;
    DoubleVector b;
    DoubleVector c;
    DoubleVector d;
    DoubleVector x;

    ht = -ht;

    for (unsigned int k1=0; k1<=M; k1++)
    {
        unsigned int k = M-k1;

        if (k==M)
        {
            for (unsigned int j=0; j<=N2; j++)
            {
                for (unsigned i=0; i<=N1; i++)
                {
                    u0[j][i] = -2.0*(uT[j][i] - U[j][i]);/*u(i*h1, j*h2, 1.0)*/;
                }
            }
        }
        else
        {
            // Approximation to x1 direction
            alpha1 = (a1*ht)/(2.0*h1*h1);
            beta1  = 1.0 - (a1*ht)/(h1*h1);
            alpha2 = -(a2*ht)/(2.0*h2*h2);
            beta2  = 1.0 + (a2*ht)/(h2*h2);

            a.resize(N1-1);
            b.resize(N1-1);
            c.resize(N1-1);
            d.resize(N1-1);
            x.resize(N1-1);

            for (unsigned int j=1; j<N2; j++)
            {
                for (unsigned int i=1; i<N1; i++)
                {
                    a[i-1] = alpha1;
                    b[i-1] = beta1;
                    c[i-1] = alpha1;
                    d[i-1] = alpha2*u0[j-1][i] + beta2*u0[j][i] + alpha2*u0[j+1][i] - (ht/2.0) * f[k*(N1+1)*(N2+1)+j*(N1+1)+i]/*fxt(i*h1, j*h2, k*ht)*/;
                }

                a[0]     = 0.0;
                c[N1-2]  = 0.0;
                d[0]    -= alpha1 * m1(h2*j, ht*(k+0.5));
                d[N1-2] -= alpha1 * m2(h2*j, ht*(k+0.5));

                TomasAlgorithm(a, b, c, d, x);

                u1[j][0]  = m1(h2*j, ht*(k+0.5));
                for (unsigned int i=1; i<N1; i++)
                {
                    u1[j][i] = x[i-1];
                }
                u1[j][N1] = m2(h2*j, ht*(k+0.5));
            }

            for (unsigned int i=0; i<=N1; i++)
            {
                u1[0][i]  = m3(h1*i, ht*(k+0.5));
                u1[N2][i] = m4(h1*i, ht*(k+0.5));
            }

            a.clear();
            b.clear();
            c.clear();
            d.clear();
            x.clear();

            // Approximation to x2 direction
            alpha1 = +(a2*ht)/(2.0*h2*h2);
            beta1  = 1.0 - (a2*ht)/(h2*h2);
            alpha2 = -(a1*ht)/(2.0*h1*h1);
            beta2  = 1.0 + (a1*ht)/(h1*h1);

            a.resize(N2-1);
            b.resize(N2-1);
            c.resize(N2-1);
            d.resize(N2-1);
            x.resize(N2-1);

            for (unsigned int i=1; i<N1; i++)
            {
                for (unsigned int j=1; j<N2; j++)
                {
                    a[j-1] = alpha1;
                    b[j-1] = beta1;
                    c[j-1] = alpha1;
                    d[j-1] = alpha2*u1[j][i-1] + beta2*u1[j][i] + alpha2*u1[j][i+1] - (ht/2.0) * f[k*(N1+1)*(N2+1)+i*(N2+1)+j]/*fxt(i*h1, j*h2, k*ht)*/;
                }
                a[0]     = 0.0;
                c[N2-2]  = 0.0;
                d[0]    -= alpha1 * m3(h1*i, ht*(k));
                d[N2-2] -= alpha1 * m4(h1*i, ht*(k));
                TomasAlgorithm(a, b, c, d, x);

                u0[0][i]  = m3(h1*i, ht*(k));
                for (unsigned int j=1; j<N2; j++)
                {
                    u0[j][i] = x[j-1];
                }
                u0[N2][i] = m4(h1*i, ht*(k));
            }

            for (unsigned int j=0; j<=N2; j++)
            {
                u0[j][0]  = m1(h2*j, ht*(k));
                u0[j][N1] = m2(h2*j, ht*(k));
            }

            a.clear();
            b.clear();
            c.clear();
            d.clear();
            x.clear();
        }
    }

    ht = -ht;

    Printer::printMatrix(u0, N2/10, N1/10);
}

double HeatControl2D::u(double x1, double x2, double t)
{
    return x1*x1 + x2*x2 + t*t;
}

double HeatControl2D::fi(double x1, double x2)
{
    return u(x1, x2, t0);
}

double HeatControl2D::m1(double x2, double t)
{
    return u(x10, x2, t);
}

double HeatControl2D::m2(double x2, double t)
{
    return u(x11, x2, t);
}

double HeatControl2D::m3(double x1, double t)
{
    return u(x1, x20, t);
}

double HeatControl2D::m4(double x1, double t)
{
    return u(x1, x21, t);
}

double HeatControl2D::fxt(double x1, double x2, double t)
{
    return 2.0*t - 2.0*a1 - 2.0*a2;
}

void HeatControl2D::main()
{
    /* Function */
    HeatControl2D hc;

    DoubleVector f;
    f.resize(hc.C);
    for (unsigned int i=0; i<hc.C; i++) f[i] = (i/((hc.N1+1)*(hc.N2+1)))*hc.ht-4.0;

//    DoubleVector g;
//    puts("---");
//    hc.calculateU(f);
//    puts("---");
//    hc.calculateP(f, g);

    /* Minimization */
//    SteepestDescentGradient g1;
//    g1.setFunction(&hc);
//    g1.setEpsilon1(0.0000001);
//    g1.setEpsilon2(0.0000001);
//    g1.setGradientStep(0.000001);
//    g1.setR1MinimizeEpsilon(0.01, 0.0000001);
//    g1.setPrinter(new HeatControlPrinter);
//    g1.calculate(f0);

    /* Minimization */
    ConjugateGradient g2;
    g2.setFunction(&hc);
    g2.setEpsilon1(0.0000001);
    g2.setEpsilon2(0.0000001);
    g2.setGradientStep(0.000001);
    g2.setR1MinimizeEpsilon(0.01, 0.0000001);
    g2.setPrinter(new HeatControl2DPrinter);
    g2.setNormalize(false);
    g2.calculate(f);

    for (unsigned int c=0; c<f.size(); c++)
    {
        unsigned int k = c/((hc.N1+1)*(hc.N2+1));
        printf("%d %d %d\n", k, k*((hc.N1+1)*(hc.N2+1)), (k+1)*((hc.N1+1)*(hc.N2+1)));

//        for (unsigned int j1=k*((hc.N1+1)*(hc.N2+1)); j1<(k+1)*((hc.N1+1)*(hc.N2+1)); j1++)
//        {
//            DoubleMatrix m;
//            m.resize((hc.N2+1)*(hc.N1+1));

//            for (unsigned j=0; j<=hc.N2; j++)
//            {
//                for (unsigned i=0; i<=hc.N1; i++)
//                {
//                    m[j][i] = f[j1];
//                    Printer::printMatrix(m, hc.N2/10, hc.N1/10);
//                    printf("---");
//                }
//            }
//        }
    }

//    for (unsigned int j=0; j<=hc.M; j++)
//    {
//        if (j%(hc.M/10)==0)
//        {
//            printf("%6d|", j);
//            for (unsigned int i=0; i<=hc.N; i++)
//            {
//                if (i%(hc.N/10)==0)
//                    printf("%12.6f", f0[j*(hc.N+1)+i]);
//            }
//            puts("");
//        }
//    }
}

void HeatControl2DPrinter::print(unsigned int i, const DoubleVector &f0, const DoubleVector &s, double a, RnFunction *f) const
{
    printf("J: %.16f\n", f->fx(f0));
}


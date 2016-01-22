#include "parabolicequation.h"
#include "tomasmethod.h"

void IParabolicEquation::calculateU(DoubleVector &u, double hx, double ht, unsigned int N, unsigned int M, double a) const
{
    u.clear();
    u.resize(N+1);

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

    double alpha = -(a*ht)/(hx*hx);
    double beta  = 1.0 + (2.0*a*ht)/(hx*hx);

    for (unsigned int j=0; j<=M; j++)
    {
        if (j == 0)
        {
            for (unsigned int i=0; i<=N; i++)
            {
                u[i] = fi(i);
            }
        }
        else
        {
            for (unsigned int i=1; i<=N-1; i++)
            {
                a1[i-1] = alpha;
                b1[i-1] = beta;
                c1[i-1] = alpha;
                d1[i-1] = u[i] + ht * f(i, j);
            }

            a1[0]   = 0.0;
            c1[N-2] = 0.0;
            d1[0]   -= alpha * m1(j);
            d1[N-2] -= alpha * m2(j);

            TomasAlgorithm(a1, b1, c1, d1, x1);

            u[0] = m1(j);
            for (unsigned int i=1; i<=N-1; i++)
            {
                u[i] = x1[i-1];
            }
            u[N] = m2(j);
        }
    }

    a1.clear();
    b1.clear();
    c1.clear();
    d1.clear();
    x1.clear();
}

void IParabolicEquation::calculateU(DoubleMatrix &u, double hx, double ht, unsigned int N, unsigned int M, double a) const
{
    for (unsigned int i=0; i<u.size(); i++)
        u[i].clear();
    u.clear();

    u.resize(M+1);
    for (unsigned int i=0; i<u.size(); i++)
        u[i].resize(N+1);

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

    double alpha = -(a*ht)/(hx*hx);
    double beta  = 1.0 + (2.0*a*ht)/(hx*hx);

    for (unsigned int j=0; j<=M; j++)
    {
        if (j == 0)
        {
            for (unsigned int i=0; i<=N; i++)
            {
                u[j][i] = fi(i);
            }
        }
        else
        {
            for (unsigned int i=1; i<=N-1; i++)
            {
                a1[i-1] = alpha;
                b1[i-1] = beta;
                c1[i-1] = alpha;
                d1[i-1] = u[j-1][i] + ht * f(i, j);
            }

            a1[0]   = 0.0;
            c1[N-2] = 0.0;
            d1[0]   -= alpha * m1(j);
            d1[N-2] -= alpha * m2(j);

            TomasAlgorithm(a1, b1, c1, d1, x1);

            u[j][0] = m1(j);
            for (unsigned int i=1; i<=N-1; i++)
            {
                u[j][i] = x1[i-1];
            }
            u[j][N] = m2(j);
        }
    }

    a1.clear();
    b1.clear();
    c1.clear();
    d1.clear();
    x1.clear();
}

void IParabolicEquation2D::calculateU(DoubleMatrix &u, double h1, double h2, double ht, unsigned int N1, unsigned int N2, unsigned int M, double a1, double a2) const
{
    DoubleMatrix u0;
    DoubleMatrix u1;

    u0.resize(N2+1); for (unsigned int j=0; j<=N2; j++) u0[j].resize(N1+1);
    u1.resize(N2+1); for (unsigned int j=0; j<=N2; j++) u1[j].resize(N1+1);

    DoubleVector a;
    DoubleVector b;
    DoubleVector c;
    DoubleVector d;
    DoubleVector x;

    // Approximation to x1 direction
    double alpha11 = -(a2*ht)/(2.0*h2*h2);
    double beta11  = 1.0 + (a2*ht)/(h2*h2);
    double alpha12 = (a1*ht)/(2.0*h1*h1);
    double beta12  = 1.0 - (a1*ht)/(h1*h1);

    // Approximation to x2 direction
    double alpha21 = -(a1*ht)/(2.0*h1*h1);
    double beta21  = 1.0 + (a1*ht)/(h1*h1);
    double alpha22 = (a2*ht)/(2.0*h2*h2);
    double beta22  = 1.0 - (a2*ht)/(h2*h2);

    for (unsigned int k=0; k<=M; k++)
    {
        if (k==0)
        {
            for (unsigned int j=0; j<=N2; j++)
            {
                for (unsigned i=0; i<=N1; i++)
                {
                    u0[j][i] = fi(i, j);
                }
            }
        }
        else
        {
            // Approximation to x1 direction
            a.resize(N1-1);
            b.resize(N1-1);
            c.resize(N1-1);
            d.resize(N1-1);
            x.resize(N1-1);

            for (unsigned int i=1; i<N1; i++)
            {
                for (unsigned int j=1; j<N2; j++)
                {
                    a[j-1] = alpha11;
                    b[j-1] = beta11;
                    c[j-1] = alpha11;
                    d[j-1] = alpha12*u0[j][i-1] + beta12*u0[j][i] + alpha12*u0[j][i+1] + (ht/2.0) * f(i, j, k);
                }

                a[0]     = 0.0;
                c[N2-2]  = 0.0;
                d[0]    -= alpha11 * m3(i, k);
                d[N2-2] -= alpha11 * m4(i, k);

                TomasAlgorithm(a, b, c, d, x);

                u1[0][i]  = m3(i, k-0.5);
                for (unsigned int j=1; j<N2; j++)
                {
                    u1[j][i] = x[j-1];
                }
                u1[N2][i] = m4(i, k-0.5);
            }

            for (unsigned int j=0; j<=N2; j++)
            {
                u1[j][0]  = m1(j, k-0.5);
                u1[j][N1] = m2(j, k-0.5);
            }

            a.clear();
            b.clear();
            c.clear();
            d.clear();
            x.clear();

            // Approximation to x2 direction
            a.resize(N2-1);
            b.resize(N2-1);
            c.resize(N2-1);
            d.resize(N2-1);
            x.resize(N2-1);

            for (unsigned int j=1; j<N2; j++)
            {
                for (unsigned int i=1; i<N1; i++)
                {
                    a[i-1] = alpha21;
                    b[i-1] = beta21;
                    c[i-1] = alpha21;
                    d[i-1] = alpha22*u1[j-1][i] + beta22*u1[j][i] + alpha22*u1[j+1][i] + (ht/2.0) * f(i, j, k);
                }
                a[0]     = 0.0;
                c[N1-2]  = 0.0;
                d[0]    -= alpha21 * m1(j, k);
                d[N1-2] -= alpha21 * m2(j, k);

                TomasAlgorithm(a, b, c, d, x);

                u0[j][0]  = m1(j, k);
                for (unsigned int i=1; i<N1; i++)
                {
                    u0[j][i] = x[i-1];
                }
                u0[j][N1] = m2(j, k);
            }

            for (unsigned int i=0; i<=N1; i++)
            {
                u0[0][i]  = m3(i, k);
                u0[N2][i] = m4(i, k);
            }

            a.clear();
            b.clear();
            c.clear();
            d.clear();
            x.clear();
        }
    }

    u = u0;
}

//void ParabolicEquation2D::calculateBack(DoubleMatrix& u)
//{
//    u.Clear();

//    DoubleMatrix u0;
//    DoubleMatrix u1;

//    u0.resize(N2+1); for (unsigned int j=0; j<=N2; j++) u0[j].resize(N1+1);
//    u1.resize(N2+1); for (unsigned int j=0; j<=N2; j++) u1[j].resize(N1+1);

//    double alpha1 = 0;
//    double beta1  = 0;
//    double alpha2 = 0;
//    double beta2  = 0;

//    DoubleVector a;
//    DoubleVector b;
//    DoubleVector c;
//    DoubleVector d;
//    DoubleVector x;

//    for (unsigned int k1=0; k1<=M; k1++)
//    {
//        unsigned int k=M-k1;

//        if (k==M)
//        {
//            for (unsigned int j=0; j<=N2; j++)
//            {
//                for (unsigned i=0; i<=N1; i++)
//                {
//                    u0[j][i] = fi(i, j);
//                }
//            }
//        }
//        else
//        {
//            // Approximation to x1 direction
//            alpha1 = -(a2*ht)/(2.0*h2*h2);
//            beta1  = 1.0 + (a2*ht)/(h2*h2);
//            alpha2 = (a1*ht)/(2.0*h1*h1);
//            beta2  = 1.0 - (a1*ht)/(h1*h1);

//            a.resize(N1-1);
//            b.resize(N1-1);
//            c.resize(N1-1);
//            d.resize(N1-1);
//            x.resize(N1-1);

//            for (unsigned int i=1; i<N1; i++)
//            {
//                for (unsigned int j=1; j<N2; j++)
//                {
//                    a[j-1] = alpha1;
//                    b[j-1] = beta1;
//                    c[j-1] = alpha1;
//                    d[j-1] = alpha2*u0[j][i-1] + beta2*u0[j][i] + alpha2*u0[j][i+1] + (ht/2.0) * f(i, j, k);
//                }

//                a[0]     = 0.0;
//                c[N2-2]  = 0.0;
//                d[0]    -= alpha1 * m3(i, (k+0.5));
//                d[N2-2] -= alpha1 * m4(i, (k+0.5));

//                TomasAlgorithm(a, b, c, d, x);

//                u1[0][i]  = m3(i, k+0.5);
//                for (unsigned int j=1; j<N2; j++)
//                {
//                    u1[j][i] = x[j-1];
//                }
//                u1[N2][i] = m4(i, k+0.5);
//            }

//            for (unsigned int j=0; j<=N2; j++)
//            {
//                u1[j][0]  = m1(j, k+0.5);
//                u1[j][N1] = m2(j, k+0.5);
//            }

//            a.clear();
//            b.clear();
//            c.clear();
//            d.clear();
//            x.clear();

//            // Approximation to x2 direction
//            alpha1 = -(a1*ht)/(2.0*h1*h1);
//            beta1  = 1.0 + (a1*ht)/(h1*h1);
//            alpha2 = (a2*ht)/(2.0*h2*h2);
//            beta2  = 1.0 - (a2*ht)/(h2*h2);

//            a.resize(N2-1);
//            b.resize(N2-1);
//            c.resize(N2-1);
//            d.resize(N2-1);
//            x.resize(N2-1);

//            for (unsigned int j=1; j<N2; j++)
//            {
//                for (unsigned int i=1; i<N1; i++)
//                {
//                    a[i-1] = alpha1;
//                    b[i-1] = beta1;
//                    c[i-1] = alpha1;
//                    d[i-1] = alpha2*u1[j-1][i] + beta2*u1[j][i] + alpha2*u1[j+1][i] + (ht/2.0) * f(i, j, k);
//                }
//                a[0]     = 0.0;
//                c[N1-2]  = 0.0;
//                d[0]    -= alpha1 * m1(j, k);
//                d[N1-2] -= alpha1 * m2(j, k);
//                TomasAlgorithm(a, b, c, d, x);

//                u0[j][0]  = m1(j, k);
//                for (unsigned int i=1; i<N1; i++)
//                {
//                    u0[j][i] = x[i-1];
//                }
//                u0[j][N1] = m2(j, k);
//            }

//            for (unsigned int i=0; i<=N1; i++)
//            {
//                u0[0][i]  = m3(i, k);
//                u0[N2][i] = m4(i, k);
//            }

//            a.clear();
//            b.clear();
//            c.clear();
//            d.clear();
//            x.clear();
//        }
//    }

//    u = u0;
//}

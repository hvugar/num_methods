#include "borderparabolic2d.h"
#include <math.h>
#include <cmethods.h>

#define SAMPLE1

void BorderParabolic2D::Main(int argc UNUSED_PARAM, char **argv UNUSED_PARAM)
{
    BorderParabolic2D bp;
    bp.ht = 0.01;
    bp.h1 = 0.01;
    bp.h2 = 0.01;
    bp.M = 100;
    bp.N1 = 100;
    bp.N2 = 100;
    bp.a1 = bp.a2 = 1.0;

    DoubleMatrix m;
    bp.calculate(m, bp.h1, bp.h2, bp.ht, bp.N1, bp.N2, bp.M, bp.a1, bp.a2);
}

double BorderParabolic2D::u(unsigned int i UNUSED_PARAM, unsigned int j UNUSED_PARAM, unsigned int k UNUSED_PARAM) const
{
    double x1 UNUSED_PARAM = i*h1;
    double x2 UNUSED_PARAM = j*h2;
    double t  UNUSED_PARAM = 0.5*k*ht;
    return x1*x1 + x2*x2 + t*t;
}

double BorderParabolic2D::initial(unsigned int i UNUSED_PARAM, unsigned int j UNUSED_PARAM) const
{
#ifdef SAMPLE1
    return u(i, j, 0);
#endif
#ifdef SAMPLE2
    return 4.0;
#endif
#ifdef SAMPLE3
    return 0.0;
#endif
}

double BorderParabolic2D::boundary(unsigned int i UNUSED_PARAM, unsigned int j UNUSED_PARAM, unsigned int k UNUSED_PARAM) const
{
#ifdef SAMPLE1
    return u(i, j, k);
#endif
#ifdef SAMPLE2
    return initial(i,j)+3.0*t;
#endif
#ifdef SAMPLE3
    return 1.0;
#endif
}

double BorderParabolic2D::f(unsigned int i UNUSED_PARAM, unsigned int j UNUSED_PARAM, unsigned int k UNUSED_PARAM) const
{
    double x1 UNUSED_PARAM = i*h1;
    double x2 UNUSED_PARAM = j*h2;
    double t  UNUSED_PARAM = 0.5*k*ht;

    //static double sgm1 = 3.0*h1;
    //static double sgm2 = 3.0*h2;
    //static double gause_a = 1.0/(2.0*M_PI*sgm1*sgm2);
    //static double gause_b = 2.0*sgm1*sgm2;

#ifdef SAMPLE1
    return 2.0*t - 2.0*a1 - 2.0*a2;
#endif

#ifdef SAMPLE2
    double sum = 0.0;
    sum += (10*t) * gause_a * exp(-((x1-0.5)*(x1-0.5) + (x2-0.5)*(x2-0.5))/gause_b);
    //    if (i==N1/2 && j==N2/2) return (100*t)/(h1*h2);
#endif

#ifdef SAMPLE3
    //double sum = 0.0;
    //double power = 2000.0;
    //double delta1 = (1.0/h1);//(1.0/(sqrt(2.0*M_PI)*sgm1)) * exp(-((x1-k*ht)*(x1-k*ht))/(2.0*sgm1*sgm1));
    //double delta2 = (1.0/h2);//(1.0/(sqrt(2.0*M_PI)*sgm2)) * exp(-((x1-j*h2)*(x1-j*h2))/(2.0*sgm2*sgm2));
    //if ( (i==N1/2) && (j==(N2/2)) && (k <= 20) ) sum = power * delta1 * delta2;
    //if ( k <= 20) sum = 10000.0;
    return 0.0;
#endif
    return 0.0;
}

void BorderParabolic2D::calculate(DoubleMatrix &u, double h1, double h2, double ht, unsigned int N1, unsigned int N2, unsigned int M, double a1, double a2) const
{
    //cleaning matrix
    u.resize(N2+1, N1+1);

    DoubleMatrix uh(N2+1, N1+1);

    DoubleVector da1(N1-1);
    DoubleVector db1(N1-1);
    DoubleVector dc1(N1-1);
    DoubleVector dd1(N1-1);
    DoubleVector rx1(N1-1);

    DoubleVector da2(N2-1);
    DoubleVector db2(N2-1);
    DoubleVector dc2(N2-1);
    DoubleVector dd2(N2-1);
    DoubleVector rx2(N2-1);

    double x1_a = -(a1*a1*ht)/(2.0*h1*h1);
    double x1_b  = 1.0 + (a1*a1*ht)/(h1*h1);
    double x1_c = (a2*a2*ht)/(2.0*h2*h2);

    double x2_a = -(a2*a2*ht)/(2.0*h2*h2);
    double x2_b  = 1.0 + (a2*a2*ht)/(h2*h2);
    double x2_c = (a1*a1*ht)/(2.0*h1*h1);

    for (unsigned int k=0; k<=M; k++)
    {
        if (k==0)
        {
            for (unsigned int j=0; j<=N2; j++)
            {
                for (unsigned int i=0; i<=N1; i++)
                {
                    u[j][i] = initial(i, j);
                }
            }
        }
        else
        {
            // Approximation to x1 direction
            for (unsigned int j=1; j<N2; j++)
            {
                for (unsigned int i=1; i<N1; i++)
                {
                    da1[i-1] = x1_a;
                    db1[i-1] = x1_b;
                    dc1[i-1] = x1_a;
                    dd1[i-1] = x1_c*(u[j-1][i] - 2.0*u[j][i] + u[j+1][i]) + u[j][i] + (ht/2.0) * f(i, j, 2*k-1);
                }

                da1[0]     = 0.0;
                dc1[N1-2]  = 0.0;

                uh[j][0]  = boundary(0, j, 2*k-1);
                uh[j][N1] = boundary(N1, j, 2*k-1);

                dd1[0]    -= x1_a * uh[j][0];
                dd1[N1-2] -= x1_a * uh[j][N1];

                tomasAlgorithm(da1.data(), db1.data(), dc1.data(), dd1.data(), rx1.data(), rx1.length());

                for (unsigned int i=1; i<N1; i++)
                {
                    uh[j][i] = rx1[i-1];
                }
            }

            for (unsigned int i=0; i<=N1; i++)
            {
                uh[0][i]  = boundary(i, 0, 2*k-1);
                uh[N2][i] = boundary(i, N2, 2*k-1);
            }

            // Approximation to x2 direction
            for (unsigned int i=1; i<N1; i++)
            {
                for (unsigned int j=1; j<N2; j++)
                {
                    da2[j-1] = x2_a;
                    db2[j-1] = x2_b;
                    dc2[j-1] = x2_a;
                    dd2[j-1] = x2_c*(uh[j][i-1] - 2.0*uh[j][i] + uh[j][i+1]) + uh[j][i] + (ht/2.0) * f(i, j, 2*k);
                }
                da2[0]     = 0.0;
                dc2[N2-2]  = 0.0;

                u[0][i]  = boundary(i, 0, 2*k);
                u[N2][i] = boundary(i, N2, 2*k);

                dd2[0]    -= x2_a * u[0][i];
                dd2[N2-2] -= x2_a * u[N2][i];

                tomasAlgorithm(da2.data(), db2.data(), dc2.data(), dd2.data(), rx2.data(), rx2.length());

                for (unsigned int j=1; j<N2; j++)
                {
                    u[j][i] = rx2[j-1];
                }
            }

            for (unsigned int j=0; j<=N2; j++)
            {
                u[j][0]  = boundary(0, j, 2*k);
                u[j][N1] = boundary(N1, j, 2*k);
            }
        }

        IPrinter::printSeperatorLine();
        IPrinter::printMatrix(14,10,u);
    }

    da1.clear();
    db1.clear();
    dc1.clear();
    dd1.clear();
    rx1.clear();

    da2.clear();
    db2.clear();
    dc2.clear();
    dd2.clear();
    rx2.clear();
}

#include "borderhyperbolic2d.h"

#define SAMPLE1

void BorderHyperbolic2D::Main(int argc UNUSED_PARAM, char **argv UNUSED_PARAM)
{
    BorderHyperbolic2D bp;
    bp.h1 = 0.01;
    bp.h2 = 0.01;
    bp.ht = 0.00125;
    bp.N1 = 100;
    bp.N2 = 100;
    bp.M  = 800;
    bp.a1 = 1.0;
    bp.a2 = 1.0;

    DoubleMatrix m1;
    bp.calculate(m1, bp.h1, bp.h2, bp.ht, bp.N1, bp.N2, bp.M, bp.a1, bp.a2);
    IPrinter::printMatrix(14,10,m1);
    IPrinter::printSeperatorLine();

    DoubleMatrix m2;
    bp.calculateMVD(m2, bp.h1, bp.h2, bp.ht, bp.N1, bp.N2, bp.M, bp.a1, bp.a2);
    IPrinter::printMatrix(14,10,m2);
    IPrinter::printSeperatorLine();
}

double BorderHyperbolic2D::u(unsigned int i UNUSED_PARAM, unsigned int j UNUSED_PARAM, unsigned int k UNUSED_PARAM) const
{
    double x1 UNUSED_PARAM = i*h1;
    double x2 UNUSED_PARAM = j*h2;
    double t  UNUSED_PARAM = k*ht;
#ifdef SAMPLE1
    return x1*x1 + x2*x2 + t*t;
#endif
#ifdef SAMPLE2
    return x1*x1 + x1*x2 + sin(x1) + cos(x2) + exp(t);
#endif
#ifdef SAMPLE3
    return sin(x1)*sin(x1) + cos(x2) + x1*x2 + t*t*t;
#endif
#ifdef SAMPLE4
    return sin(t) + exp(x1) + x2*x2;
#endif
#ifdef SAMPLE5
    //return x1*x1 + x1*x2 + 5.0*sin(10.0*x1) + 3.0*cos(5.0*x2) + 4.0*cos(5.0*x1*t);
    return x1*x1 + x1*x2 + k1*sin(k2*x1) + k3*cos(k4*x2) + k5*cos(k6*x1*t);
#endif
#ifdef SAMPLE6
    return 2.0*sin(3.0*x1) + 4.0*cos(2.0*x2) + sin(2.0*t);
#endif
}

double BorderHyperbolic2D::initial1(unsigned int i UNUSED_PARAM, unsigned int j UNUSED_PARAM) const
{
    return u(i, j, 0);
}

double BorderHyperbolic2D::initial2(unsigned int i UNUSED_PARAM, unsigned int j UNUSED_PARAM) const
{
    C_UNUSED(i);
    C_UNUSED(j);
#ifdef SAMPLE1
    return 0.0;
#endif
#ifdef SAMPLE2
    return 1.0;
#endif
#ifdef SAMPLE3
    return 0.0;
#endif
#ifdef SAMPLE4
    return 1.0;
#endif
#ifdef SAMPLE5
    return 0.0;
#endif
#ifdef SAMPLE6
    return 2.0;
#endif
}

double BorderHyperbolic2D::boundary(unsigned int i, unsigned int j, unsigned int k) const
{
    return u(i, j, k);
}

double BorderHyperbolic2D::f(unsigned int i, unsigned int j, unsigned int k) const
{
    double x1 UNUSED_PARAM = i*h1;
    double x2 UNUSED_PARAM = j*h2;
    double t  UNUSED_PARAM = k*ht;

#ifdef SAMPLE1
    return 2.0 - 2.0*(a1*a1) - 2.0*(a2*a2);
#endif
#ifdef SAMPLE2
    return exp(t) + (a1*a1)*(sin(x1) - 2.0) + (a2*a2)*cos(x2);// + qamma*exp(t);
#endif
#ifdef SAMPLE3
    return 6*t - (a1*a1)*(2.0*cos(2.0*x1)) + (a2*a2)*cos(x2) ;
#endif
#ifdef SAMPLE4
    return -sin(t) - (a1*a1)*exp(x1) - (a2*a2)*2.0;
#endif
#ifdef SAMPLE5
    //return -100.0*x1*x1*cos(5.0*x1*t) + a1*a1*(500.0*sin(10.0*x1) + 100.0*t*t*cos(5.0*x1*t) - 2.0) + a2*a2*(75.0*cos(5.0*x2));
    return -k5*k6*k6*x1*x1*cos(k6*x1*t)
            - a1*a1*(2.0-k1*k2*k2*sin(k2*x1) - k5*k6*k6*t*t*cos(k6*x1*t))
            - a2*a2*(-k3*k4*k4*cos(k4*x2));
#endif
#ifdef SAMPLE6
    return -4.0*sin(2.0*t) + a1*a1*(18.0*sin(3.0*x1)) + a2*a2*(16.0*cos(2.0*x2)) + qamma*2.0*cos(2.0*t);
#endif
}

void BorderHyperbolic2D::calculate(DoubleMatrix &u, double h1, double h2, double ht, unsigned int N1, unsigned int N2, unsigned int M, double a1, double a2) const
{
    u.resize(N2+1, N1+1);

    DoubleMatrix u0(N2+1, N1+1);
    DoubleMatrix u1(N2+1, N1+1);

    DoubleVector da(N1-1);
    DoubleVector db(N1-1);
    DoubleVector dc(N1-1);
    DoubleVector dd(N1-1);
    DoubleVector rx(N1-1);

    double x1_a = -(a1*a1*ht*ht)/(h1*h1);
    double x1_b  = 1.0 + (2.0*a1*a1*ht*ht)/(h1*h1);
    double x1_c = (a2*a2*ht*ht)/(h2*h2);

    double x2_a = -(a2*a2*ht*ht)/(h2*h2);
    double x2_b  = 1.0 + (2.0*a2*a2*ht*ht)/(h2*h2);
    double x2_c = (a1*a1*ht*ht)/(h1*h1);

    for (unsigned int k=1; k<=M; k++)
    {
        if (k==1)
        {
            for (unsigned int j=0; j<=N2; j++)
            {
                for (unsigned int i=0; i<=N1; i++)
                {
                    u0[j][i] = initial1(i, j);
                    u1[j][i] = u0[j][i] + ht*initial2(i, j);
                }
            }

//            IPrinter::printMatrix(14,10,u0);
//            IPrinter::printSeperatorLine();
//            IPrinter::printMatrix(14,10,u1);
//            IPrinter::printSeperatorLine();
        }
        else
        {
            if ((k % 2) == 0)
            {
                for (unsigned int i=0; i<=N1; i++)
                {
                    u[0][i]  = boundary(i, 0, k);
                    u[N2][i] = boundary(i, N2, k);
                }

                // Approximation to x1 direction
                for (unsigned int j=1; j<N2; j++)
                {
                    for (unsigned int i=1; i<N1; i++)
                    {
                        da[i-1] = x1_a;
                        db[i-1] = x1_b;
                        dc[i-1] = x1_a;
                        dd[i-1] = x1_c*(u1[j-1][i] - 2.0*u1[j][i] + u1[j+1][i]) + 2.0*u1[j][i] - u0[j][i] + (ht*ht) * f(i, j, k);
                    }

                    da[0]     = 0.0;
                    dc[N1-2]  = 0.0;

                    u[j][0]  = boundary(0, j, k);
                    u[j][N1] = boundary(N1, j, k);

                    dd[0]    -= x1_a * u[j][0];
                    dd[N1-2] -= x1_a * u[j][N1];

                    tomasAlgorithm(da.data(), db.data(), dc.data(), dd.data(), rx.data(), rx.size());

                    for (unsigned int i=1; i<N1; i++)
                    {
                        u[j][i] = rx[i-1];
                    }
                }

            }
            else
            {
                for (unsigned int j=0; j<=N2; j++)
                {
                    u[j][0]  = boundary(0, j, k);
                    u[j][N1] = boundary(N1, j, k);
                }

                // Approximation to x2 direction
                for (unsigned int i=1; i<N1; i++)
                {
                    for (unsigned int j=1; j<N2; j++)
                    {
                        da[j-1] = x2_a;
                        db[j-1] = x2_b;
                        dc[j-1] = x2_a;
                        dd[j-1] = x2_c*(u1[j][i-1] - 2.0*u1[j][i] + u1[j][i+1]) + 2.0*u1[j][i] - u0[j][i] + (ht*ht) * f(i, j, k);
                    }
                    da[0]     = 0.0;
                    dc[N2-2]  = 0.0;

                    u[0][i]  = boundary(i, 0, k);
                    u[N2][i] = boundary(i, N2, k);

                    dd[0]    -= x2_a * u[0][i];
                    dd[N2-2] -= x2_a * u[N2][i];

                    tomasAlgorithm(da.data(), db.data(), dc.data(), dd.data(), rx.data(), rx.size());

                    for (unsigned int j=1; j<N2; j++)
                    {
                        u[j][i] = rx[j-1];
                    }
                }
            }

            for (unsigned int j=0; j<=N2; j++)
            {
                for (unsigned int i=0; i<=N1; i++)
                {
                    u0[j][i] = u1[j][i];
                    u1[j][i] = u[j][i];
                }
            }

//            IPrinter::printMatrix(14,10,u1);
//            IPrinter::printSeperatorLine();
        }
    }

    da.clear();
    db.clear();
    dc.clear();
    dd.clear();
    rx.clear();

//    da2.clear();
//    db2.clear();
//    dc2.clear();
//    dd2.clear();
//    rx2.clear();
}

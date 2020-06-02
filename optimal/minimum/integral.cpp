#include "integral.h"

double Trapesium(R1Function *f, size_t N, double a, double b)
{
    double sum = 0.0;
    double h = (a-b)/N;
    for (size_t i=0; i<=N-1; i++)
    {
        double f1 = f->fx(a+i*h);
        double f2 = f->fx(a+(i+1)*h);
        sum = sum + (f1+f2);
    }
    sum = (h/2.0)*sum;
    return sum;
}

double Trapesium(R2Function *f, size_t N1, size_t N2, double a1, double b1, double a2, double b2)
{
    double sum = 0.0;
    double h1 = (a1-b1)/N1;
    double h2 = (a2-b2)/N2;

    for (size_t j=0; j<=N2-1; j++)
    {
        for (size_t i=0; i<=N1-1; i++)
        {
            double f11 = f->fx(a1+(i+0)*h1, a2+(j+0)*h2);
            double f12 = f->fx(a1+(i+0)*h1, a2+(j+1)*h2);
            double f21 = f->fx(a1+(i+1)*h1, a2+(j+0)*h2);
            double f22 = f->fx(a1+(i+1)*h1, a2+(j+1)*h2);
            sum = sum + (f11+f12+f21+f22);
        }
    }
    sum = (h1*h2*0.25)*sum;
    return sum;
}

double Trapesium(R3Function *f, size_t N1, size_t N2, size_t N3, double a1, double b1, double a2, double b2, double a3, double b3)
{
    double sum = 0.0;
    double h1 = (a1-b1)/N1;
    double h2 = (a2-b2)/N2;
    double h3 = (a3-b3)/N3;

    for (size_t k=0; k<=N3-1; k++)
    {
        for (size_t j=0; j<=N2-1; j++)
        {
            for (size_t i=0; i<=N1-1; i++)
            {
                double f111 = f->fx(a1+(i+0)*h1, a2+(j+0)*h2, a3+(k+0)*h3);
                double f112 = f->fx(a1+(i+0)*h1, a2+(j+0)*h2, a3+(k+1)*h3);
                double f121 = f->fx(a1+(i+0)*h1, a2+(j+1)*h2, a3+(k+0)*h3);
                double f122 = f->fx(a1+(i+0)*h1, a2+(j+1)*h2, a3+(k+1)*h3);
                double f211 = f->fx(a1+(i+1)*h1, a2+(j+0)*h2, a3+(k+0)*h3);
                double f212 = f->fx(a1+(i+1)*h1, a2+(j+0)*h2, a3+(k+1)*h3);
                double f221 = f->fx(a1+(i+1)*h1, a2+(j+1)*h2, a3+(k+0)*h3);
                double f222 = f->fx(a1+(i+1)*h1, a2+(j+1)*h2, a3+(k+1)*h3);
                sum = sum + (f111+f112+f121+f122+f211+f212+f221+f222);
            }
        }
    }
    sum = (h1*h2*h3*0.125)*sum;
    return sum;
}

double trapesium(R1Function *f, size_t N, double a, double b)
{
    double h = (b - a) / N;
    double sum = 0;
    sum = sum + f->fx(a) + f->fx(b);
    for (size_t i=1; i<=N-1; i++)
    {
        sum += f->fx(a + i*h);
    }
    return sum;
}

double trapesium2D(const DoubleMatrix &m, double hx, double hy)
{
    //    size_t rows = m.rows();
    //    size_t cols = m.cols();

    double sum = 0.0;
    //    sum += 0.25 * m[0][0];
    //    udiff = u1[0][N]; usum += 0.25 * udiff * udiff;// * mu(N, 0);
    //    udiff = u1[M][0]; usum += 0.25 * udiff * udiff;// * mu(0, M);
    //    udiff = u1[M][N]; usum += 0.25 * udiff * udiff;// * mu(N, M);

    //    for (size_t n=1; n<=N-1; n++)
    //    {
    //        udiff = u1[0][n]; usum += 0.5 * udiff * udiff;// * mu(n, 0);
    //        udiff = u1[M][n]; usum += 0.5 * udiff * udiff;// * mu(n, M);
    //    }

    //    for (size_t m=1; m<=M-1; m++)
    //    {
    //        udiff = u1[m][0]; usum += 0.5 * udiff * udiff;// * mu(0, m);
    //        udiff = u1[m][N]; usum += 0.5 * udiff * udiff;// * mu(N, m);
    //    }

    //    for (size_t m=1; m<=M-1; m++)
    //    {
    //        for (size_t n=1; n<=N-1; n++)
    //        {
    //            udiff = u1[m][n]; usum += udiff * udiff;// * mu(n, m);
    //        }
    //    }


    return sum;
}

double Integral::rectangle(const DoubleVector &v)
{
    return 0.0;
}

double Integral::trapezoidal(const DoubleVector &v)
{
    return 0.0;
}

double Integral::simpsons(const DoubleVector &v)
{
    return 0.0;
}

double Integral::rectangle(const DoubleMatrix &m)
{
    return 0.0;
}

double Integral::trapezoidal(const DoubleMatrix &m, double hx, double hy)
{
    size_t M = m.rows() - 1;
    size_t N = m.cols() - 1;

    double sum = 0.0;

    sum += 0.25*m[0][0];
    sum += 0.25*m[M][0];
    sum += 0.25*m[M][N];
    sum += 0.25*m[0][N];

    for (size_t i=0; i<=N; i++)
    {
        sum += 0.5*m[0][i];
        sum += 0.5*m[M][i];
    }

    for (size_t j=0; j<=M; j++)
    {
        sum += 0.5*m[j][0];
        sum += 0.5*m[j][N];
    }

    for (size_t j=1; j<M; j++)
    {
        for (size_t i=0; i<N; i++)
        {
            sum += m[j][i];
        }
    }

    return sum*hx*hy;
}

double Integral::simpsons(const DoubleMatrix &m, double hx, double hy)
{
    size_t M = m.rows() - 1;
    size_t N = m.cols() - 1;

    double sum = 0.0;
    for (size_t j=0; j<=M; j++)
    {
        for (size_t i=0; i<=N; i++)
        {
            double fx = m[j][i];
            if (i==0 || i==N) { fx *= 1.0; } else { if (i%2==1) { fx *= 4.0; } else { fx *= 2.0; } }
            if (j==0 || j==M) { fx *= 1.0; } else { if (j%2==1) { fx *= 4.0; } else { fx *= 2.0; } }
            sum += fx;
        }
    }
    sum *= ((hx*hy)/9.0);

    return sum;
}

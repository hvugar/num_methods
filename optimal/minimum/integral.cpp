#include "integral.h"

double Trapesium(R1Function *f, unsigned int N, double a, double b)
{
    double sum = 0.0;
    double h = (a-b)/N;
    for (unsigned int i=0; i<=N-1; i++)
    {
        double f1 = f->fx(a+i*h);
        double f2 = f->fx(a+(i+1)*h);
        sum = sum + (f1+f2);
    }
    sum = (h/2.0)*sum;
    return sum;
}

double Trapesium(R2Function *f, unsigned int N1, unsigned int N2, double a1, double b1, double a2, double b2)
{
    double sum = 0.0;
    double h1 = (a1-b1)/N1;
    double h2 = (a2-b2)/N2;

    for (unsigned int j=0; j<=N2-1; j++)
    {
        for (unsigned int i=0; i<=N1-1; i++)
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

double Trapesium(R3Function *f, unsigned int N1, unsigned int N2, unsigned int N3, double a1, double b1, double a2, double b2, double a3, double b3)
{
    double sum = 0.0;
    double h1 = (a1-b1)/N1;
    double h2 = (a2-b2)/N2;
    double h3 = (a3-b3)/N3;

    for (unsigned int k=0; k<=N3-1; k++)
    {
        for (unsigned int j=0; j<=N2-1; j++)
        {
            for (unsigned int i=0; i<=N1-1; i++)
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

double trapesium(R1Function *f, unsigned int N, double a, double b)
{
    double h = (b - a) / N;
    double sum = 0;
    sum = sum + f->fx(a) + f->fx(b);
    for (unsigned int i=1; i<=N-1; i++)
    {
        sum += f->fx(a + i*h);
    }
    return sum;
}

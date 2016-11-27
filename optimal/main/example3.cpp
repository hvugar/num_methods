#include "example3.h"

#define SAMPLE_1
#define SAMPLE_2
#define SAMPLE_3

Example3::Example3()
{
    calculate1();
}

double Example3::a(unsigned int i, unsigned int j, unsigned int k UNUSED_PARAM) const
{
    if (i==1)
    {
        if (j==1) return 2.0;
        if (j==2) return 3.0;
    }
    if (i==2)
    {
        if (j==1) return 3.0;
        if (j==2) return 4.0;
    }
    return 0.0;
}

double Example3::b(unsigned int i, unsigned int k) const
{
    double t = k*h;
    if (i==1)
    {
#ifdef SAMPLE_1
        //return 2.0*t - 2.0*t*t - 6.0*sin(t);
#endif
#ifdef SAMPLE_2
        //return -3.0*t*t*t - 2.0*t*t - 4.0*t;
        return 2.0*cos(t)-3.0*cos(2.0*t)-4.0*sin(t)-4.0*t-2.0*t*t;
    }
#endif
    if (i==2)
    {
        //return 2.0*cos(t) - 3.0*t*t - 8.0*sin(t);
        //return 2.0 - 4.0*t*t*t - 8.0*t;
        return 2.0-2.0*sin(2.0*t)-6.0*sin(t)-4.0*cos(2.0*t)-3.0*t*t-8.0*t;
    }
    return 0.0;
}

double Example3::X1(unsigned int k) const
{
    double t = k*h;
    //return t*t;
    //return t*t;
    return 2.0*sin(t)+t*t;
}

double Example3::X2(unsigned int k) const
{
    double t = k*h;
    //return 2.0*sin(t);
    //return t*t*t+2.0*t;
    return cos(2.0*t)+2.0*t;
}

void Example3::calculate()
{
    h = 0.01;
    N = 100;

    DoubleVector x01(N+1);
    DoubleVector x02(N+2);
    for (unsigned int k=0; k<=N; k++)
    {
        x01.at(k) = X1(k);
        x02.at(k) = X2(k);
    }
    IPrinter::printVector(14,10,x01);
    IPrinter::printVector(14,10,x02);

    DoubleVector x1(N+1);
    DoubleVector x2(N+2);

    x1.at(N-0) = X1(N-0); x2.at(N-0) = X2(N-0);
    x1.at(N-1) = X1(N-1); x2.at(N-1) = X2(N-1);
    x1.at(N-2) = X1(N-2); x2.at(N-2) = X2(N-2);
    x1.at(N-3) = X1(N-3); x2.at(N-3) = X2(N-3);

    x1.at(0) = X1(0); x2.at(0) = X2(0);
    x1.at(1) = X1(1); x2.at(1) = X2(1);
    x1.at(2) = X1(2); x2.at(2) = X2(2);
    x1.at(3) = X1(3); x2.at(3) = X2(3);

    //    for (unsigned int k=0; k<=N; k++)
    //    {
    //        x1.at(k) = X1(k);
    //        x2.at(k) = X2(k);
    //    }
    //    IPrinter::printVector(14,10,x1);
    //    IPrinter::printVector(14,10,x2);

    //    for (unsigned int k=0; k<=N-4; k++)
    //    {
    //        x1.at(k+4) = (16.0/3.0)*x1.at(k+3) + (-12.0)*x1.at(k+2) + (16.0)*x1.at(k+1)  + (-25.0/3.0-4.0*h*a(1,1,k))*x1.at(k) + (-4.0*h*a(1,2,k))*x2.at(k) + (-4.0*h*b(1,k));
    //        x2.at(k+4) = (16.0/3.0)*x2.at(k+3) + (-12.0)*x2.at(k+2) + (16.0)*x2.at(k+1)  + (-25.0/3.0-4.0*h*a(2,2,k))*x2.at(k) + (-4.0*h*a(2,1,k))*x1.at(k) + (-4.0*h*b(2,k));
    //    }

    //    for (unsigned int k=N; k>=4; k--)
    //    {
    //        x1.at(k-4) = (16.0/3.0)*x1.at(k-3) + (-12.0)*x1.at(k-2) + (16.0)*x1.at(k-1)  + (-25.0/3.0+4.0*h*a(1,1,k))*x1.at(k) + (+4.0*h*a(1,2,k))*x2.at(k) + (+4.0*h*b(1,k));
    //        x2.at(k-4) = (16.0/3.0)*x2.at(k-3) + (-12.0)*x2.at(k-2) + (16.0)*x2.at(k-1)  + (-25.0/3.0+4.0*h*a(2,2,k))*x2.at(k) + (+4.0*h*a(2,1,k))*x1.at(k) + (+4.0*h*b(2,k));
    //    }

    puts("---");
    for (unsigned int k=N-4; k!=UINT32_MAX; k--)
    {
        unsigned int k1 = k+1;
        x1.at(k) = (48.0/25.0)*x1.at(k+1) + (-36.0/25.0)*x1.at(k+2) + (16.0/25.0)*x1.at(k+3) + (-3.0/25.0)*x1.at(k+4)
                + ((-12.0/25.0)*h*a(1,1,k1))*x1.at(k1) + ((-12.0/25.0)*h*a(1,2,k1))*x2.at(k1) + ((-12.0/25.0)*h*b(1,k1));
        x2.at(k) = (48.0/25.0)*x2.at(k+1) + (-36.0/25.0)*x2.at(k+2) + (16.0/25.0)*x2.at(k+3) + (-3.0/25.0)*x2.at(k+4)
                + ((-12.0/25.0)*h*a(2,1,k1))*x1.at(k1) + ((-12.0/25.0)*h*a(2,2,k1))*x2.at(k1) + ((-12.0/25.0)*h*b(2,k1));
    }
    IPrinter::printVector(14,10,x1);
    IPrinter::printVector(14,10,x2);
    puts("---");
    for (unsigned int k=N-4; k!=UINT32_MAX; k--)
    {
        unsigned int k1 = k+2;
        x1.at(k) = (48.0/25.0)*x1.at(k+1) + (-36.0/25.0)*x1.at(k+2) + (16.0/25.0)*x1.at(k+3) + (-3.0/25.0)*x1.at(k+4)
                + ((-12.0/25.0)*h*a(1,1,k1))*x1.at(k1) + ((-12.0/25.0)*h*a(1,2,k1))*x2.at(k1) + ((-12.0/25.0)*h*b(1,k1));
        x2.at(k) = (48.0/25.0)*x2.at(k+1) + (-36.0/25.0)*x2.at(k+2) + (16.0/25.0)*x2.at(k+3) + (-3.0/25.0)*x2.at(k+4)
                + ((-12.0/25.0)*h*a(2,1,k1))*x1.at(k1) + ((-12.0/25.0)*h*a(2,2,k1))*x2.at(k1) + ((-12.0/25.0)*h*b(2,k1));
    }
    IPrinter::printVector(14,10,x1);
    IPrinter::printVector(14,10,x2);
    puts("---");
    for (unsigned int k=N-4; k!=UINT32_MAX; k--)
    {
        unsigned int k1 = k+3;
        x1.at(k) = (48.0/25.0)*x1.at(k+1) + (-36.0/25.0)*x1.at(k+2) + (16.0/25.0)*x1.at(k+3) + (-3.0/25.0)*x1.at(k+4)
                + ((-12.0/25.0)*h*a(1,1,k1))*x1.at(k1) + ((-12.0/25.0)*h*a(1,2,k1))*x2.at(k1) + ((-12.0/25.0)*h*b(1,k1));
        x2.at(k) = (48.0/25.0)*x2.at(k+1) + (-36.0/25.0)*x2.at(k+2) + (16.0/25.0)*x2.at(k+3) + (-3.0/25.0)*x2.at(k+4)
                + ((-12.0/25.0)*h*a(2,1,k1))*x1.at(k1) + ((-12.0/25.0)*h*a(2,2,k1))*x2.at(k1) + ((-12.0/25.0)*h*b(2,k1));
    }
    IPrinter::printVector(14,10,x1);
    IPrinter::printVector(14,10,x2);
    puts("---");
}

void Example3::calculate1()
{
    h = 0.0001;
    N = 10000;

    DoubleVector x01(N+1);
    DoubleVector x02(N+2);
    for (unsigned int k=0; k<=N; k++)
    {
        x01.at(k) = X1(k);
        x02.at(k) = X2(k);
    }
    IPrinter::printVector(14,10,x01);
    IPrinter::printVector(14,10,x02);
    puts("---");

    DoubleVector x1(N+1);
    DoubleVector x2(N+2);
    x1.at(0) = X1(0); x2.at(0) = X2(0);
    for (unsigned int k=0; k<=N-1; k++)
    {
        x1.at(k+1) = x1.at(k) + (h*a(1,1,k))*x1.at(k) + (h*a(1,2,k))*x2.at(k) + (h*b(1,k));
        x2.at(k+1) = x2.at(k) + (h*a(2,1,k))*x1.at(k) + (h*a(2,2,k))*x2.at(k) + (h*b(2,k));
    }
    IPrinter::printVector(14,10,x1);
    IPrinter::printVector(14,10,x2);
    puts("---");

    //DoubleVector x1(N+1);
    //DoubleVector x2(N+2);
    x1.at(N-0) = X1(N-0); x2.at(N-0) = X2(N-0);
    for (unsigned int k=N; k>=1; k--)
    {
        x1.at(k-1) = x1.at(k) - ((h*a(1,1,k))*x1.at(k) + (h*a(1,2,k))*x2.at(k) + (h*b(1,k)));
        x2.at(k-1) = x2.at(k) - ((h*a(2,1,k))*x1.at(k) + (h*a(2,2,k))*x2.at(k) + (h*b(2,k)));
    }
    IPrinter::printVector(14,10,x1);
    IPrinter::printVector(14,10,x2);
    puts("---");
}

#include "example3.h"

void Example3::Main(int argc, char *argv[])
{
    C_UNUSED(argc);
    C_UNUSED(argv);

    Example3 e3;
}

Example3::Example3()
{
    DoubleVector x;
    x << 2.50 << 2.60; //k
    x << 2.89 << 2.81; //z
    x << 0.40 << 0.70; //x

    sx = &x;

    DoubleMatrix u;
    calculateU(u, x);
    V = u.row(M);

    fx(x);
}

double Example3::fx(const DoubleVector &x)
{
    px = &x;

    DoubleMatrix u;
    calculateU(u, x);

    double sum = 0.0;
    for (unsigned int n=0; n<N; n++)
    {
        unsigned int n1 = n+0;
        unsigned int n2 = n+1;

        double f1 = mu(n1)*(u.at(M, n1)-V.at(n1))*(u.at(M, n1)-V.at(n1));
        double f2 = mu(n2)*(u.at(M, n2)-V.at(n2))*(u.at(M, n2)-V.at(n2));
        sum = sum + (f1 + f2);
    }
    sum = 0.5*dx*sum;

    DoubleVector k = x.mid(0, 1);
    DoubleVector z = x.mid(2, 3);
    DoubleVector e = x.mid(4, 5);

    return alpha0*sum + alpha1*k.EuclideanNorm() + alpha2*z.EuclideanNorm() + alpha3*e.EuclideanNorm();
}

void Example3::gradient(const DoubleVector &x, DoubleVector &g)
{}

void Example3::print(unsigned int iteration, const DoubleVector &x, const DoubleVector &gradient, double alpha, RnFunction *fn) const
{}

void Example3::project(DoubleVector &x, int index)
{}

void Example3::calculateU(DoubleMatrix &u, const DoubleVector &x)
{
    u.clear();
    u.resize(M+1, N+1);

    for (unsigned int m=0; m<=M; m++)
    {
        if (m==0)
        {
            for (unsigned int n=0; n<=N; n++)
            u.at(0,n) = initial(n);
        }
        else
        {
            DoubleVector da(N+1);
            DoubleVector db(N+1);
            DoubleVector dc(N+1);
            DoubleVector dd(N+1);
        }
    }
}

double Example3::initial(unsigned int i UNUSED_PARAM) const
{
    return Ti;
}

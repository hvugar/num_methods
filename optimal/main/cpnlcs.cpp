#include "cpnlcs.h"

CauchyProblemNonLocalContions::CauchyProblemNonLocalContions()
{
    times << 0.0 << 0.5 << 1.0;

    alpa0.resize(n,n);
    alpa0[0][0] = 1.0; alpa0[0][1] = 3.0;
    alpa0[1][0] = 2.0; alpa0[1][1] = 1.0;

    alpa1.resize(n,n);
    alpa1[0][0] = 2.0; alpa1[0][1] = 1.0;
    alpa1[1][0] = 0.0; alpa1[1][1] = 0.0;

    alpa2.resize(n,n);
    alpa2[0][0] = 5.0; alpa2[0][1] = 4.0;
    alpa2[1][0] = 8.0; alpa0[1][1] = 1.0;

    betta.resize(n);
    betta[0] = alpa0[0][0]*x1(0)   + alpa0[0][1]*x2(0) +
               alpa1[0][0]*x1(N/2) + alpa1[0][1]*x2(N/2) +
               alpa2[0][0]*x1(N)   + alpa2[0][1]*x2(N);
    betta[1] = alpa0[1][0]*x1(0)   + alpa0[1][1]*x2(0) +
               alpa1[1][0]*x1(N/2) + alpa1[1][1]*x2(N/2) +
               alpa2[1][0]*x1(N)   + alpa2[1][1]*x2(N);
}

double CauchyProblemNonLocalContions::A(unsigned int k, unsigned int i, unsigned int j) const
{
    double t = k*h;
    if (i==0)
    {
        if (j==0) { return t; }
        if (j==1) { return -1.0; }
    }
    if (i==1)
    {
        if (j==0) { return +2.0*t+3.0; }
        if (j==1) { return -2.0; }
    }
    return 0.0;
}

double CauchyProblemNonLocalContions::B(unsigned int k, unsigned int i) const
{
    double t = k*h;
    if (i==0)
    {
        return -t;
    }
    if (i==1)
    {
        return -6.0*t-11.0;
    }
    return 0.0;
}

double CauchyProblemNonLocalContions::x1(unsigned int k) const
{
    double t = k*h;
    return t*t + 4.0;
}

double CauchyProblemNonLocalContions::x2(unsigned int k) const
{
    double t = k*h;
    return t*t*t + t;
}

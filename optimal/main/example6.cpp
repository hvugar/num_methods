#include "example6.h"

Example6::Example6()
{
    s[0].push_back(0); s[0].push_back(2*F);   s[0].push_back(5*F);   s[0].push_back(8*F);   s[0].push_back(9*F); s[0].push_back(10*e.F);
    s[1].push_back(0); s[1].push_back(1);     s[1].push_back(2);     s[1].push_back(3);     s[1].push_back(4);
    s[2].push_back(0); s[2].push_back(1);     s[2].push_back(2);     s[2].push_back(3);     s[2].push_back(4);
    s[3].push_back(0); s[3].push_back(1);     s[3].push_back(2);     s[3].push_back(3);     s[3].push_back(4);
}

double Example6::a(unsigned int i, unsigned int j, unsigned int k) const
{
    if (i==1 && j==1) return +2.0;
    if (i==1 && j==2) return +3.0;
    if (i==1 && j==3) return -1.0;

    if (i==2 && j==1) return +4.0;
    if (i==2 && j==2) return +6.0;
    if (i==2 && j==3) return -2.0;

    if (i==3 && j==1) return -1.0;
    if (i==3 && j==2) return +1.0;
    if (i==3 && j==3) return -1.0;

    return NAN;
}

double Example6::b(unsigned int i, unsigned int k) const
{
    double t = k*h;
    if (i==1) return -(+2.0*sin(20.0*t*t) + 3.0*(cos(10.0*t) - sin(20.0*t)) - 1.0*(t*t*t - sin(8.0*t)*sin(8.0*t))) + (+40.0*t*cos(20.0*t*t));
    if (i==2) return -(+4.0*sin(20.0*t*t) + 6.0*(cos(10.0*t) - sin(20.0*t)) - 2.0*(t*t*t - sin(8.0*t)*sin(8.0*t))) + (-10.0*sin(10.0*t) - 20.0*cos(20.0*t));
    if (i==3) return -(-1.0*sin(20.0*t*t) + 1.0*(cos(10.0*t) - sin(20.0*t)) - 1.0*(t*t*t - sin(8.0*t)*sin(8.0*t))) + (+3.0*t*t - 16.0*cos(8.0*t)*sin(8.0*t));
    return NAN;
}

double Example6::fx(unsigned int i, unsigned int k) const
{
    if (i==1) return sin(20.0*t*t);
    if (i==2) return cos(10.0*t) - sin(20.0*t);
    if (i==3) return t*t*t - sin(8.0*t)*sin(8.0*t);
    return NAN;
}

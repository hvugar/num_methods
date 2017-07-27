#include "problem4.h"

Problem4::Problem4(const Dimension &time) : mtime(time)
{
}

void Problem4::calculate(DoubleVector &x)
{

}

double Problem4::X(unsigned int k) const
{
    double x = mtime.step()*k;
    return x*x;
}

double Problem4::a(const TimeNode &tn) const
{
    return 3.0*tn.t;
}

double Problem4::b(const TimeNode &tn) const
{
    double t = tn.t;
    return 2.0*t - 3.0*t*t*t*t - 0.00390625*t - 0.177978515625*t*t;
}

double Problem4::b1(const TimeNode &tn) const
{
    return tn.t;
}

double Problem4::b2(const TimeNode &tn) const
{
    return tn.t*tn.t;
}


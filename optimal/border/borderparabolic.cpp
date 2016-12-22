#include "borderparabolic.h"

void BorderParabolic::Main(int argc UNUSED_PARAM, char *argv[] UNUSED_PARAM)
{
    DoubleMatrix u;
    BorderParabolic bp;
    bp.calculateU(u, bp.hx, bp.ht, bp.N, bp.M, bp.a);
    IPrinter::printMatrix(u, 10, 10, NULL);
}

double BorderParabolic::initial(unsigned int i UNUSED_PARAM) const
{
    double x = i*hx;
    return x*x;
}

double BorderParabolic::boundary(Boundary type UNUSED_PARAM, unsigned int j UNUSED_PARAM) const
{
    double t = j*ht;
    if (type == Left) return t*t;
    if (type == Right) return t*t + 1.0;
    return 0.0;
}

double BorderParabolic::f(unsigned int i UNUSED_PARAM, unsigned int j UNUSED_PARAM) const
{
    double t = j*ht;
    return 2.0*t - 2.0*a*a;
}

void BorderParabolic::calculate()
{

}

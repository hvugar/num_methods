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

void BorderParabolic::calculate(DoubleMatrix &u)
{
    u.resize(M+1, N+1);
    for (unsigned int i=0; i<=N; i++) u.at(0,i) = initial(i);
    for (unsigned int j=1; j<=M; j++) { u.at(j,0) = boundary(Left,j); u.at(j,N) = boundary(Right, j); }


    double alpha = (ht*a*a)/(24.0*hx*hx);
}

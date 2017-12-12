#include "cproblem2backward2d.h"

double CProblem2Backward2D::initial(const SpaceNodePDE &sn) const
{
    unsigned int i = sn.i;
    unsigned int j = sn.j;
    return -2.0 * mu[j][i] * (uT[j][i]-U[j][i]) + h(sn);
}

double CProblem2Backward2D::boundary(const SpaceNodePDE &sn UNUSED_PARAM, const TimeNodePDE &tn UNUSED_PARAM, BoundaryType boundary UNUSED_PARAM) const
{
    return NAN;
}

double CProblem2Backward2D::f(const SpaceNodePDE &sn, const TimeNodePDE &tn) const
{
    double x = sn.x;
    double y = sn.y;
    double t = tn.t;

    double res = 1.0 + 4.0*a*a - lambda0*P(x,y,t);

    double W = 0.0;
    for (unsigned int j=0; j<setting.Lo; j++)
    {
        double _delta = delta(sn, setting.xi[j], 0);
        if (checkDelta(_delta))
        {
            double vi = 0.0;
            for (unsigned int i=0; i<setting.Lc; i++)
            {
                vi += setting.k[i][j] * P(setting.eta[i].x, setting.eta[i].y, t);
            }
            W += vi * _delta;
        }
    }
    res += W;

    return res;
}

double CProblem2Backward2D::g1(const SpaceNodePDE &sn UNUSED_PARAM, const TimeNodePDE &tn UNUSED_PARAM) const
{
    return -lambda*(sn.y*sn.y+tn.t);
}

double CProblem2Backward2D::g2(const SpaceNodePDE &sn UNUSED_PARAM, const TimeNodePDE &tn UNUSED_PARAM) const
{
    return 2.0 - lambda*(1.0+sn.y*sn.y+tn.t);
}

double CProblem2Backward2D::g3(const SpaceNodePDE &sn UNUSED_PARAM, const TimeNodePDE &tn UNUSED_PARAM) const
{
    return -lambda*(sn.x*sn.x+tn.t);
}

double CProblem2Backward2D::g4(const SpaceNodePDE &sn UNUSED_PARAM, const TimeNodePDE &tn UNUSED_PARAM) const
{
    return 2.0 - lambda*(sn.x*sn.x+1.0+tn.t);
}

double CProblem2Backward2D::h(const SpaceNodePDE &sn) const
{
    double x = sn.x; unsigned int i = sn.i;
    double y = sn.y; unsigned int j = sn.j;

    return (x*x + y*y + 1.0) + 2.0 * mu[j][i] * (uT[j][i] - U[j][i]);
}

double CProblem2Backward2D::P(double x, double y, double t) const
{
    return x*x + y*y + t;
}

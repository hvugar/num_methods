#include "cproblem2forward2d.h"

double CProblem2Forward2D::initial(const SpaceNodePDE &sn) const
{
    double x = sn.x;
    double y = sn.y;
    return x*x + y*y;
}

double CProblem2Forward2D::boundary(const SpaceNodePDE &sn UNUSED_PARAM, const TimeNodePDE &tn UNUSED_PARAM, BoundaryType boundary UNUSED_PARAM) const
{
    return NAN;
}

double CProblem2Forward2D::f(const SpaceNodePDE &sn, const TimeNodePDE &tn) const
{
    double x = sn.x;
    double y = sn.y;
    double t = tn.t;

    double res = 1.0 - 4.0*a*a + lambda0*(U(x,y,t) - theta);

//    std::vector<ControlDeltaNode> cndeltaNodes;
//    for (unsigned int i=0; i<setting.Lc; i++) extendContrlDeltaPoint(setting.eta[i], cndeltaNodes, i);

//    double W = 0.0;
//    for (unsigned int cni=0; cni<cndeltaNodes.size(); cni++)
//    {
//        const ControlDeltaNode &cn = cndeltaNodes.at(cni);
//        if (sn.i == cn.n && sn.j == cn.m)
//        {
//            double vi = 0.0;
//            for (unsigned int j=0; j<setting.Lo; j++)
//            {
//                vi += setting.k[cn.i][j] * (U(setting.xi[j].x, setting.xi[j].y, t) - setting.z[cn.i][j]);
//            }
//            W += vi * cn.w;
//        }
//    }
//    res -= W;
//    cndeltaNodes.clear();

    double W = 0.0;
    for (unsigned int i=0; i<setting.Lc; i++)
    {
        double _delta = delta(sn, setting.eta[i], i, 0);
        if (checkDelta(_delta))
        {
            double vi = 0.0;
            for (unsigned int j=0; j<setting.Lo; j++)
            {
                vi += setting.k[i][j] * ( U(setting.xi[j].x, setting.xi[j].y, t) - setting.z[i][j]);
            }
            W += vi * _delta;
        }
    }
    res -= W;

    return res;
}

double CProblem2Forward2D::g1(const SpaceNodePDE &sn, const TimeNodePDE &tn) const
{
    double y = sn.y;
    double t = tn.t;
    return lambda*(y*y+t - theta);
}

double CProblem2Forward2D::g2(const SpaceNodePDE &sn, const TimeNodePDE &tn) const
{
    double y = sn.y;
    double t = tn.t;
    //return 2.0 + lambda*(U(1.0, y, t) - theta);
    return 2.0 + lambda*(1.0 + y*y + t - theta);
}

double CProblem2Forward2D::g3(const SpaceNodePDE &sn, const TimeNodePDE &tn) const
{
    double x = sn.x;
    double t = tn.t;
    return lambda*(x*x + t - theta);
}

double CProblem2Forward2D::g4(const SpaceNodePDE &sn, const TimeNodePDE &tn) const
{
    double x = sn.x;
    double t = tn.t;
    return 2.0 + lambda*(1.0 + x*x + t - theta);
}

double CProblem2Forward2D::U(double x, double y, double t) const
{
    return x*x + y*y + t;
}

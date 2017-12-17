#include "cproblem2forward2d.h"

void CProblem2Forward2D::Main(int argc UNUSED_PARAM, char *argv[] UNUSED_PARAM)
{
    CProblem2Forward2D cpfp2d;
    cpfp2d.setTimeDimension(Dimension(0.01, 0,  100));
    cpfp2d.addSpaceDimension(Dimension(0.01, 0, 100));
    cpfp2d.addSpaceDimension(Dimension(0.01, 0, 100));
    cpfp2d.a = 1.0;
    cpfp2d.lambda0 = 0.01;
    cpfp2d.lambda = 0.1;
    cpfp2d.theta = 10.0;

    Parameter setting;
    setting.Lc = 4;
    setting.Lo = 4;

    setting.eta.resize(setting.Lc);
    //setting.eta[0].x = 0.50; setting.eta[0].y = 0.50;

    setting.eta.resize(setting.Lc);
    setting.eta[0].x = 0.3345; setting.eta[0].y = 0.3055;
    setting.eta[1].x = 0.3314; setting.eta[1].y = 0.6041;
    setting.eta[2].x = 0.6625; setting.eta[2].y = 0.6555;
    setting.eta[3].x = 0.6684; setting.eta[3].y = 0.3514;

    setting.xi.resize(setting.Lo);
    setting.xi[0].x = 0.25;  setting.xi[0].y = 0.25;
    setting.xi[1].x = 0.25;  setting.xi[1].y = 0.75;
    setting.xi[2].x = 0.75;  setting.xi[2].y = 0.75;
    setting.xi[3].x = 0.75;  setting.xi[3].y = 0.25;

    setting.k.resize(setting.Lc, setting.Lo);
    setting.z.resize(setting.Lc, setting.Lo);
    for (unsigned int i=0; i<setting.Lc; i++)
    {
        for (unsigned int j=0; j<setting.Lo; j++)
        {
            setting.k[i][j] = -0.1;
            setting.z[i][j] = +10.0;
        }
    }
    cpfp2d.setParamter(setting);

    DoubleMatrix u;
    vector<ExtendedSpaceNode2D> info;
    cpfp2d.calculateMVD(u,info,false);
    IPrinter::printSeperatorLine();
    IPrinter::printMatrix(u);
    IPrinter::printSeperatorLine();
}

double CProblem2Forward2D::initial(const SpaceNodePDE &sn) const
{
    return sn.x*sn.x + sn.y*sn.y;
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
    for (unsigned int i=0; i<mParameter.Lc; i++)
    {
        double _delta = delta(sn, mParameter.eta[i], i, 0);
        if (checkDelta(_delta))
        {
            double vi = 0.0;
            for (unsigned int j=0; j<mParameter.Lo; j++)
            {
                vi += mParameter.k[i][j] * ( U(mParameter.xi[j].x, mParameter.xi[j].y, t) - mParameter.z[i][j]);
            }
            W += vi * _delta;
        }
    }
    res -= W;

    return res;
}

double CProblem2Forward2D::g1(const SpaceNodePDE &sn, const TimeNodePDE &tn) const
{
    return lambda*(sn.y*sn.y + tn.t - theta);
}

double CProblem2Forward2D::g2(const SpaceNodePDE &sn, const TimeNodePDE &tn) const
{
    return 2.0 + lambda*(1.0 + sn.y*sn.y + tn.t - theta);
}

double CProblem2Forward2D::g3(const SpaceNodePDE &sn, const TimeNodePDE &tn) const
{
    return lambda*(sn.x*sn.x + tn.t - theta);
}

double CProblem2Forward2D::g4(const SpaceNodePDE &sn, const TimeNodePDE &tn) const
{
    return 2.0 + lambda*(sn.x*sn.x + 1.0 + tn.t - theta);
}

double CProblem2Forward2D::U(double x, double y, double t) const
{
    return x*x + y*y + t;
}

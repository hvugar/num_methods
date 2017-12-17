#include "cproblem2backward2d.h"

void CProblem2Backward2D::Main(int argc UNUSED_PARAM, char *argv[] UNUSED_PARAM)
{
    CProblem2Backward2D cpbp2d;
    cpbp2d.setTimeDimension(Dimension(0.01, 0,  100));
    cpbp2d.addSpaceDimension(Dimension(0.01, 0, 100));
    cpbp2d.addSpaceDimension(Dimension(0.01, 0, 100));
    cpbp2d.a = 1.0;
    cpbp2d.lambda0 = 0.01;
    cpbp2d.lambda = 0.1;
    cpbp2d.theta = 10.0;

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
    cpbp2d.setParamter(setting);

    DoubleMatrix u;
    vector<ExtendedSpaceNode2D> info;
    cpbp2d.calculateMVD(u,info,false);
    IPrinter::printSeperatorLine();
    IPrinter::printMatrix(u);
    IPrinter::printSeperatorLine();
}

double CProblem2Backward2D::initial(const SpaceNodePDE &sn) const
{
    return sn.x*sn.x + sn.y*sn.y + 1.0;
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
    for (unsigned int j=0; j<mParameter.Lo; j++)
    {
        double _delta = delta(sn, mParameter.xi[j], 0);
        if (checkDelta(_delta))
        {
            double vi = 0.0;
            for (unsigned int i=0; i<mParameter.Lc; i++)
            {
                vi += mParameter.k[i][j] * P(mParameter.eta[i].x, mParameter.eta[i].y, t);
            }
            W += vi * _delta;
        }
    }
    res += W;

    return res;
}

double CProblem2Backward2D::g1(const SpaceNodePDE &sn UNUSED_PARAM, const TimeNodePDE &tn UNUSED_PARAM) const
{
    return -lambda*(sn.y*sn.y + tn.t);
}

double CProblem2Backward2D::g2(const SpaceNodePDE &sn UNUSED_PARAM, const TimeNodePDE &tn UNUSED_PARAM) const
{
    return 2.0 - lambda*(1.0 + sn.y*sn.y + tn.t);
}

double CProblem2Backward2D::g3(const SpaceNodePDE &sn UNUSED_PARAM, const TimeNodePDE &tn UNUSED_PARAM) const
{
    return -lambda*(sn.x*sn.x + tn.t);
}

double CProblem2Backward2D::g4(const SpaceNodePDE &sn UNUSED_PARAM, const TimeNodePDE &tn UNUSED_PARAM) const
{
    return 2.0 - lambda*(sn.x*sn.x + 1.0 + tn.t);
}

double CProblem2Backward2D::P(double x, double y, double t) const
{
    return x*x + y*y + t;
}

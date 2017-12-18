#include "cproblem2forward2d.h"

void CProblem2Forward2D::Main(int argc UNUSED_PARAM, char *argv[] UNUSED_PARAM)
{
    CProblem2Forward2D cpfp2d;
    cpfp2d.setTimeDimension(Dimension(0.005, 0,  200));
    cpfp2d.addSpaceDimension(Dimension(0.005, 0, 200));
    cpfp2d.addSpaceDimension(Dimension(0.005, 0, 200));
    cpfp2d.setEquationParameters(1.0, 0.01, 0.1, 10.0);

    Parameter prm(4,4);

    prm.eta[0].setPoint(0.3345, 0.3055);
    prm.eta[1].setPoint(0.3314, 0.6041);
    prm.eta[2].setPoint(0.6625, 0.6555);
    prm.eta[3].setPoint(0.6684, 0.3514);

    prm.xi[0].setPoint(0.2534, 0.2534);
    prm.xi[1].setPoint(0.2578, 0.7544);
    prm.xi[2].setPoint(0.7545, 0.7556);
    prm.xi[3].setPoint(0.7543, 0.2589);

    for (unsigned int i=0; i<prm.Lc; i++)
    {
        for (unsigned int j=0; j<prm.Lo; j++)
        {
            prm.k[i][j] = -0.1;
            prm.z[i][j] = +10.0;
        }
    }
    cpfp2d.setParamter(prm);

    DoubleMatrix u;
    vector<ExtendedSpaceNode2D> info;
    clock_t t1 = clock();
    cpfp2d.calculateMVD(u,info,false);
    clock_t t2 = clock();
    printf("%f\n", (double)(t2-t1)/CLOCKS_PER_SEC);
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

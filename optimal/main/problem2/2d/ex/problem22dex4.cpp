#include "problem22dex4.h"

void Problem22DEx4::Main(int argc UNUSED_PARAM, char *argv[] UNUSED_PARAM)
{
    Problem22DEx4 p22Dex4;

    P2Setting setting;
    setting.Lc = 2;
    setting.Lo = 3;

    setting.eta.resize(setting.Lc);
    setting.eta[0].x = 0.20; setting.eta[0].y = 0.30;
    setting.eta[1].x = 0.75; setting.eta[1].y = 0.85;

    setting.xi.resize(setting.Lo);
    setting.xi[0].x = 0.15; setting.xi[0].y = 0.83;
    setting.xi[1].x = 0.74; setting.xi[1].y = 0.25;
    setting.xi[2].x = 0.50; setting.xi[1].y = 0.50;

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

    DoubleVector prm;
    setting.toVector(prm);
    DoubleVector g(prm.length());
    p22Dex4.gradient(prm,g);
}

//-------------------------------------------------------------------------------------------------------//

Problem22DEx4::Problem22DEx4() : AbstactProblem22D ()
{
    Dimension timeDimension = Dimension(0.01, 0, 100);
    Dimension spaceDimensionX = Dimension(0.01, 0, 100);
    Dimension spaceDimensionY = Dimension(0.01, 0, 100);

    forward = new Problem2Forward2DEx4;
    forward->setTimeDimension(timeDimension);
    forward->addSpaceDimension(spaceDimensionX);
    forward->addSpaceDimension(spaceDimensionY);
    //forward->setSettings(s);
    forward->fi = 0.0;

    backward = new Problem2Backward2DEx4;
    backward->setTimeDimension(timeDimension);
    backward->addSpaceDimension(spaceDimensionX);
    backward->addSpaceDimension(spaceDimensionY);
    //backward->setSettings(s);
    backward->p22dEx4 = this;

    forward->a = backward->a = 1.0;
    forward->lambda0 = backward->lambda0 = 0.01;
    forward->lambda = backward->lambda = 0.1;
    forward->theta = backward->theta = 10.0;

    U.resize(spaceDimensionY.sizeN()+1, spaceDimensionX.sizeN()+1, 0.0);
}

Problem22DEx4::~Problem22DEx4() {}

void Problem22DEx4::compareNandAGradients()
{

}

//-------------------------------------------------------------------------------------------------------//

double Problem2Forward2DEx4::initial(const SpaceNodePDE &) const
{
    return fi;
}

double Problem2Forward2DEx4::boundary(const SpaceNodePDE &, const TimeNodePDE &, BoundaryType) const { return NAN; }
double Problem2Forward2DEx4::f(const SpaceNodePDE &, const TimeNodePDE &) const { return 0.0; }
double Problem2Forward2DEx4::g1(const SpaceNodePDE &, const TimeNodePDE &) const { return 0.0; }
double Problem2Forward2DEx4::g2(const SpaceNodePDE &, const TimeNodePDE &) const { return 0.0; }
double Problem2Forward2DEx4::g3(const SpaceNodePDE &, const TimeNodePDE &) const { return 0.0; }
double Problem2Forward2DEx4::g4(const SpaceNodePDE &, const TimeNodePDE &) const { return 0.0; }

void Problem2Forward2DEx4::layerInfo(const DoubleMatrix &, unsigned int) const {}

//-------------------------------------------------------------------------------------------------------//

double Problem2Backward2DEx4::initial(const SpaceNodePDE & sn) const
{
    return -2.0 * p22dEx4->mu(sn.x, sn.y) * (u->at(sn.j,sn.i) - p22dEx4->U.at(sn.j,sn.i));
}

double Problem2Backward2DEx4::boundary(const SpaceNodePDE &, const TimeNodePDE &, BoundaryType) const { return NAN; }
double Problem2Backward2DEx4::f(const SpaceNodePDE &, const TimeNodePDE &) const { return 0.0; }
double Problem2Backward2DEx4::h(const SpaceNodePDE &) const { return 0.0; }
double Problem2Backward2DEx4::g1(const SpaceNodePDE &, const TimeNodePDE &) const { return 0.0; }
double Problem2Backward2DEx4::g2(const SpaceNodePDE &, const TimeNodePDE &) const { return 0.0; }
double Problem2Backward2DEx4::g3(const SpaceNodePDE &, const TimeNodePDE &) const { return 0.0; }
double Problem2Backward2DEx4::g4(const SpaceNodePDE &, const TimeNodePDE &) const { return 0.0; }

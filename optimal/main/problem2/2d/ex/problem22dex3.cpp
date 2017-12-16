#include "problem22dex3.h"
#include <QPixmap>

void Problem22DEx3::Main(int argc UNUSED_PARAM, char *argv[] UNUSED_PARAM)
{
    P2Setting s;
    //s.a = 0.1;
    //s.lambda = 0.1;
    //s.lambda0 = 0.1;
    //s.theta = 10.0;
    s.Lc = 1;
    s.Lo = 4;

    s.eta.resize(s.Lc);
    s.eta[0].x = 0.50; s.eta[0].y = 0.50;

    s.xi.resize(s.Lo);
    s.xi[0].x = 0.25;  s.xi[0].y = 0.25;
    s.xi[1].x = 0.25;  s.xi[1].y = 0.75;
    s.xi[2].x = 0.75;  s.xi[2].y = 0.75;
    s.xi[3].x = 0.75;  s.xi[3].y = 0.25;

    s.k.resize(s.Lc, s.Lo);
    s.z.resize(s.Lc, s.Lo);
    for (unsigned int i=0; i<s.Lc; i++)
    {
        for (unsigned int j=0; j<s.Lo; j++)
        {
            s.k[i][j] = -0.1;
            s.z[i][j] = +10.0;
        }
    }

    Dimension timeDimension = Dimension(0.001, 0, 1000);
    Dimension spaceDimensionX = Dimension(0.001, 0, 1000);
    Dimension spaceDimensionY = Dimension(0.001, 0, 1000);

    Problem2Forward2DEx3 *forward = new Problem2Forward2DEx3;
    forward->setSetting(s);
    forward->setTimeDimension(timeDimension);
    forward->addSpaceDimension(spaceDimensionX);
    forward->addSpaceDimension(spaceDimensionY);

    Problem2Backward2DEx3 *backward = NULL;// = new Problem2Backward2DEx3;
    backward->setSetting(s);
    backward->setTimeDimension(timeDimension);
    backward->addSpaceDimension(spaceDimensionX);
    backward->addSpaceDimension(spaceDimensionY);

    Problem22DEx3 p22dEx3;
    //p22dEx3.setForward(new Problem2Forward2DEx3);
    //p22dEx3.setBackward(new Problem2Backward2DEx3);
    p22dEx3.setGridParameters(timeDimension, spaceDimensionX, spaceDimensionY);
    p22dEx3.setP2Setting(s);

    backward->p22dEx3 = &p22dEx3;

    p22dEx3.U.resize(spaceDimensionY.sizeN()+1, spaceDimensionX.sizeN()+1, 10.0);

    DoubleVector prm;
    s.toVector(prm);
    DoubleVector g(prm.length());
    p22dEx3.gradient(prm,g);
}

double Problem2Forward2DEx3::initial(const SpaceNodePDE &) const
{
    return 0.0;
}

double Problem2Backward2DEx3::initial(const SpaceNodePDE & sn) const
{
    const DoubleMatrix &U = p22dEx3->U;
    const DoubleMatrix &u = p22dEx3->U;
    return -2.0 * p22dEx3->mu(sn.x, sn.y) * (u[sn.j][sn.i]-u[sn.j][sn.i]);
}



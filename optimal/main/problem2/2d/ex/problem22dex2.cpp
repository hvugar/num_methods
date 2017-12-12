#include "problem22dex2.h"
#include <QPixmap>

void Problem22DEx2::Main(int argc UNUSED_PARAM, char *argv[] UNUSED_PARAM)
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

    Problem22DEx2 p22dEx2;
    p22dEx2.setForward(new Problem2Forward2DEx2);
    p22dEx2.setBackward(new Problem2Backward2DEx2);
    p22dEx2.setGridParameters(timeDimension, spaceDimensionX, spaceDimensionY);
    p22dEx2.setP2Setting(s);

    Problem2Forward2DEx2* forward = new Problem2Forward2DEx2;
    forward->setSettings(s);
    forward->setTimeDimension(timeDimension);
    forward->addSpaceDimension(spaceDimensionX);
    forward->addSpaceDimension(spaceDimensionY);
    DoubleMatrix u;

    vector<ExtendedSpaceNode2D> info;
    forward->calculateMVD(u, info);

    FILE *file1 = fopen("pic.txt", "w");
    IPrinter::printMatrix(u, u.rows(), u.cols(), NULL, file1);
    fclose(file1);
}

void Problem2Forward2DEx2::layerInfo(const DoubleMatrix &u UNUSED_PARAM, unsigned int layerNumber UNUSED_PARAM) const
{
    if (layerNumber == 1000)
    {
        QPixmap pixmap;
        visualizeMatrixHeat(u, 0.0, u.max(), pixmap, u.cols(), u.rows());
        pixmap.save(QString("pics/pic%1.png").arg(layerNumber), "PNG");
    }
    printf("Layer %u. Painted\n", layerNumber);
}

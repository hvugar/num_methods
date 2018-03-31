#include "problem22dex1.h"
#include <QPixmap>

void Problem22DEx1::Main(int argc UNUSED_PARAM, char *argv[] UNUSED_PARAM)
{
    Parameter s;
    s.Lc = 4;
    s.Lo = 5;

    s.eta.resize(s.Lc);
    s.eta[0].x = 0.3345; s.eta[0].y = 0.3055;
    s.eta[1].x = 0.3314; s.eta[1].y = 0.6041;
    s.eta[2].x = 0.6625; s.eta[2].y = 0.6555;
    s.eta[3].x = 0.6684; s.eta[3].y = 0.3514;

    s.xi.resize(s.Lo);
    s.xi[0].x = 0.25;  s.xi[0].y = 0.25;
    s.xi[1].x = 0.25;  s.xi[1].y = 0.75;
    s.xi[2].x = 0.75;  s.xi[2].y = 0.75;
    s.xi[3].x = 0.75;  s.xi[3].y = 0.25;
    s.xi[4].x = 0.50;  s.xi[4].y = 0.50;

//    s.eta.resize(s.Lc);
//    s.eta[0].x = 0.40; s.eta[0].y = 0.50; s.eta[0].i = 400; s.eta[0].j = 500;
//    s.eta[1].x = 0.60; s.eta[1].y = 0.50; s.eta[1].i = 600; s.eta[1].j = 500;

//    s.xi.resize(s.Lo);
//    s.xi[0].x = 0.25;  s.xi[0].y = 0.25;
//    s.xi[1].x = 0.25;  s.xi[1].y = 0.75;
//    s.xi[2].x = 0.75;  s.xi[2].y = 0.75;
//    s.xi[3].x = 0.75;  s.xi[3].y = 0.25;

    s.k.resize(s.Lc, s.Lo);
    s.z.resize(s.Lc, s.Lo);
    for (unsigned int i=0; i<s.Lc; i++)
    {
        for (unsigned int j=0; j<s.Lo; j++)
        {
            s.k[i][j] = -5.1;//fabs(sin((i+1)*(j+2)));
            s.z[i][j] = +10.0;
        }
    }

    Dimension timeDimension = Dimension(0.001, 0, 1000);
    Dimension spaceDimensionX = Dimension(0.001, 0, 1000);
    Dimension spaceDimensionY = Dimension(0.001, 0, 1000);

    Problem22DEx1 p22dEx1;
    //p22dEx1.setForward(new Problem2Forward2DEx1);
    //p22dEx1.setBackward(new Problem2Backward2DEx1);
    p22dEx1.setGridParameters(timeDimension, spaceDimensionX, spaceDimensionY);
    p22dEx1.setParameter(s);

    Problem2Forward2DEx1* forward = new Problem2Forward2DEx1;
    forward->setEquationParameters(0.1, 0.1, 0.1);
    forward->setEnvTemperature(10.0);
    forward->setParameter(s);
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

void Problem2Forward2DEx1::layerInfo(const DoubleMatrix &u UNUSED_PARAM, unsigned int layerNumber UNUSED_PARAM) const
{
    if (layerNumber == 1000)
    {
        QPixmap pixmap;
        visualizeMatrixHeat(u, 0.0, u.max(), pixmap, u.cols(), u.rows());
        pixmap.save(QString("pics/pic%1.png").arg(layerNumber), "PNG");
    }
    printf("Layer %u. Painted\n", layerNumber);
}

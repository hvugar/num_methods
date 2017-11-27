#include "problem22dex1.h"
#include <QPixmap>

void Problem22DEx1::Main(int argc, char *argv[])
{
    P2Setting s;
    s.a = 0.1;
    s.lambda = 0.1;
    s.lambda0 = 0.1;
    s.theta = 10.0;
    s.Lc = 4;
    s.Lo = 5;

    s.eta.resize(s.Lc);
    s.eta[0].x = 0.3345; s.eta[0].y = 0.305;
    s.eta[1].x = 0.3314; s.eta[1].y = 0.604;
    s.eta[2].x = 0.6625; s.eta[2].y = 0.655;
    s.eta[3].x = 0.6684; s.eta[3].y = 0.351;

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

    Dimension timeDimension = Dimension(0.01, 0, 100);
    Dimension spaceDimensionX = Dimension(0.001, 0, 1000);
    Dimension spaceDimensionY = Dimension(0.001, 0, 1000);

    Problem22DEx1 p22dEx1;
    p22dEx1.setForward(new Problem2Forward2DEx1);
    p22dEx1.setBackward(new Problem2Backward2DEx1);
    p22dEx1.setGridParameters(timeDimension, spaceDimensionX, spaceDimensionY);
    p22dEx1.setP2Setting(s);

    Problem2Forward2DEx1* forward = new Problem2Forward2DEx1;
    forward->setSettings(s);
    forward->setTimeDimension(timeDimension);
    forward->addSpaceDimension(spaceDimensionX);
    forward->addSpaceDimension(spaceDimensionY);
    std::vector<DoubleMatrix> u;
    forward->calculateMVD(u);

    FILE *file1 = fopen("pic.txt", "w");

    DoubleMatrix uT = u[u.size()-1];
    IPrinter::printMatrix(uT, uT.rows(), uT.cols(), NULL, file1);

    fclose(file1);

    //    p22dEx1.gradient();
}

void Problem2Forward2DEx1::layerInfo(const DoubleMatrix &u, unsigned int layerNumber) const
{
//    if (layerNumber == 100)
//    {
//        QPixmap pixmap;
//        visualizeMatrixHeat(u, -50.0, 50.0, pixmap, u.cols(), u.rows());
//        pixmap.save(QString("pics/pic%1.png").arg(layerNumber), "PNG");
//    }
}

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
    setting.xi[2].x = 0.50; setting.xi[2].y = 0.50;

    setting.k.resize(setting.Lc, setting.Lo);
    setting.z.resize(setting.Lc, setting.Lo);
    for (unsigned int i=0; i<setting.Lc; i++)
    {
        for (unsigned int j=0; j<setting.Lo; j++)
        {
            setting.k[i][j] = -fabs(sin(0.1*(i+1))+cos((j-1)*0.2));
            setting.z[i][j] = +10.0;//+Random::value(-2,2,5);
        }
    }

    p22Dex4.setP2Setting(setting);

    DoubleMatrix U;
    vector<ExtendedSpaceNode2D> info;
    p22Dex4.forward->calculateMVD(U,info);
    p22Dex4.U = U;

    setting.eta[0].x = 0.22; setting.eta[0].y = 0.39;
    setting.eta[1].x = 0.73; setting.eta[1].y = 0.86;

    setting.xi[0].x = 0.12; setting.xi[0].y = 0.80;
    setting.xi[1].x = 0.78; setting.xi[1].y = 0.27;
    setting.xi[2].x = 0.55; setting.xi[2].y = 0.45;

    for (unsigned int i=0; i<setting.Lc; i++)
    {
        for (unsigned int j=0; j<setting.Lo; j++)
        {
            setting.k[i][j] = +fabs(cos(0.1*(i+1))+sin((j-1)*0.2));
            setting.z[i][j] = +10.0+Random::value(-2,2,5);
        }
    }

//    QPixmap pixmap;
//    visualizeMatrixHeat(U, U.min(), U.max(), pixmap, U.cols(), U.rows());
//    pixmap.save(QString("image1.png"), "PNG");
//    return;

//    IPrinter::printSeperatorLine();
//    IPrinter::printMatrix(p22Dex4.U ,10,10);
//    IPrinter::printSeperatorLine();


    DoubleVector prm;
    setting.toVector(prm);
    DoubleVector g(prm.length());
    IPrinter::print(prm, prm.length(), 10, 6);

    p22Dex4.gradient(prm,g);

    g.L2Normalize();
    IPrinter::print(g, g.length(), 10, 6);

    DoubleVector ng(prm.length());
    ng.resize(prm.length(), 0.0);

    IGradient::Gradient(&p22Dex4, 0.001, prm, ng);

//    for (unsigned int i=0; i<setting.Lc; i++)
//    {
//        for (unsigned int j=0; j<setting.Lo; j++)
//        {
//            unsigned int index = i*setting.Lo + j;
//            DoubleVector x2 = prm; x2[index] += 0.0001; double f2 = p22Dex4.fx(x2);
//            DoubleVector x1 = prm; x1[index] -= 0.0001; double f1 = p22Dex4.fx(x1);
//            ng[index] = (f2 - f1)/0.0002;
//        }
//    }

//    for (unsigned int i=0; i<setting.Lc; i++)
//    {
//        for (unsigned int j=0; j<setting.Lo; j++)
//        {
//            unsigned int index = setting.Lc*setting.Lo + i*setting.Lo + j;
//            DoubleVector x2 = prm; x2[index] += 0.0001; double f2 = p22Dex4.fx(x2);
//            DoubleVector x1 = prm; x1[index] -= 0.0001; double f1 = p22Dex4.fx(x1);
//            ng[index] = (f2 - f1)/0.0002;
//        }
//    }

//    for (unsigned int i=0; i<setting.Lc; i++)
//    {
//        unsigned int indexX = 2*setting.Lc*setting.Lo + 2*i + 0;
//        DoubleVector x2 = prm; x2[indexX] += 0.0001; double f2 = p22Dex4.fx(x2);
//        DoubleVector x1 = prm; x1[indexX] -= 0.0001; double f1 = p22Dex4.fx(x1);
//        ng[indexX] = (f2 - f1)/0.0002;

//        unsigned int indexY = 2*setting.Lc*setting.Lo + 2*i + 1;
//        DoubleVector x4 = prm; x4[indexY] += 0.0001; double f4 = p22Dex4.fx(x4);
//        DoubleVector x3 = prm; x3[indexY] -= 0.0001; double f3 = p22Dex4.fx(x3);
//        ng[indexY] = (f4 - f3)/0.0002;
//    }

//    for (unsigned int j=0; j<setting.Lo; j++)
//    {
//        unsigned int indexX = 2*setting.Lc*setting.Lo + 2*setting.Lc + 2*j + 0;
//        DoubleVector x2 = prm; x2[indexX] += 0.0001; double f2 = p22Dex4.fx(x2);
//        DoubleVector x1 = prm; x1[indexX] -= 0.0001; double f1 = p22Dex4.fx(x1);
//        ng[indexX] = (f2 - f1)/0.0002;

//        unsigned int indexY = 2*setting.Lc*setting.Lo + 2*setting.Lc + 2*j + 1;
//        DoubleVector x4 = prm; x4[indexY] += 0.0001; double f4 = p22Dex4.fx(x4);
//        DoubleVector x3 = prm; x3[indexY] -= 0.0001; double f3 = p22Dex4.fx(x3);
//        ng[indexY] = (f4 - f3)/0.0002;
//    }

    ng.L2Normalize();
    IPrinter::print(ng, ng.length(), 10, 6);
}

//-------------------------------------------------------------------------------------------------------//

Problem22DEx4::Problem22DEx4() : AbstactProblem22D ()
{
    Dimension timeDimension   = Dimension(0.005, 0, 200);
    Dimension spaceDimensionX = Dimension(0.005, 0, 200);
    Dimension spaceDimensionY = Dimension(0.005, 0, 200);

    setGridParameters(timeDimension, spaceDimensionX, spaceDimensionY);

    U.resize(spaceDimensionY.sizeN()+1, spaceDimensionX.sizeN()+1, 0.0);
}

Problem22DEx4::~Problem22DEx4() {}

void Problem22DEx4::compareNandAGradients()
{

}


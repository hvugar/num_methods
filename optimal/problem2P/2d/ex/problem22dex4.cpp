#include "problem22dex4.h"

void Problem22DEx4::Main(int argc UNUSED_PARAM, char *argv[] UNUSED_PARAM)
{
    Problem22DEx4 p22Dex4;

    Parameter p(2, 3);

    p.eta[0].setPoint(0.20, 0.30);
    p.eta[1].setPoint(0.75, 0.85);

    p.xi[0].setPoint(0.15, 0.83);
    p.xi[1].setPoint(0.74, 0.25);
    p.xi[2].setPoint(0.50, 0.50);

    for (unsigned int i=0; i<p.Lc; i++)
    {
        for (unsigned int j=0; j<p.Lo; j++)
        {
            p.k[i][j] = -fabs(sin(0.1*(i+1))+cos((j-1)*0.2));
            p.z[i][j] = +10.0+cos(0.1*(i+1))+sin((j-1)*0.2);//+Random::value(-2,2,5);
        }
    }

    p22Dex4.setParameter(p);

    DoubleMatrix U;
    vector<ExtendedSpaceNode2D> info;
    p22Dex4.forward->calculateMVD(U,info);
    p22Dex4.U = U;

    p.eta[0].setPoint(0.22, 0.35);
    p.eta[1].setPoint(0.73, 0.86);

    p.xi[0].setPoint(0.12, 0.80);
    p.xi[1].setPoint(0.78, 0.27);
    p.xi[2].setPoint(0.55, 0.45);

    for (unsigned int i=0; i<p.Lc; i++)
    {
        for (unsigned int j=0; j<p.Lo; j++)
        {
            p.k[i][j] = +Random::value(0,1,5);
            p.z[i][j] = +10.0+Random::value(-2,2,5);
        }
    }

//    QPixmap pixmap;
//    visualizeMatrixHeat(U, U.min(), U.max(), pixmap, U.cols(), U.rows());
//    pixmap.save(QString("image1.png"), "PNG");
//    return;

//    IPrinter::printSeperatorLine();
//    IPrinter::printMatrix(p22Dex4.U ,10,10);
//    IPrinter::printSeperatorLine();


    DoubleVector pv;
    p.toVector(pv);
    DoubleVector ag(pv.length());
    IPrinter::print(pv, pv.length(), 10, 6);

    p22Dex4.gradient(pv,ag);

    //ag.L2Normalize();
    //IPrinter::print(ag, ag.length(), 10, 6);

    DoubleVector ng(pv.length());
    ng.resize(pv.length(), 0.0);

    IGradient::Gradient(&p22Dex4, 0.01, pv, ng);

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

    //ng.L2Normalize();
    //IPrinter::print(ng, ng.length(), 10, 6);

    //------------------------------------------------------//
    DoubleVector pk = pv.mid(0, p.Lc*p.Lo-1);
    DoubleVector ak = ag.mid(0, p.Lc*p.Lo-1); ak.L2Normalize();
    DoubleVector nk = ng.mid(0, p.Lc*p.Lo-1); nk.L2Normalize();

    IPrinter::print(pk);
    IPrinter::print(ak);
    IPrinter::print(nk);
    IPrinter::printSeperatorLine();

    //------------------------------------------------------//
    DoubleVector pz = pv.mid(p.Lc*p.Lo, 2*p.Lc*p.Lo-1);
    DoubleVector az = ag.mid(p.Lc*p.Lo, 2*p.Lc*p.Lo-1); az.L2Normalize();
    DoubleVector nz = ng.mid(p.Lc*p.Lo, 2*p.Lc*p.Lo-1); nz.L2Normalize();

    IPrinter::print(pz);
    IPrinter::print(az);
    IPrinter::print(nz);
    IPrinter::printSeperatorLine();

    //------------------------------------------------------//
    DoubleVector pe = pv.mid(2*p.Lc*p.Lo, 2*p.Lc*p.Lo+2*p.Lc-1);
    DoubleVector ae = ag.mid(2*p.Lc*p.Lo, 2*p.Lc*p.Lo+2*p.Lc-1); ae.L2Normalize();
    DoubleVector ne = ng.mid(2*p.Lc*p.Lo, 2*p.Lc*p.Lo+2*p.Lc-1); ne.L2Normalize();

    IPrinter::print(pe);
    IPrinter::print(ae);
    IPrinter::print(ne);
    IPrinter::printSeperatorLine();

    //------------------------------------------------------//
    DoubleVector px = pv.mid(2*p.Lc*p.Lo+2*p.Lc, 2*p.Lc*p.Lo+2*p.Lc+2*p.Lo-1);
    DoubleVector ax = ag.mid(2*p.Lc*p.Lo+2*p.Lc, 2*p.Lc*p.Lo+2*p.Lc+2*p.Lo-1); ax.L2Normalize();
    DoubleVector nx = ng.mid(2*p.Lc*p.Lo+2*p.Lc, 2*p.Lc*p.Lo+2*p.Lc+2*p.Lo-1); nx.L2Normalize();

    IPrinter::print(px);
    IPrinter::print(ax);
    IPrinter::print(nx);
    IPrinter::printSeperatorLine();
}

//-------------------------------------------------------------------------------------------------------//

Problem22DEx4::Problem22DEx4() : AbstactProblem22D ()
{
    Dimension timeDimension   = Dimension(0.01, 0, 100);
    Dimension spaceDimensionX = Dimension(0.01, 0, 100);
    Dimension spaceDimensionY = Dimension(0.01, 0, 100);

    setGridParameters(timeDimension, spaceDimensionX, spaceDimensionY);

    U.resize(spaceDimensionY.sizeN()+1, spaceDimensionX.sizeN()+1, 0.0);
}

Problem22DEx4::~Problem22DEx4() {}

void Problem22DEx4::compareNandAGradients()
{

}


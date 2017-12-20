#include "problem22dex5.h"

void Problem22DEx5::Main(int argc UNUSED_PARAM, char *argv[] UNUSED_PARAM)
{
    Problem22DEx5 p22d5;
    p22d5.espilon = 0.0;

    {
        Parameter rpm0(2, 1);

        rpm0.eta[0].setPoint(0.40, 0.60);
        rpm0.eta[1].setPoint(0.70, 0.40);

        rpm0.xi[0].setPoint(0.50, 0.30);
        //rpm0.xi[1].setPoint(0.80, 0.70);

        rpm0.k[0][0] = -2.0; //rpm0.k[0][1] = -4.0;
        rpm0.k[1][0] = -2.8; //rpm0.k[1][1] = +3.2;

        rpm0.z[0][0] = 10.5; //rpm0.z[0][1] = 11.4;
        rpm0.z[1][0] = 11.7; //rpm0.z[1][1] = 12.5;

        //p22d5.calculateU(rpm0);

        IPrinter::printSeperatorLine("U");
        IPrinter::printMatrix(p22d5.U);

        QPixmap px;
        visualizeMatrixHeat(p22d5.U, 0.1, p22d5.U.max(), px);
        px.save("U.png", "PNG");

        //-----------------------------------------------//
        IPrinter::printSeperatorLine();
        DoubleVector pv1;
        rpm0.toVector(pv1);
        IPrinter::print(pv1, pv1.length(), 10, 6);
        IPrinter::printSeperatorLine();
        //-----------------------------------------------//
    }

    Parameter rpm(2, 1);
    rpm.eta[0].setPoint(0.40, 0.60);
    rpm.eta[1].setPoint(0.70, 0.40);

    rpm.xi[0].setPoint(0.50, 0.30);
    //rpm.xi[1].setPoint(0.80, 0.70);

    rpm.k[0][0] = -3.0; //rpm.k[0][1] = -5.0;
    rpm.k[1][0] = -2.8; //rpm.k[1][1] = +3.2;

    rpm.z[0][0] = 10.5; //rpm.z[0][1] = 11.4;
    rpm.z[1][0] = 11.7; //rpm.z[1][1] = 12.5;

    p22d5.setParameter(rpm);

    IPrinter::printSeperatorLine();
    DoubleVector pv;
    rpm.toVector(pv);
    DoubleVector ag(pv.length());
    IPrinter::print(pv, pv.length(), 10, 6);
    IPrinter::printSeperatorLine();

    p22d5.gradient(pv,ag);

    double functional = p22d5.fx(pv);
    printf("Functional: %f\n", functional);

    DoubleVector ng(pv.length());
    ng.resize(pv.length(), 0.0);

    IGradient::Gradient(&p22d5, 0.01, pv, ng);

    //------------------------------------------------------//
    DoubleVector pk = pv.mid(0, rpm.Lc*rpm.Lo-1);
    DoubleVector ak = ag.mid(0, rpm.Lc*rpm.Lo-1);
    DoubleVector nk = ng.mid(0, rpm.Lc*rpm.Lo-1);

    IPrinter::print(pk);
    IPrinter::print(ak); ak.L2Normalize();
    IPrinter::print(nk); nk.L2Normalize();
    IPrinter::print(ak);
    IPrinter::print(nk);
    IPrinter::printSeperatorLine();

    //------------------------------------------------------//
    DoubleVector pz = pv.mid(rpm.Lc*rpm.Lo, 2*rpm.Lc*rpm.Lo-1);
    DoubleVector az = ag.mid(rpm.Lc*rpm.Lo, 2*rpm.Lc*rpm.Lo-1);
    DoubleVector nz = ng.mid(rpm.Lc*rpm.Lo, 2*rpm.Lc*rpm.Lo-1);

    IPrinter::print(pz);
    IPrinter::print(az); az.L2Normalize();
    IPrinter::print(nz); nz.L2Normalize();
    IPrinter::print(az);
    IPrinter::print(nz);
    IPrinter::printSeperatorLine();

    //------------------------------------------------------//
    DoubleVector pe = pv.mid(2*rpm.Lc*rpm.Lo, 2*rpm.Lc*rpm.Lo+2*rpm.Lc-1);
    DoubleVector ae = ag.mid(2*rpm.Lc*rpm.Lo, 2*rpm.Lc*rpm.Lo+2*rpm.Lc-1);
    DoubleVector ne = ng.mid(2*rpm.Lc*rpm.Lo, 2*rpm.Lc*rpm.Lo+2*rpm.Lc-1);

    IPrinter::print(pe);
    IPrinter::print(ae); ae.L2Normalize();
    IPrinter::print(ne); ne.L2Normalize();
    IPrinter::print(ae);
    IPrinter::print(ne);
    IPrinter::printSeperatorLine();

    //------------------------------------------------------//
    DoubleVector px = pv.mid(2*rpm.Lc*rpm.Lo+2*rpm.Lc, 2*rpm.Lc*rpm.Lo+2*rpm.Lc+2*rpm.Lo-1);
    DoubleVector ax = ag.mid(2*rpm.Lc*rpm.Lo+2*rpm.Lc, 2*rpm.Lc*rpm.Lo+2*rpm.Lc+2*rpm.Lo-1);
    DoubleVector nx = ng.mid(2*rpm.Lc*rpm.Lo+2*rpm.Lc, 2*rpm.Lc*rpm.Lo+2*rpm.Lc+2*rpm.Lo-1);

    IPrinter::print(px);
    IPrinter::print(ax); ax.L2Normalize();
    IPrinter::print(nx); nx.L2Normalize();
    IPrinter::print(ax);
    IPrinter::print(nx);
    IPrinter::printSeperatorLine();
}

Problem22DEx5::Problem22DEx5() : AbstactProblem22D()
{
    forward->setEquationParameters(1.0, 0.001, 1.0, 6.3);
    backward->setEquationParameters(1.0, 0.001, 1.0, 6.3);
    forward->fi = 0.1;

    Dimension time = Dimension(0.005, 0, 200);
    Dimension dimX = Dimension(0.01, 0, 100);
    Dimension dimY = Dimension(0.01, 0, 100);

    setGridParameters(time, dimX, dimY);

    U.resize(dimY.sizeN()+1, dimX.sizeN()+1, 10.0);
}

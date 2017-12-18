#include "problem22dex5.h"

void Problem22DEx5::Main(int argc UNUSED_PARAM, char *argv[] UNUSED_PARAM)
{
    Problem22DEx5 p22d5;
    p22d5.espilon = 0.0;

    {
        Parameter rpm0(2, 2);

        rpm0.eta[0].setPoint(0.20, 0.80);
        rpm0.eta[1].setPoint(0.70, 0.40);

        rpm0.xi[0].setPoint(0.30, 0.20);
        rpm0.xi[1].setPoint(0.80, 0.70);

        rpm0.k[0][0] = -2.0; rpm0.k[0][1] = -4.0;
        rpm0.k[1][0] = -2.8; rpm0.k[1][1] = +3.2;

        rpm0.z[0][0] = 10.5; rpm0.z[0][1] = 11.4;
        //rpm0.z[1][0] = 9.7;  rpm0.z[1][1] = 8.5;

        p22d5.calculateU(rpm0);

        IPrinter::printSeperatorLine("U");
        IPrinter::printMatrix(p22d5.U);

        QPixmap px;
        visualizeMatrixHeat(p22d5.U,0.0, p22d5.U.max(), px);
        px.save("U.png", "PNG");


        //-----------------------------------------------//
        IPrinter::printSeperatorLine();
        DoubleVector pv1;
        rpm0.toVector(pv1);
        IPrinter::print(pv1, pv1.length(), 10, 6);
        IPrinter::printSeperatorLine();
        //-----------------------------------------------//
    }

    Parameter rpm(2, 2);
    rpm.eta[0].setPoint(0.22, 0.85);
    rpm.eta[1].setPoint(0.73, 0.45);

    rpm.xi[0].setPoint(0.28, 0.25);
    rpm.xi[1].setPoint(0.78, 0.77);

    rpm.k[0][0] = -1.0; rpm.k[0][1] = -2.6;
    rpm.k[1][0] = -1.8; rpm.k[1][1] = -6.7;

    rpm.z[0][0] = 11.1; rpm.z[0][1] = 14.4;
    rpm.z[1][0] = 10.7; rpm.z[1][1] = 12.5;

    p22d5.setParameter(rpm);

    IPrinter::printSeperatorLine();
    DoubleVector pv;
    rpm.toVector(pv);
    DoubleVector ag(pv.length());
    IPrinter::print(pv, pv.length(), 10, 6);
    IPrinter::printSeperatorLine();puts("OK");

    p22d5.gradient(pv,ag);
    printf("Functional: %f\n", p22d5.fx(pv));

    DoubleVector ng(pv.length());
    ng.resize(pv.length(), 0.0);

    IGradient::Gradient(&p22d5, 0.001, pv, ng);

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
    forward->setEquationParameters(1.0, 0.001, 0.1, 6.3);
    backward->setEquationParameters(1.0, 0.001, 0.1, 6.3);
    forward->fi = 0.1;

    Dimension time = Dimension(0.01, 0, 100);
    Dimension dimX = Dimension(0.01, 0, 100);
    Dimension dimY = Dimension(0.01, 0, 100);

    setGridParameters(time, dimX, dimY);

    U.resize(dimY.sizeN()+1, dimX.sizeN()+1, 10.0);
}

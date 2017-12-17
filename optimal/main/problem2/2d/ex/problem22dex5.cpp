#include "problem22dex5.h"

void Problem22DEx5::Main(int argc UNUSED_PARAM, char *argv[] UNUSED_PARAM)
{
    Problem22DEx5 p22d5;

    Parameter p(2, 2);

    p.eta[0].setPoint(0.20, 0.80);
    p.eta[1].setPoint(0.70, 0.40);

    p.xi[0].setPoint(0.30, 0.20);
    p.xi[1].setPoint(0.80, 0.70);

    p.k[0][0] = -2.0; p.k[0][1] = -4.0;
    p.k[1][0] = -2.8; p.k[1][1] = -3.2;

    p.z[0][0] = 10.5; p.z[0][1] = 11.4;
    p.z[1][0] = 9.7;  p.z[1][1] = 8.5;

    IPrinter::printSeperatorLine();
    DoubleVector pv1;
    p.toVector(pv1);
    IPrinter::print(pv1, pv1.length(), 10, 6);
    IPrinter::printSeperatorLine();

    p22d5.setParameter(p);

    DoubleMatrix u;
    vector<ExtendedSpaceNode2D> info;
    p22d5.forward->calculateMVD(u,info,false);
    IPrinter::printMatrix(u);
    p22d5.U = u;

    //p.eta[0].setPoint(0.22, 0.85);
    //p.eta[1].setPoint(0.73, 0.45);

    //p.xi[0].setPoint(0.28, 0.25);
    //p.xi[1].setPoint(0.78, 0.77);

    //p.k[0][0] = -3.0; p.k[0][1] = -2.6;
    //p.k[1][0] = -1.8; p.k[1][1] = -6.7;

    //p.z[0][0] = 11.1; p.z[0][1] = 14.4;
    //p.z[1][0] = 10.7; p.z[1][1] = 12.5;

    p22d5.setParameter(p);

    IPrinter::printSeperatorLine();
    DoubleVector pv;
    p.toVector(pv);
    DoubleVector ag(pv.length());
    IPrinter::print(pv, pv.length(), 10, 6);
    IPrinter::printSeperatorLine();

    p22d5.gradient(pv,ag);
    printf("Functional: %f\n", p22d5.fx(pv));

    DoubleVector ng(pv.length());
    ng.resize(pv.length(), 0.0);

    IGradient::Gradient(&p22d5, 0.001, pv, ng);

    //------------------------------------------------------//
    DoubleVector pk = pv.mid(0, p.Lc*p.Lo-1);
    DoubleVector ak = ag.mid(0, p.Lc*p.Lo-1); //ak.L2Normalize();
    DoubleVector nk = ng.mid(0, p.Lc*p.Lo-1); //nk.L2Normalize();

    IPrinter::print(pk);
    IPrinter::print(ak);
    IPrinter::print(nk);
    IPrinter::printSeperatorLine();

    //------------------------------------------------------//
    DoubleVector pz = pv.mid(p.Lc*p.Lo, 2*p.Lc*p.Lo-1);
    DoubleVector az = ag.mid(p.Lc*p.Lo, 2*p.Lc*p.Lo-1); //az.L2Normalize();
    DoubleVector nz = ng.mid(p.Lc*p.Lo, 2*p.Lc*p.Lo-1); //nz.L2Normalize();

    IPrinter::print(pz);
    IPrinter::print(az);
    IPrinter::print(nz);
    IPrinter::printSeperatorLine();

    //------------------------------------------------------//
    DoubleVector pe = pv.mid(2*p.Lc*p.Lo, 2*p.Lc*p.Lo+2*p.Lc-1);
    DoubleVector ae = ag.mid(2*p.Lc*p.Lo, 2*p.Lc*p.Lo+2*p.Lc-1); //ae.L2Normalize();
    DoubleVector ne = ng.mid(2*p.Lc*p.Lo, 2*p.Lc*p.Lo+2*p.Lc-1); //ne.L2Normalize();

    IPrinter::print(pe);
    IPrinter::print(ae);
    IPrinter::print(ne);
    IPrinter::printSeperatorLine();

    //------------------------------------------------------//
    DoubleVector px = pv.mid(2*p.Lc*p.Lo+2*p.Lc, 2*p.Lc*p.Lo+2*p.Lc+2*p.Lo-1);
    DoubleVector ax = ag.mid(2*p.Lc*p.Lo+2*p.Lc, 2*p.Lc*p.Lo+2*p.Lc+2*p.Lo-1); //ax.L2Normalize();
    DoubleVector nx = ng.mid(2*p.Lc*p.Lo+2*p.Lc, 2*p.Lc*p.Lo+2*p.Lc+2*p.Lo-1); //nx.L2Normalize();

    IPrinter::print(px);
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

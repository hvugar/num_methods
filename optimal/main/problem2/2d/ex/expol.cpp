#include "expol.h"

void ExpOptimalLetters::Main(int argc UNUSED_PARAM, char *argv[] UNUSED_PARAM)
{
    //Table1Y2();
    //Table2Y1();
    Table3Y2();
    //test();
}

void ExpOptimalLetters::Table1Y1()
{
    JFunctional jfunc;
    unsigned int Lc = 2;
    unsigned int Lo = 2;

    jfunc.setGridParameters(Dimension(0.005, 0, 200), Dimension(0.010, 0, 100), Dimension(0.010, 0, 100));
    jfunc.U.resize(101, 101, 10.0);

    DoubleVector fis; fis << +0.2 << +0.3 << +0.5;
    DoubleVector p_fis(fis.length(), 1.0/fis.length());

    DoubleVector thetas; thetas << +6.3 << +6.4 << +6.5;
    DoubleVector p_thetas(thetas.length(), 1.0/thetas.length());

    jfunc.setInitTemperatures(fis, p_fis);
    jfunc.setEnvrTemperatures(thetas, p_thetas);

    jfunc.setEquationParameters(1.0, 0.01, 0.01);
    jfunc.setEpsilon(0.0);
    jfunc.setPenaltyCoefficient(1.0);
    jfunc.setPenaltyLimits(DoubleVector(Lc, -5.0), DoubleVector(Lc, +20.0));

    Parameter prm0(Lc, Lo);
    prm0.k[0][0] = -1.100; prm0.k[1][0] = -1.128;
    prm0.k[0][1] = -1.110; prm0.k[1][1] = -1.104;
    prm0.z[0][0] = +10.50; prm0.z[1][0] = +10.70;
    prm0.z[0][1] = +12.40; prm0.z[1][1] = +10.50;
    prm0.eta[0].setPoint(0.3000, 0.6000);
    prm0.eta[1].setPoint(0.6000, 0.2000);
    prm0.xi[0].setPoint(0.5000, 0.8000);
    prm0.xi[1].setPoint(0.2500, 0.3000);
    jfunc.setParameter0(prm0);

    Parameter prm(Lc, Lo);
    prm.k[0][0] = -1.86; prm.k[1][0] = -1.98;
    prm.k[0][1] = -1.93; prm.k[1][1] = -1.82;
    prm.z[0][0] = +10.5; prm.z[1][0] = +11.7;
    prm.z[0][1] = +11.4; prm.z[1][1] = +10.5;
    prm.eta[0].setPoint(0.3164,0.6854);
    prm.eta[1].setPoint(0.5847,0.3524);
    prm.xi[0].setPoint(0.7341,0.8248);
    prm.xi[1].setPoint(0.2116,0.2329);
    jfunc.setParameter(prm);

    DoubleMatrix u;
    vector<ExtendedSpaceNode2D> info;
    jfunc.setIntTemperature(fis[0]);
    jfunc.setEnvTemperature(thetas[0]);
    jfunc.forward->calculateMVD(u, info, false);
    QPixmap pxm;
    visualizeMatrixHeat(u, u.min(), u.max(), pxm, 101, 101);
    pxm.save("imageU.png", "PNG");
    DoubleVector v;
    prm.toVector(v);
    printf("%f %f %f\n", u.min(), u.max(), jfunc.fx(v));
    //return;

    IPrinter::printSeperatorLine();
    DoubleVector pv;
    prm.toVector(pv);
    DoubleVector ag(pv.length());
    IPrinter::print(pv, pv.length(), 6, 4);
    IPrinter::printSeperatorLine();

    puts("Calculating gradients....");
    jfunc.gradient(pv,ag);
    puts("Gradients are calculated.");

    double functional = jfunc.fx(pv);
    printf("Functional: %f\n", functional);

    DoubleVector ng(pv.length());
    ng.resize(pv.length(), 0.0);

    puts("Calculating numerical gradients....");
    IGradient::Gradient(&jfunc, 0.010, pv, ng, 0*prm.Lc*prm.Lo,          1*prm.Lc*prm.Lo-1);
    IGradient::Gradient(&jfunc, 0.010, pv, ng, 1*prm.Lc*prm.Lo,          2*prm.Lc*prm.Lo-1);
    IGradient::Gradient(&jfunc, 0.010, pv, ng, 2*prm.Lc*prm.Lo+0*Lc,     2*prm.Lc*prm.Lo+2*prm.Lc-1);
    IGradient::Gradient(&jfunc, 0.010, pv, ng, 2*prm.Lc*prm.Lo+2*prm.Lc, 2*prm.Lc*prm.Lo+2*prm.Lc+2*prm.Lo-1);
    puts("Numerical gradients are calculated.");

    //k------------------------------------------------------//
    DoubleVector pk = pv.mid(0, prm.Lc*prm.Lo-1);
    DoubleVector ak = ag.mid(0, prm.Lc*prm.Lo-1);
    DoubleVector nk = ng.mid(0, prm.Lc*prm.Lo-1);

    IPrinter::print(pk,pk.length(),14,4);
    IPrinter::print(ak,ak.length(),14,4); ak.L2Normalize();
    IPrinter::print(nk,nk.length(),14,4); nk.L2Normalize();
    IPrinter::print(ak,ak.length(),14,4);
    IPrinter::print(nk,nk.length(),14,4);
    IPrinter::printSeperatorLine();

    //z------------------------------------------------------//
    DoubleVector pz = pv.mid(prm.Lc*prm.Lo, 2*prm.Lc*prm.Lo-1);
    DoubleVector az = ag.mid(prm.Lc*prm.Lo, 2*prm.Lc*prm.Lo-1);
    DoubleVector nz = ng.mid(prm.Lc*prm.Lo, 2*prm.Lc*prm.Lo-1);

    IPrinter::print(pz,pz.length(),14,4);
    IPrinter::print(az,az.length(),14,4); az.L2Normalize();
    IPrinter::print(nz,nz.length(),14,4); nz.L2Normalize();
    IPrinter::print(az,az.length(),14,4);
    IPrinter::print(nz,nz.length(),14,4);
    IPrinter::printSeperatorLine();

    //eta------------------------------------------------------//
    DoubleVector pe = pv.mid(2*prm.Lc*prm.Lo, 2*prm.Lc*prm.Lo+2*prm.Lc-1);
    DoubleVector ae = ag.mid(2*prm.Lc*prm.Lo, 2*prm.Lc*prm.Lo+2*prm.Lc-1);
    DoubleVector ne = ng.mid(2*prm.Lc*prm.Lo, 2*prm.Lc*prm.Lo+2*prm.Lc-1);

    IPrinter::print(pe,pe.length(),14,4);
    IPrinter::print(ae,ae.length(),14,4); ae.L2Normalize();
    IPrinter::print(ne,ne.length(),14,4); ne.L2Normalize();
    IPrinter::print(ae,ae.length(),14,4);
    IPrinter::print(ne,ne.length(),14,4);
    IPrinter::printSeperatorLine();

    //xi------------------------------------------------------//
    DoubleVector px = pv.mid(2*prm.Lc*prm.Lo+2*prm.Lc, 2*prm.Lc*prm.Lo+2*prm.Lc+2*prm.Lo-1);
    DoubleVector ax = ag.mid(2*prm.Lc*prm.Lo+2*prm.Lc, 2*prm.Lc*prm.Lo+2*prm.Lc+2*prm.Lo-1);
    DoubleVector nx = ng.mid(2*prm.Lc*prm.Lo+2*prm.Lc, 2*prm.Lc*prm.Lo+2*prm.Lc+2*prm.Lo-1);

    IPrinter::print(px,px.length(),14,4);
    IPrinter::print(ax,ax.length(),14,4); ax.L2Normalize();
    IPrinter::print(nx,nx.length(),14,4); nx.L2Normalize();
    IPrinter::print(ax,ax.length(),14,4);
    IPrinter::print(nx,nx.length(),14,4);
    IPrinter::printSeperatorLine();
}

void ExpOptimalLetters::Table1Y2()
{
    JFunctional jfunc;
    unsigned int Lc = 2;
    unsigned int Lo = 2;

    jfunc.setGridParameters(Dimension(0.005, 0, 200), Dimension(0.010, 0, 100), Dimension(0.010, 0, 100));
    jfunc.U.resize(101, 101, 10.0);

    DoubleVector fis; fis << +0.2 << +0.3 << +0.5;
    DoubleVector p_fis(fis.length(), 1.0/fis.length());

    DoubleVector thetas; thetas << +6.3 << +6.4 << +6.5;
    DoubleVector p_thetas(thetas.length(), 1.0/thetas.length());

    jfunc.setInitTemperatures(fis, p_fis);
    jfunc.setEnvrTemperatures(thetas, p_thetas);

    jfunc.setEquationParameters(1.0, 0.01, 0.01);
    jfunc.setEpsilon(0.0);
    jfunc.setPenaltyCoefficient(1.0);
    jfunc.setPenaltyLimits(DoubleVector(Lc, -5.0), DoubleVector(Lc, +20.0));

    Parameter prm0(Lc, Lo);
    prm0.k[0][0] = -1.100; prm0.k[1][0] = -1.128;
    prm0.k[0][1] = -1.110; prm0.k[1][1] = -1.104;
    prm0.z[0][0] = +10.50; prm0.z[1][0] = +10.70;
    prm0.z[0][1] = +12.40; prm0.z[1][1] = +10.50;
    prm0.eta[0].setPoint(0.3000, 0.6000);
    prm0.eta[1].setPoint(0.6000, 0.2000);
    prm0.xi[0].setPoint(0.5000, 0.8000);
    prm0.xi[1].setPoint(0.2500, 0.3000);
    jfunc.setParameter0(prm0);

    Parameter prm(Lc, Lo);
    prm.k[0][0] = -1.30; prm.k[1][0] = -1.15;
    prm.k[0][1] = -1.65; prm.k[1][1] = -1.34;
    prm.z[0][0] = +12.54; prm.z[1][0] = +12.58;
    prm.z[0][1] = +14.44; prm.z[1][1] = +13.48;
    prm.eta[0].setPoint(0.7574,0.6384);
    prm.eta[1].setPoint(0.3685,0.1275);
    prm.xi[0].setPoint(0.2514, 0.5254);
    prm.xi[1].setPoint(0.6147, 0.2487);

    jfunc.setParameter(prm);

    DoubleMatrix u;
    vector<ExtendedSpaceNode2D> info;
    jfunc.setIntTemperature(fis[0]);
    jfunc.setEnvTemperature(thetas[0]);
    jfunc.forward->calculateMVD(u, info, false);
    QPixmap pxm;
    visualizeMatrixHeat(u, u.min(), u.max(), pxm, 101, 101);
    pxm.save("imageU.png", "PNG");
    DoubleVector v;
    prm.toVector(v);
    printf("%f %f %f\n", u.min(), u.max(), jfunc.fx(v));
    //return;

    IPrinter::printSeperatorLine();
    DoubleVector pv;
    prm.toVector(pv);
    DoubleVector ag(pv.length());
    IPrinter::print(pv, pv.length(), 6, 4);
    IPrinter::printSeperatorLine();

    puts("Calculating gradients....");
    jfunc.gradient(pv,ag);
    puts("Gradients are calculated.");

    double functional = jfunc.fx(pv);
    printf("Functional: %f\n", functional);

    DoubleVector ng(pv.length());
    ng.resize(pv.length(), 0.0);

    puts("Calculating numerical gradients....");
    IGradient::Gradient(&jfunc, 0.010, pv, ng, 0*prm.Lc*prm.Lo,          1*prm.Lc*prm.Lo-1);
    IGradient::Gradient(&jfunc, 0.010, pv, ng, 1*prm.Lc*prm.Lo,          2*prm.Lc*prm.Lo-1);
    IGradient::Gradient(&jfunc, 0.010, pv, ng, 2*prm.Lc*prm.Lo+0*Lc,     2*prm.Lc*prm.Lo+2*prm.Lc-1);
    IGradient::Gradient(&jfunc, 0.010, pv, ng, 2*prm.Lc*prm.Lo+2*prm.Lc, 2*prm.Lc*prm.Lo+2*prm.Lc+2*prm.Lo-1);
    puts("Numerical gradients are calculated.");

    //k------------------------------------------------------//
    DoubleVector pk = pv.mid(0, prm.Lc*prm.Lo-1);
    DoubleVector ak = ag.mid(0, prm.Lc*prm.Lo-1);
    DoubleVector nk = ng.mid(0, prm.Lc*prm.Lo-1);

    IPrinter::print(pk,pk.length(),14,4);
    IPrinter::print(ak,ak.length(),14,4); ak.L2Normalize();
    IPrinter::print(nk,nk.length(),14,4); nk.L2Normalize();
    IPrinter::print(ak,ak.length(),14,4);
    IPrinter::print(nk,nk.length(),14,4);
    IPrinter::printSeperatorLine();

    //z------------------------------------------------------//
    DoubleVector pz = pv.mid(prm.Lc*prm.Lo, 2*prm.Lc*prm.Lo-1);
    DoubleVector az = ag.mid(prm.Lc*prm.Lo, 2*prm.Lc*prm.Lo-1);
    DoubleVector nz = ng.mid(prm.Lc*prm.Lo, 2*prm.Lc*prm.Lo-1);

    IPrinter::print(pz,pz.length(),14,4);
    IPrinter::print(az,az.length(),14,4); az.L2Normalize();
    IPrinter::print(nz,nz.length(),14,4); nz.L2Normalize();
    IPrinter::print(az,az.length(),14,4);
    IPrinter::print(nz,nz.length(),14,4);
    IPrinter::printSeperatorLine();

    //eta------------------------------------------------------//
    DoubleVector pe = pv.mid(2*prm.Lc*prm.Lo, 2*prm.Lc*prm.Lo+2*prm.Lc-1);
    DoubleVector ae = ag.mid(2*prm.Lc*prm.Lo, 2*prm.Lc*prm.Lo+2*prm.Lc-1);
    DoubleVector ne = ng.mid(2*prm.Lc*prm.Lo, 2*prm.Lc*prm.Lo+2*prm.Lc-1);

    IPrinter::print(pe,pe.length(),14,4);
    IPrinter::print(ae,ae.length(),14,4); ae.L2Normalize();
    IPrinter::print(ne,ne.length(),14,4); ne.L2Normalize();
    IPrinter::print(ae,ae.length(),14,4);
    IPrinter::print(ne,ne.length(),14,4);
    IPrinter::printSeperatorLine();

    //xi------------------------------------------------------//
    DoubleVector px = pv.mid(2*prm.Lc*prm.Lo+2*prm.Lc, 2*prm.Lc*prm.Lo+2*prm.Lc+2*prm.Lo-1);
    DoubleVector ax = ag.mid(2*prm.Lc*prm.Lo+2*prm.Lc, 2*prm.Lc*prm.Lo+2*prm.Lc+2*prm.Lo-1);
    DoubleVector nx = ng.mid(2*prm.Lc*prm.Lo+2*prm.Lc, 2*prm.Lc*prm.Lo+2*prm.Lc+2*prm.Lo-1);

    IPrinter::print(px,px.length(),14,4);
    IPrinter::print(ax,ax.length(),14,4); ax.L2Normalize();
    IPrinter::print(nx,nx.length(),14,4); nx.L2Normalize();
    IPrinter::print(ax,ax.length(),14,4);
    IPrinter::print(nx,nx.length(),14,4);
    IPrinter::printSeperatorLine();
}

void ExpOptimalLetters::Table2Y1()
{
    JFunctional jfunc;
    unsigned int Lc = 2;
    unsigned int Lo = 2;

    jfunc.setGridParameters(Dimension(0.005, 0, 200), Dimension(0.010, 0, 100), Dimension(0.010, 0, 100));
    jfunc.U.resize(101, 101, 10.0);

    DoubleVector fis; fis << +0.2;// << +0.3 << +0.5;
    DoubleVector p_fis(fis.length(), 1.0/fis.length());

    DoubleVector thetas; thetas << +6.3;// << +6.4 << +6.5;
    DoubleVector p_thetas(thetas.length(), 1.0/thetas.length());

    jfunc.setInitTemperatures(fis, p_fis);
    jfunc.setEnvrTemperatures(thetas, p_thetas);

    jfunc.setEquationParameters(1.0, 0.01, 0.01);
    jfunc.setEpsilon(0.0);
    jfunc.setPenaltyCoefficient(20.0);
    jfunc.setPenaltyLimits(DoubleVector(Lc, -5.0), DoubleVector(Lc, +20.0));

    Parameter prm0(Lc, Lo);
    prm0.k[0][0] = -1.100; prm0.k[1][0] = -1.128;
    prm0.k[0][1] = -1.110; prm0.k[1][1] = -1.104;
    prm0.z[0][0] = +10.50; prm0.z[1][0] = +10.70;
    prm0.z[0][1] = +12.40; prm0.z[1][1] = +10.50;
    prm0.eta[0].setPoint(0.3000, 0.6000);
    prm0.eta[1].setPoint(0.6000, 0.2000);
    prm0.xi[0].setPoint(0.5000, 0.8000);
    prm0.xi[1].setPoint(0.2500, 0.3000);
    jfunc.setParameter0(prm0);

    Parameter prm(Lc, Lo);
    prm.k[0][0] = -1.86; prm.k[1][0] = -1.98;
    prm.k[0][1] = -1.93; prm.k[1][1] = -1.82;
    prm.z[0][0] = +10.5; prm.z[1][0] = +11.7;
    prm.z[0][1] = +11.4; prm.z[1][1] = +10.5;
    prm.eta[0].setPoint(0.3164,0.6854);
    prm.eta[1].setPoint(0.5847,0.3524);
    prm.xi[0].setPoint(0.7341,0.8248);
    prm.xi[1].setPoint(0.2116,0.2329);
    jfunc.setParameter(prm);

    ConjugateGradient g;
    g.setFunction(&jfunc);
    g.setGradient(&jfunc);
    g.setPrinter(&jfunc);
    g.setProjection(&jfunc);
    g.setEpsilon1(0.000);
    g.setEpsilon2(0.000);
    g.setEpsilon3(0.000);
    g.setR1MinimizeEpsilon(1.0, 0.001);
    g.setNormalize(true);
    g.showEndMessage(true);
    //g.setResetIteration(false);

    DoubleVector x; prm.toVector(x);
    g.calculate(x);
}

void ExpOptimalLetters::Table3Y2()
{
    JFunctional jfunc;
    unsigned int Lc = 2;
    unsigned int Lo = 2;

    jfunc.setGridParameters(Dimension(0.005, 0, 200), Dimension(0.010, 0, 100), Dimension(0.010, 0, 100));
    jfunc.U.resize(101, 101, 10.0);

    DoubleVector fis; fis << +0.2;// << +0.3 << +0.5;
    DoubleVector p_fis(fis.length(), 1.0/fis.length());

    DoubleVector thetas; thetas << +6.3;// << +6.4 << +6.5;
    DoubleVector p_thetas(thetas.length(), 1.0/thetas.length());

    jfunc.setInitTemperatures(fis, p_fis);
    jfunc.setEnvrTemperatures(thetas, p_thetas);

    jfunc.setEquationParameters(1.0, 0.01, 0.01);
    jfunc.setEpsilon(0.0);
    jfunc.setPenaltyCoefficient(20.0);
    jfunc.setPenaltyLimits(DoubleVector(Lc, -5.0), DoubleVector(Lc, +20.0));

    Parameter prm0(Lc, Lo);
    prm0.k[0][0] = -1.100; prm0.k[1][0] = -1.128;
    prm0.k[0][1] = -1.110; prm0.k[1][1] = -1.104;
    prm0.z[0][0] = +10.50; prm0.z[1][0] = +10.70;
    prm0.z[0][1] = +12.40; prm0.z[1][1] = +10.50;
    prm0.eta[0].setPoint(0.3000, 0.6000);
    prm0.eta[1].setPoint(0.6000, 0.2000);
    prm0.xi[0].setPoint(0.5000, 0.8000);
    prm0.xi[1].setPoint(0.2500, 0.3000);
    jfunc.setParameter0(prm0);

    Parameter prm(Lc, Lo);
    prm.k[0][0] = -1.30; prm.k[1][0] = -1.15;
    prm.k[0][1] = -1.65; prm.k[1][1] = -1.34;
    prm.z[0][0] = +12.54; prm.z[1][0] = +12.58;
    prm.z[0][1] = +14.44; prm.z[1][1] = +13.48;
    prm.eta[0].setPoint(0.7574,0.6384);
    prm.eta[1].setPoint(0.3685,0.1275);
    prm.xi[0].setPoint(0.2514, 0.5254);
    prm.xi[1].setPoint(0.6147, 0.2487);
    jfunc.setParameter(prm);

    ConjugateGradient g;
    g.setFunction(&jfunc);
    g.setGradient(&jfunc);
    g.setPrinter(&jfunc);
    g.setProjection(&jfunc);
    g.setEpsilon1(0.000);
    g.setEpsilon2(0.000);
    g.setEpsilon3(0.000);
    g.setR1MinimizeEpsilon(1.0, 0.001);
    g.setNormalize(true);
    g.showEndMessage(true);
    //g.setResetIteration(false);

    DoubleVector x; prm.toVector(x);
    g.calculate(x);
}

void ExpOptimalLetters::test()
{
    JFunctional jfunc;
    unsigned int Lc = 2;
    unsigned int Lo = 2;

    jfunc.setGridParameters(Dimension(0.005, 0, 200), Dimension(0.010, 0, 100), Dimension(0.010, 0, 100));
    jfunc.U.resize(101, 101, 10.0);

    DoubleVector fis; fis << +0.2;// << +0.3 << +0.5;
    DoubleVector p_fis(fis.length(), 1.0/fis.length());

    DoubleVector thetas; thetas << +6.3;// << +6.4 << +6.5;
    DoubleVector p_thetas(thetas.length(), 1.0/thetas.length());

    jfunc.setInitTemperatures(fis, p_fis);
    jfunc.setEnvrTemperatures(thetas, p_thetas);

    jfunc.setEquationParameters(1.0, 0.01, 0.01);
    jfunc.setEpsilon(0.0);
    jfunc.setPenaltyCoefficient(20.0);
    jfunc.setPenaltyLimits(DoubleVector(Lc, -5.0), DoubleVector(Lc, +20.0));

    Parameter prm0(Lc, Lo);
//    prm0.k[0][0] = -1.100; prm0.k[1][0] = -1.128;
//    prm0.k[0][1] = -1.110; prm0.k[1][1] = -1.104;
//    prm0.z[0][0] = +10.50; prm0.z[1][0] = +10.70;
//    prm0.z[0][1] = +12.40; prm0.z[1][1] = +10.50;
//    prm0.eta[0].setPoint(0.3000, 0.6000);
//    prm0.eta[1].setPoint(0.6000, 0.2000);
//    prm0.xi[0].setPoint(0.5000, 0.8000);
//    prm0.xi[1].setPoint(0.2500, 0.3000);
    jfunc.setParameter0(prm0);

//    Parameter prm(Lc, Lo);
//    prm.k[0][0] = -1.8600;  prm.k[1][0] = -1.9800;
//    prm.k[0][1] = -1.9300;  prm.k[1][1] = -1.8200;
//    prm.z[0][0] = +10.5000; prm.z[1][0] = +11.7000;
//    prm.z[0][1] = +11.4000; prm.z[1][1] = +10.5000;
//    prm.eta[0].setPoint(0.3164,0.6854);
//    prm.eta[1].setPoint(0.5847,0.3524);
//    prm.xi[0].setPoint(0.7341, 0.8248);
//    prm.xi[1].setPoint(0.2116, 0.2329);
//    jfunc.setParameter(prm);

    Parameter prm(Lc, Lo);
    prm.k[0][0] = -0.8000;  prm.k[1][0] = -0.7751;
    prm.k[0][1] = -0.7661;  prm.k[1][1] = -0.7906;
    prm.z[0][0] = +10.2793; prm.z[1][0] = +11.4658;
    prm.z[0][1] = +11.1711; prm.z[1][1] = +10.2847;
    prm.eta[0].setPoint(0.6829,0.1908);
    prm.eta[1].setPoint(0.3947,0.6322);
    prm.xi[0].setPoint(0.6471, 0.3204);
    prm.xi[1].setPoint(0.5307, 0.5904);
    jfunc.setParameter(prm);

//    DoubleMatrix u;
//    vector<ExtendedSpaceNode2D> info;
//    jfunc.setIntTemperature(fis[0]);
//    jfunc.setEnvTemperature(thetas[0]);
//    jfunc.forward->calculateMVD(u, info, false);
//    QPixmap pxm;
//    visualizeMatrixHeat(u, u.min(), u.max(), pxm, 101, 101);
//    pxm.save("imageU.png", "PNG");
//    DoubleVector v;
//    prm.toVector(v);
//    printf("%f %f %f\n", u.min(), u.max(), jfunc.fx(v));
//    return;

    IPrinter::printSeperatorLine();
    DoubleVector pv;
    prm.toVector(pv);
    DoubleVector ag(pv.length());
    IPrinter::print(pv, pv.length(), 6, 4);
    IPrinter::printSeperatorLine();

    puts("Calculating gradients....");
    jfunc.gradient(pv,ag);
    puts("Gradients are calculated.");

    double functional = jfunc.fx(pv);
    printf("Functional: %f\n", functional);

    DoubleVector ng(pv.length());
    ng.resize(pv.length(), 0.0);

    puts("Calculating numerical gradients....");
    //IGradient::Gradient(&jfunc, 0.010, pv, ng, 0*prm.Lc*prm.Lo,          1*prm.Lc*prm.Lo-1);
    //IGradient::Gradient(&jfunc, 0.010, pv, ng, 1*prm.Lc*prm.Lo,          2*prm.Lc*prm.Lo-1);
    //IGradient::Gradient(&jfunc, 0.010, pv, ng, 2*prm.Lc*prm.Lo+0*Lc,     2*prm.Lc*prm.Lo+2*prm.Lc-1);
    //IGradient::Gradient(&jfunc, 0.010, pv, ng, 2*prm.Lc*prm.Lo+2*prm.Lc, 2*prm.Lc*prm.Lo+2*prm.Lc+2*prm.Lo-1);
    puts("Numerical gradients are calculated.");

    //k------------------------------------------------------//
    DoubleVector pk = pv.mid(0, prm.Lc*prm.Lo-1);
    DoubleVector ak = ag.mid(0, prm.Lc*prm.Lo-1);
    DoubleVector nk = ng.mid(0, prm.Lc*prm.Lo-1);

    IPrinter::print(pk,pk.length(),14,4);
    IPrinter::print(ak,ak.length(),14,4); ak.L2Normalize();
    IPrinter::print(nk,nk.length(),14,4); nk.L2Normalize();
    IPrinter::print(ak,ak.length(),14,4);
    IPrinter::print(nk,nk.length(),14,4);
    IPrinter::printSeperatorLine();

    //z------------------------------------------------------//
    DoubleVector pz = pv.mid(prm.Lc*prm.Lo, 2*prm.Lc*prm.Lo-1);
    DoubleVector az = ag.mid(prm.Lc*prm.Lo, 2*prm.Lc*prm.Lo-1);
    DoubleVector nz = ng.mid(prm.Lc*prm.Lo, 2*prm.Lc*prm.Lo-1);

    IPrinter::print(pz,pz.length(),14,4);
    IPrinter::print(az,az.length(),14,4); az.L2Normalize();
    IPrinter::print(nz,nz.length(),14,4); nz.L2Normalize();
    IPrinter::print(az,az.length(),14,4);
    IPrinter::print(nz,nz.length(),14,4);
    IPrinter::printSeperatorLine();

    //eta------------------------------------------------------//
    DoubleVector pe = pv.mid(2*prm.Lc*prm.Lo, 2*prm.Lc*prm.Lo+2*prm.Lc-1);
    DoubleVector ae = ag.mid(2*prm.Lc*prm.Lo, 2*prm.Lc*prm.Lo+2*prm.Lc-1);
    DoubleVector ne = ng.mid(2*prm.Lc*prm.Lo, 2*prm.Lc*prm.Lo+2*prm.Lc-1);

    IPrinter::print(pe,pe.length(),14,4);
    IPrinter::print(ae,ae.length(),14,4); ae.L2Normalize();
    IPrinter::print(ne,ne.length(),14,4); ne.L2Normalize();
    IPrinter::print(ae,ae.length(),14,4);
    IPrinter::print(ne,ne.length(),14,4);
    IPrinter::printSeperatorLine();

    //xi------------------------------------------------------//
    DoubleVector px = pv.mid(2*prm.Lc*prm.Lo+2*prm.Lc, 2*prm.Lc*prm.Lo+2*prm.Lc+2*prm.Lo-1);
    DoubleVector ax = ag.mid(2*prm.Lc*prm.Lo+2*prm.Lc, 2*prm.Lc*prm.Lo+2*prm.Lc+2*prm.Lo-1);
    DoubleVector nx = ng.mid(2*prm.Lc*prm.Lo+2*prm.Lc, 2*prm.Lc*prm.Lo+2*prm.Lc+2*prm.Lo-1);

    IPrinter::print(px,px.length(),14,4);
    IPrinter::print(ax,ax.length(),14,4); ax.L2Normalize();
    IPrinter::print(nx,nx.length(),14,4); nx.L2Normalize();
    IPrinter::print(ax,ax.length(),14,4);
    IPrinter::print(nx,nx.length(),14,4);
    IPrinter::printSeperatorLine();
}

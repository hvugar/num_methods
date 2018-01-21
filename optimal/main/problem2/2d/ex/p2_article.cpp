#include "p2_article.h"

void Problem2Article::Main(int argc, char *argv[])
{
    //Table1_Y1();
    //Table2_Y1();
    //Table3_Y1();
    Table23_Y1();

    //Table1_Y2();
    //Table2_Y2();
    //Table3_Y2();

    //Table2Y1();
    //Table2Y2();
}

void Problem2Article::Table1_Y1()
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
    jfunc.setEpsilon(0.001);
    jfunc.setPenaltyCoefficient(500.0);
    jfunc.setPenaltyLimits(DoubleVector(Lc, -5.0), DoubleVector(Lc, +20.0));

    Parameter prm0(Lc, Lo);
    prm0.k[0][0] = -2.004; prm0.k[0][1] = +0.704;
    prm0.k[1][0] = +0.714; prm0.k[1][1] = -2.384;
    prm0.z[0][0] = +7.964; prm0.z[0][1] = +5.833;
    prm0.z[1][0] = +7.677; prm0.z[1][1] = +9.717;
    prm0.eta[0].setPoint(0.4149,0.7549);
    prm0.eta[1].setPoint(0.4052,0.7077);
    prm0.xi[0].setPoint(0.0501,0.0501);
    prm0.xi[1].setPoint(0.9500,0.0751);
    jfunc.setParameter0(prm0);

    DoubleVector hx; prm0.toVector(hx);
    IPrinter::print(hx, hx.length(), 6, 4);

    Parameter prm(Lc, Lo);
    prm.k[0][0] = -1.12; prm.k[0][1] = -1.24;
    prm.k[1][0] = -1.38; prm.k[1][1] = -1.58;
    prm.z[0][0] = +4.50; prm.z[0][1] = +3.40;
    prm.z[1][0] = +2.70; prm.z[1][1] = +3.50;
    prm.eta[0].setPoint(0.4574,0.8614);
    prm.eta[1].setPoint(0.2375,0.2347);
    prm.xi[0].setPoint(0.6911,0.5511);
    prm.xi[1].setPoint(0.8244,0.6700);
    jfunc.setParameter(prm);

    //DoubleMatrix u;
    //vector<ExtendedSpaceNode2D> info;
    //jfunc.setIntTemperature(fis[0]);
    //jfunc.setEnvTemperature(thetas[0]);
    //jfunc.forward->calculateMVD(u, info, false);
    //QPixmap pxm;
    //visualizeMatrixHeat(u, u.min(), u.max(), pxm, 101, 101);
    //pxm.save("imageU.png", "PNG");
    //DoubleVector v;
    //prm.toVector(v);
    //printf("%f %f %f\n", u.min(), u.max(), jfunc.fx(v));
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
    IGradient::Gradient(&jfunc, 0.001, pv, ng, 0*prm.Lc*prm.Lo,          1*prm.Lc*prm.Lo-1);
    IGradient::Gradient(&jfunc, 0.001, pv, ng, 1*prm.Lc*prm.Lo,          2*prm.Lc*prm.Lo-1);
    IGradient::Gradient(&jfunc, 0.001, pv, ng, 2*prm.Lc*prm.Lo+0*Lc,     2*prm.Lc*prm.Lo+2*prm.Lc-1);
    IGradient::Gradient(&jfunc, 0.001, pv, ng, 2*prm.Lc*prm.Lo+2*prm.Lc, 2*prm.Lc*prm.Lo+2*prm.Lc+2*prm.Lo-1);
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

void Problem2Article::Table2_Y1()
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
    jfunc.setEpsilon(0.001);
    jfunc.setPenaltyCoefficient(20.0);
    jfunc.setPenaltyLimits(DoubleVector(Lc, -5.0), DoubleVector(Lc, +20.0));

    Parameter prm0(Lc, Lo);
    prm0.k[0][0] = -2.004; prm0.k[0][1] = +0.704;
    prm0.k[1][0] = +0.714; prm0.k[1][1] = -2.384;
    prm0.z[0][0] = +7.964; prm0.z[0][1] = +5.833;
    prm0.z[1][0] = +7.677; prm0.z[1][1] = +9.717;
    prm0.eta[0].setPoint(0.4149,0.7549);
    prm0.eta[1].setPoint(0.4052,0.7077);
    prm0.xi[0].setPoint(0.0501,0.0501);
    prm0.xi[1].setPoint(0.9500,0.0751);
    jfunc.setParameter0(prm0);

    DoubleVector hx; prm0.toVector(hx);
    IPrinter::print(hx, hx.length(), 6, 4);

    Parameter prm(Lc, Lo);
    prm.k[0][0] = -1.12; prm.k[0][1] = -1.24;
    prm.k[1][0] = -1.38; prm.k[1][1] = -1.58;
    prm.z[0][0] = +4.50; prm.z[0][1] = +3.40;
    prm.z[1][0] = +2.70; prm.z[1][1] = +3.50;
    prm.eta[0].setPoint(0.4574,0.8614);
    prm.eta[1].setPoint(0.2375,0.2347);
    prm.xi[0].setPoint(0.6911,0.5511);
    prm.xi[1].setPoint(0.8244,0.6700);
    jfunc.setParameter(prm);

    ConjugateGradient g;
    g.setFunction(&jfunc);
    g.setGradient(&jfunc);
    g.setPrinter(&jfunc);
    g.setProjection(&jfunc);
    g.setEpsilon1(0.0000001);
    g.setEpsilon2(0.0000001);
    g.setEpsilon3(0.0000001);
    g.setR1MinimizeEpsilon(0.1, 0.001);
    g.setNormalize(true);
    g.showEndMessage(true);
    //g.setResetIteration(false);

    DoubleVector x; prm.toVector(x);
    g.calculate(x);
}

void Problem2Article::Table3_Y1()
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
    jfunc.setEpsilon(0.001);
    jfunc.setPenaltyCoefficient(500.0);
    jfunc.setPenaltyLimits(DoubleVector(Lc, -5.0), DoubleVector(Lc, +20.0));

    Parameter prm0(Lc, Lo);
    prm0.k[0][0] = -2.004; prm0.k[0][1] = +0.704;
    prm0.k[1][0] = +0.714; prm0.k[1][1] = -2.384;
    prm0.z[0][0] = +7.964; prm0.z[0][1] = +5.833;
    prm0.z[1][0] = +7.677; prm0.z[1][1] = +9.717;
    prm0.eta[0].setPoint(0.4149,0.7549);
    prm0.eta[1].setPoint(0.4052,0.7077);
    prm0.xi[0].setPoint(0.0501,0.0501);
    prm0.xi[1].setPoint(0.9500,0.0751);
    jfunc.setParameter0(prm0);

    DoubleVector hx; prm0.toVector(hx);
    IPrinter::print(hx, hx.length(), 6, 4);

    Parameter prm(Lc, Lo);
    prm.k[0][0] = -1.12; prm.k[0][1] = -1.24;
    prm.k[1][0] = -1.38; prm.k[1][1] = -1.58;
    prm.z[0][0] = +4.50; prm.z[0][1] = +3.40;
    prm.z[1][0] = +2.70; prm.z[1][1] = +3.50;
    prm.eta[0].setPoint(0.4574,0.8614);
    prm.eta[1].setPoint(0.2375,0.2347);
    prm.xi[0].setPoint(0.6911,0.5511);
    prm.xi[1].setPoint(0.8244,0.6700);
    jfunc.setParameter(prm);

    ConjugateGradient g;
    g.setFunction(&jfunc);
    g.setGradient(&jfunc);
    g.setPrinter(&jfunc);
    g.setProjection(&jfunc);
    g.setEpsilon1(0.0000001);
    g.setEpsilon2(0.0000001);
    g.setEpsilon3(0.0000001);
    g.setR1MinimizeEpsilon(0.1, 0.001);
    g.setNormalize(true);
    g.showEndMessage(true);
    //g.setResetIteration(false);

    DoubleVector x; prm.toVector(x);
    g.calculate(x);
}

void Problem2Article::Table23_Y1()
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
    jfunc.setEpsilon(1.00);
    jfunc.setPenaltyCoefficient(0.0);
    jfunc.setPenaltyLimits(DoubleVector(Lc, -5.0), DoubleVector(Lc, +20.0));

//    Parameter prm0(Lc, Lo);
//    prm0.k[0][0] = -2.336; prm0.k[0][1] = -2.191;
//    prm0.k[1][0] = -2.108; prm0.k[1][1] = -2.655;
//    prm0.z[0][0] = +10.00; prm0.z[0][1] = +10.00;
//    prm0.z[1][0] = +10.00; prm0.z[1][1] = +10.00;
//    prm0.eta[0].setPoint(0.4077,0.8433);
//    prm0.eta[1].setPoint(0.6046,0.0859);
//    prm0.xi[0].setPoint(0.2525,0.7500);
//    prm0.xi[1].setPoint(0.8882,0.5500);
//    jfunc.setParameter0(prm0);

    Parameter prm0(Lc, Lo);
    prm0.k[0][0] = -2.004; prm0.k[0][1] = +0.704;
    prm0.k[1][0] = +0.714; prm0.k[1][1] = -2.384;
    prm0.z[0][0] = +7.964; prm0.z[0][1] = +5.833;
    prm0.z[1][0] = +7.677; prm0.z[1][1] = +9.717;
    prm0.eta[0].setPoint(0.4149,0.7549);
    prm0.eta[1].setPoint(0.4052,0.7077);
    prm0.xi[0].setPoint(0.0501,0.0501);
    prm0.xi[1].setPoint(0.9500,0.0751);
    jfunc.setParameter0(prm0);

    DoubleVector hx; prm0.toVector(hx);
    IPrinter::print(hx, hx.length(), 6, 4);

//    Parameter prm(Lc, Lo);
//    prm.k[0][0] = -1.12; prm.k[0][1] = -1.64;
//    prm.k[1][0] = -1.38; prm.k[1][1] = -1.75;
//    prm.z[0][0] = +8.50; prm.z[0][1] = +9.40;
//    prm.z[1][0] = +9.70; prm.z[1][1] = +8.50;
//    prm.eta[0].setPoint(0.4574,0.8614);
//    prm.eta[1].setPoint(0.6375,0.0347);
//    prm.xi[0].setPoint(0.2911,0.7511);
//    prm.xi[1].setPoint(0.8544,0.5700);
//    jfunc.setParameter(prm);

    Parameter prm(Lc, Lo);
    prm.k[0][0] = -1.12; prm.k[0][1] = -1.24;
    prm.k[1][0] = -1.38; prm.k[1][1] = -1.58;
    prm.z[0][0] = +4.50; prm.z[0][1] = +3.40;
    prm.z[1][0] = +2.70; prm.z[1][1] = +3.50;
    prm.eta[0].setPoint(0.4574,0.8614);
    prm.eta[1].setPoint(0.2375,0.2347);
    prm.xi[0].setPoint(0.6911,0.5511);
    prm.xi[1].setPoint(0.8244,0.6700);
    jfunc.setParameter(prm);

    PFunctional pf;
    pf.jfunc = &jfunc;
    DoubleVector x; prm.toVector(x);
    pf.calculate(x);
    //printf("%f\n", jfunc.fx(x));
}

//////////////////////////////////////////////////////////////////////////////

void Problem2Article::Table1_Y2()
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
    jfunc.setEpsilon(0.001);
    jfunc.setPenaltyCoefficient(500.0);
    jfunc.setPenaltyLimits(DoubleVector(Lc, -5.0), DoubleVector(Lc, +20.0));

    Parameter prm0(Lc, Lo);
    prm0.k[0][0] = -0.85;  prm0.k[1][0] = -0.78;
    prm0.k[0][1] = -0.77;  prm0.k[1][1] = -0.79;
    prm0.z[0][0] = +10.28; prm0.z[1][0] = +11.47;
    prm0.z[0][1] = +11.17; prm0.z[1][1] = +10.29;
    prm0.eta[0].setPoint(0.6829,0.1908);
    prm0.eta[1].setPoint(0.3947,0.6322);
    prm0.xi[0].setPoint(0.6471,0.3204);
    prm0.xi[1].setPoint(0.5307,0.5904);
    jfunc.setParameter0(prm0);

    DoubleVector hx; prm0.toVector(hx);
    IPrinter::print(hx, hx.length(), 6, 4);

    Parameter prm(Lc, Lo);
    prm.k[0][0] = -1.15; prm.k[1][0] = -1.25;
    prm.k[0][1] = -1.18; prm.k[1][1] = -1.11;
    prm.z[0][0] = +8.48; prm.z[1][0] = +9.70;
    prm.z[0][1] = +7.04; prm.z[1][1] = +10.45;
    prm.eta[0].setPoint(0.3164,0.6854);
    prm.eta[1].setPoint(0.5847,0.3524);
    prm.xi[0].setPoint(0.7341,0.8248);
    prm.xi[1].setPoint(0.2116,0.2329);
    jfunc.setParameter(prm);

    //DoubleMatrix u;
    //vector<ExtendedSpaceNode2D> info;
    //jfunc.setIntTemperature(fis[0]);
    //jfunc.setEnvTemperature(thetas[0]);
    //jfunc.forward->calculateMVD(u, info, false);
    //QPixmap pxm;
    //visualizeMatrixHeat(u, u.min(), u.max(), pxm, 101, 101);
    //pxm.save("imageU.png", "PNG");
    //DoubleVector v;
    //prm.toVector(v);
    //printf("%f %f %f\n", u.min(), u.max(), jfunc.fx(v));
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
    IGradient::Gradient(&jfunc, 0.001, pv, ng, 0*prm.Lc*prm.Lo,          1*prm.Lc*prm.Lo-1);
    IGradient::Gradient(&jfunc, 0.001, pv, ng, 1*prm.Lc*prm.Lo,          2*prm.Lc*prm.Lo-1);
    IGradient::Gradient(&jfunc, 0.001, pv, ng, 2*prm.Lc*prm.Lo+0*Lc,     2*prm.Lc*prm.Lo+2*prm.Lc-1);
    IGradient::Gradient(&jfunc, 0.001, pv, ng, 2*prm.Lc*prm.Lo+2*prm.Lc, 2*prm.Lc*prm.Lo+2*prm.Lc+2*prm.Lo-1);
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

void Problem2Article::Table2_Y2()
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
    jfunc.setEpsilon(0.001);
    jfunc.setPenaltyCoefficient(20.0);
    jfunc.setPenaltyLimits(DoubleVector(Lc, -5.0), DoubleVector(Lc, +20.0));

    Parameter prm0(Lc, Lo);
    prm0.k[0][0] = -0.85;  prm0.k[1][0] = -0.78;
    prm0.k[0][1] = -0.77;  prm0.k[1][1] = -0.79;
    prm0.z[0][0] = +10.28; prm0.z[1][0] = +11.47;
    prm0.z[0][1] = +11.17; prm0.z[1][1] = +10.29;
    prm0.eta[0].setPoint(0.6829,0.1908);
    prm0.eta[1].setPoint(0.3947,0.6322);
    prm0.xi[0].setPoint(0.6471,0.3204);
    prm0.xi[1].setPoint(0.5307,0.5904);
    jfunc.setParameter0(prm0);

    DoubleVector hx; prm0.toVector(hx);
    IPrinter::print(hx, hx.length(), 6, 4);

    Parameter prm(Lc, Lo);
    prm.k[0][0] = -1.15; prm.k[1][0] = -1.25;
    prm.k[0][1] = -1.18; prm.k[1][1] = -1.11;
    prm.z[0][0] = +8.48; prm.z[1][0] = +9.70;
    prm.z[0][1] = +7.04; prm.z[1][1] = +10.45;
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
    g.setEpsilon1(0.0000001);
    g.setEpsilon2(0.0000001);
    g.setEpsilon3(0.0000001);
    g.setR1MinimizeEpsilon(0.1, 0.001);
    g.setNormalize(true);
    g.showEndMessage(true);
    //g.setResetIteration(false);

    DoubleVector x; prm.toVector(x);
    g.calculate(x);
}

void Problem2Article::Table3_Y2()
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
    jfunc.setEpsilon(0.001);
    jfunc.setPenaltyCoefficient(500.0);
    jfunc.setPenaltyLimits(DoubleVector(Lc, -5.0), DoubleVector(Lc, +20.0));

    Parameter prm0(Lc, Lo);
    prm0.k[0][0] = -0.85;  prm0.k[1][0] = -0.78;
    prm0.k[0][1] = -0.77;  prm0.k[1][1] = -0.79;
    prm0.z[0][0] = +10.28; prm0.z[1][0] = +11.47;
    prm0.z[0][1] = +11.17; prm0.z[1][1] = +10.29;
    prm0.eta[0].setPoint(0.6829,0.1908);
    prm0.eta[1].setPoint(0.3947,0.6322);
    prm0.xi[0].setPoint(0.6471,0.3204);
    prm0.xi[1].setPoint(0.5307,0.5904);
    jfunc.setParameter0(prm0);

    DoubleVector hx; prm0.toVector(hx);
    IPrinter::print(hx, hx.length(), 6, 4);

    Parameter prm(Lc, Lo);
    prm.k[0][0] = -1.15; prm.k[1][0] = -1.25;
    prm.k[0][1] = -1.18; prm.k[1][1] = -1.11;
    prm.z[0][0] = +8.48; prm.z[1][0] = +9.70;
    prm.z[0][1] = +7.04; prm.z[1][1] = +10.45;
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
    g.setEpsilon1(0.0000001);
    g.setEpsilon2(0.0000001);
    g.setEpsilon3(0.0000001);
    g.setR1MinimizeEpsilon(0.1, 0.001);
    g.setNormalize(true);
    g.showEndMessage(true);
    //g.setResetIteration(false);

    DoubleVector x; prm.toVector(x);
    g.calculate(x);
}

/////////////////////////////////////////////////////////////////////////////////

void Problem2Article::Table1Y1()
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

void Problem2Article::Table1Y2()
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

void Problem2Article::Table2Y1()
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
    jfunc.setEpsilon(0.001);
    jfunc.setPenaltyCoefficient(500.0);
    jfunc.setPenaltyLimits(DoubleVector(Lc, -5.0), DoubleVector(Lc, +20.0));

//    Parameter prm0(Lc, Lo);
//    prm0.k[0][0] = -1.100; prm0.k[1][0] = -1.128;
//    prm0.k[0][1] = -1.110; prm0.k[1][1] = -1.104;
//    prm0.z[0][0] = +10.50; prm0.z[1][0] = +10.70;
//    prm0.z[0][1] = +12.40; prm0.z[1][1] = +10.50;
//    prm0.eta[0].setPoint(0.3000, 0.6000);
//    prm0.eta[1].setPoint(0.6000, 0.2000);
//    prm0.xi[0].setPoint(0.5000, 0.8000);
//    prm0.xi[1].setPoint(0.2500, 0.3000);
//    jfunc.setParameter0(prm0);

    Parameter prm0(Lc, Lo);
    prm0.k[0][0] = -0.85;  prm0.k[1][0] = -0.78;
    prm0.k[0][1] = -0.77;  prm0.k[1][1] = -0.79;
    prm0.z[0][0] = +10.28; prm0.z[1][0] = +11.47;
    prm0.z[0][1] = +11.17; prm0.z[1][1] = +10.29;
    prm0.eta[0].setPoint(0.6829,0.1908);
    prm0.eta[1].setPoint(0.3947,0.6322);
    prm0.xi[0].setPoint(0.6471,0.3204);
    prm0.xi[1].setPoint(0.5307,0.5904);
    jfunc.setParameter0(prm0);

    DoubleVector hx; prm0.toVector(hx);
    IPrinter::print(hx, hx.length(), 6, 4);

    Parameter prm(Lc, Lo);
    prm.k[0][0] = -1.15; prm.k[1][0] = -1.25;
    prm.k[0][1] = -1.18; prm.k[1][1] = -1.11;
    prm.z[0][0] = +8.48; prm.z[1][0] = +9.70;
    prm.z[0][1] = +7.04; prm.z[1][1] = +10.45;
    prm.eta[0].setPoint(0.3164,0.6854);
    prm.eta[1].setPoint(0.5847,0.3524);
    prm.xi[0].setPoint(0.7341,0.8248);
    prm.xi[1].setPoint(0.2116,0.2329);
    jfunc.setParameter(prm);

//    Parameter prm(Lc, Lo);
//    prm.k[0][0] = -1.12; prm.k[1][0] = -1.38;
//    prm.k[0][1] = -1.24; prm.k[1][1] = -1.58;
//    prm.z[0][0] = +4.5; prm.z[1][0] = +2.7;
//    prm.z[0][1] = +3.4; prm.z[1][1] = +3.5;
//    prm.eta[0].setPoint(0.4574,0.8614);
//    prm.eta[1].setPoint(0.2375,0.2347);
//    prm.xi[0].setPoint(0.6911,0.5511);
//    prm.xi[1].setPoint(0.8244,0.6700);
//    jfunc.setParameter(prm);

    ConjugateGradient g;
    g.setFunction(&jfunc);
    g.setGradient(&jfunc);
    g.setPrinter(&jfunc);
    g.setProjection(&jfunc);
    g.setEpsilon1(0.0000001);
    g.setEpsilon2(0.0000001);
    g.setEpsilon3(0.0000001);
    g.setR1MinimizeEpsilon(0.1, 0.001);
    g.setNormalize(true);
    g.showEndMessage(true);
    //g.setResetIteration(false);

    DoubleVector x; prm.toVector(x);
    g.calculate(x);
}





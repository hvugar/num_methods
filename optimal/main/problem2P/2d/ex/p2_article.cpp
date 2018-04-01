#include "p2_article.h"

void Problem2Article::Main(int argc, char *argv[])
{
    //Table1_Y1();
    //Table2_Y1();
    //Table3_Y1();
    //Table23_Y1();

    //Table1_Y2();
    //Table2_Y2();
    //Table3_Y2();

    //Table2Y1();
    //Table2Y2();

    image1();
}

void Problem2Article::Table1_Y1()
{
    JFunctional jfunc;
    jfunc.optimizeK = jfunc.optimizeZ = jfunc.optimizeC = jfunc.optimizeO = true;
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
    jfunc.setRegEpsilon(0.0);
    jfunc.setPenaltyCoefficient(0.1);
    jfunc.setPenaltyLimits(DoubleVector(Lc, -5.0), DoubleVector(Lc, +20.0));

    Parameter prm0(Lc, Lo);
    prm0.k[0][0] = -2.00; prm0.k[0][1] = +0.70;
    prm0.k[1][0] = +0.71; prm0.k[1][1] = -2.38;
    prm0.z[0][0] = +7.96; prm0.z[0][1] = +5.83;
    prm0.z[1][0] = +7.68; prm0.z[1][1] = +9.72;
    prm0.eta[0].setPoint(0.4149,0.7549);
    prm0.eta[1].setPoint(0.4052,0.7077);
    prm0.xi[0].setPoint(0.0501,0.0501);
    prm0.xi[1].setPoint(0.9500,0.0751);
    jfunc.setParameter0(prm0);

    DoubleVector hx; jfunc.toVector(prm0, hx);
    IPrinter::print(hx, hx.length(), 6, 4);

//    Parameter prm(Lc, Lo);
//    prm.k[0][0] = -5.85; prm.k[0][1] = -3.48;
//    prm.k[1][0] = -4.74; prm.k[1][1] = -6.15;
//    prm.z[0][0] = +14.91; prm.z[0][1] = +11.45;
//    prm.z[1][0] = +16.84; prm.z[1][1] = +12.38;
//    prm.eta[0].setPoint(0.85,0.86);
//    prm.eta[1].setPoint(0.23,0.23);
//    prm.xi[0].setPoint(0.69,0.65);
//    prm.xi[1].setPoint(0.42,0.47);
//    jfunc.setParameter(prm);

    Parameter prm(Lc, Lo);
    prm.k[0][0] = -2.12; prm.k[0][1] = +1.24;
    prm.k[1][0] = -2.38; prm.k[1][1] = +2.58;
    prm.z[0][0] = +8.50; prm.z[0][1] = +7.40;
    prm.z[1][0] = +7.70; prm.z[1][1] = +9.50;
    prm.eta[0].setPoint(0.46,0.85);
    prm.eta[1].setPoint(0.24,0.24);
    prm.xi[0].setPoint(0.63,0.52);
    prm.xi[1].setPoint(0.84,0.68);
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
    jfunc.toVector(prm, pv);
    DoubleVector ag(pv.length());
    IPrinter::print(pv, pv.length(), 6, 4);
    IPrinter::printSeperatorLine();

    puts("Calculating gradients....");
    jfunc.gradient(pv,ag);
    puts("Gradients are calculated.");

    double functional = jfunc.fx(pv);
    printf("Functional: %f\n", functional);

    DoubleVector ng1(pv.length(), 0.0);
    DoubleVector ng2(pv.length(), 0.0);

    {
        puts("Calculating numerical gradients.... hx=0.01");
        IGradient::Gradient(&jfunc, 0.01, pv, ng1, 0*prm.Lc*prm.Lo,          1*prm.Lc*prm.Lo-1);
        IGradient::Gradient(&jfunc, 0.01, pv, ng1, 1*prm.Lc*prm.Lo,          2*prm.Lc*prm.Lo-1);
        IGradient::Gradient(&jfunc, 0.01, pv, ng1, 2*prm.Lc*prm.Lo+0*prm.Lc, 2*prm.Lc*prm.Lo+2*prm.Lc-1);
        IGradient::Gradient(&jfunc, 0.01, pv, ng1, 2*prm.Lc*prm.Lo+2*prm.Lc, 2*prm.Lc*prm.Lo+2*prm.Lc+2*prm.Lo-1);
        puts("Numerical gradients are calculated.");

        puts("Calculating numerical gradients.... hx=0.001");
        IGradient::Gradient(&jfunc, 0.001, pv, ng2, 0*prm.Lc*prm.Lo,          1*prm.Lc*prm.Lo-1);
        IGradient::Gradient(&jfunc, 0.001, pv, ng2, 1*prm.Lc*prm.Lo,          2*prm.Lc*prm.Lo-1);
        IGradient::Gradient(&jfunc, 0.001, pv, ng2, 2*prm.Lc*prm.Lo+0*prm.Lc, 2*prm.Lc*prm.Lo+2*prm.Lc-1);
        IGradient::Gradient(&jfunc, 0.001, pv, ng2, 2*prm.Lc*prm.Lo+2*prm.Lc, 2*prm.Lc*prm.Lo+2*prm.Lc+2*prm.Lo-1);
        puts("Numerical gradients are calculated.");

        //k------------------------------------------------------//
        DoubleVector pk0 = pv.mid(0, prm.Lc*prm.Lo-1);
        DoubleVector ak0 = ag.mid(0, prm.Lc*prm.Lo-1);
        DoubleVector nk1 = ng1.mid(0, prm.Lc*prm.Lo-1);
        DoubleVector nk2 = ng2.mid(0, prm.Lc*prm.Lo-1);

        IPrinter::print(pk0,pk0.length(),14,4);
        IPrinter::print(ak0,ak0.length(),14,4); ak0.L2Normalize();
        IPrinter::print(nk1,nk1.length(),14,4); nk1.L2Normalize();
        IPrinter::print(nk2,nk2.length(),14,4); nk2.L2Normalize();
        IPrinter::print(ak0,ak0.length(),14,4);
        IPrinter::print(nk1,nk1.length(),14,4);
        IPrinter::print(nk2,nk2.length(),14,4);
        IPrinter::printSeperatorLine();

        //z------------------------------------------------------//
        DoubleVector pz0 = pv.mid(prm.Lc*prm.Lo, 2*prm.Lc*prm.Lo-1);
        DoubleVector az0 = ag.mid(prm.Lc*prm.Lo, 2*prm.Lc*prm.Lo-1);
        DoubleVector nz1 = ng1.mid(prm.Lc*prm.Lo, 2*prm.Lc*prm.Lo-1);
        DoubleVector nz2 = ng2.mid(prm.Lc*prm.Lo, 2*prm.Lc*prm.Lo-1);

        IPrinter::print(pz0,pz0.length(),14,4);
        IPrinter::print(az0,az0.length(),14,4); az0.L2Normalize();
        IPrinter::print(nz1,nz1.length(),14,4); nz1.L2Normalize();
        IPrinter::print(nz2,nz2.length(),14,4); nz2.L2Normalize();
        IPrinter::print(az0,az0.length(),14,4);
        IPrinter::print(nz1,nz1.length(),14,4);
        IPrinter::print(nz2,nz2.length(),14,4);
        IPrinter::printSeperatorLine();

        //eta------------------------------------------------------//
        DoubleVector pe0 = pv.mid(2*prm.Lc*prm.Lo, 2*prm.Lc*prm.Lo+2*prm.Lc-1);
        DoubleVector ae0 = ag.mid(2*prm.Lc*prm.Lo, 2*prm.Lc*prm.Lo+2*prm.Lc-1);
        DoubleVector ne1 = ng1.mid(2*prm.Lc*prm.Lo, 2*prm.Lc*prm.Lo+2*prm.Lc-1);
        DoubleVector ne2 = ng2.mid(2*prm.Lc*prm.Lo, 2*prm.Lc*prm.Lo+2*prm.Lc-1);

        IPrinter::print(pe0,pe0.length(),14,4);
        IPrinter::print(ae0,ae0.length(),14,4); ae0.L2Normalize();
        IPrinter::print(ne1,ne1.length(),14,4); ne1.L2Normalize();
        IPrinter::print(ne2,ne2.length(),14,4); ne2.L2Normalize();
        IPrinter::print(ae0,ae0.length(),14,4);
        IPrinter::print(ne1,ne1.length(),14,4);
        IPrinter::print(ne2,ne2.length(),14,4);
        IPrinter::printSeperatorLine();

        //xi------------------------------------------------------//
        DoubleVector px0 = pv.mid(2*prm.Lc*prm.Lo+2*prm.Lc, 2*prm.Lc*prm.Lo+2*prm.Lc+2*prm.Lo-1);
        DoubleVector ax0 = ag.mid(2*prm.Lc*prm.Lo+2*prm.Lc, 2*prm.Lc*prm.Lo+2*prm.Lc+2*prm.Lo-1);
        DoubleVector nx1 = ng1.mid(2*prm.Lc*prm.Lo+2*prm.Lc, 2*prm.Lc*prm.Lo+2*prm.Lc+2*prm.Lo-1);
        DoubleVector nx2 = ng2.mid(2*prm.Lc*prm.Lo+2*prm.Lc, 2*prm.Lc*prm.Lo+2*prm.Lc+2*prm.Lo-1);

        IPrinter::print(px0,px0.length(),14,4);
        IPrinter::print(ax0,ax0.length(),14,4); ax0.L2Normalize();
        IPrinter::print(nx1,nx1.length(),14,4); nx1.L2Normalize();
        IPrinter::print(nx2,nx2.length(),14,4); nx2.L2Normalize();
        IPrinter::print(ax0,ax0.length(),14,4);
        IPrinter::print(nx1,nx1.length(),14,4);
        IPrinter::print(nx2,nx2.length(),14,4);
        IPrinter::printSeperatorLine();
    }
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
    jfunc.setRegEpsilon(0.001);
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

    DoubleVector hx; jfunc.toVector(prm0, hx);
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

    DoubleVector x; jfunc.toVector(prm, x);
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
    jfunc.setRegEpsilon(0.001);
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

    DoubleVector hx; jfunc.toVector(prm0, hx);
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

    DoubleVector x; jfunc.toVector(prm, x);
    g.calculate(x);
}

void Problem2Article::Table23_Y1()
{
    JFunctional jfunc;
    unsigned int Lc = 2;
    unsigned int Lo = 2;

    jfunc.setGridParameters(Dimension(0.005, 0, 200), Dimension(0.010, 0, 100), Dimension(0.010, 0, 100));
    jfunc.U.resize(101, 101, 10.0);

    DoubleVector fis; fis << +0.0;// << +0.3 << +0.5;
    DoubleVector p_fis(fis.length(), 1.0/fis.length());

    DoubleVector thetas; thetas << +6.3;// << +6.4 << +6.5;
    DoubleVector p_thetas(thetas.length(), 1.0/thetas.length());

    jfunc.setInitTemperatures(fis, p_fis);
    jfunc.setEnvrTemperatures(thetas, p_thetas);

    jfunc.setEquationParameters(1.0, 0.01, 0.01);
    jfunc.setRegEpsilon(0.00);
    jfunc.setPenaltyCoefficient(0.0);
    jfunc.setPenaltyLimits(DoubleVector(Lc, -5.0), DoubleVector(Lc, +20.0));

    Parameter prm0(Lc, Lo);
    prm0.k[0][0] = -2.00; prm0.k[0][1] = +0.70;
    prm0.k[1][0] = +0.71; prm0.k[1][1] = -2.38;
    prm0.z[0][0] = +7.96; prm0.z[0][1] = +5.83;
    prm0.z[1][0] = +7.68; prm0.z[1][1] = +9.72;
    prm0.eta[0].setPoint(0.42,0.76);
    prm0.eta[1].setPoint(0.41,0.71);
    prm0.xi[0].setPoint(0.25,0.35);
    prm0.xi[1].setPoint(0.95,0.88);
    jfunc.setParameter0(prm0);

    Parameter prm(Lc, Lo);
    prm.k[0][0] = -5.85; prm.k[0][1] = -3.48;
    prm.k[1][0] = -4.74; prm.k[1][1] = -6.15;
    prm.z[0][0] = +14.91; prm.z[0][1] = +11.45;
    prm.z[1][0] = +16.84; prm.z[1][1] = +12.38;
    prm.eta[0].setPoint(0.85,0.86);
    prm.eta[1].setPoint(0.23,0.23);
    prm.xi[0].setPoint(0.69,0.65);
    prm.xi[1].setPoint(0.42,0.47);
    jfunc.setParameter(prm);

//    Parameter prm(Lc, Lo);
//    prm.k[0][0] = -2.12; prm.k[0][1] = +1.24;
//    prm.k[1][0] = +2.38; prm.k[1][1] = +2.58;
//    prm.z[0][0] = +8.50; prm.z[0][1] = +7.40;
//    prm.z[1][0] = +7.70; prm.z[1][1] = +9.50;
//    prm.eta[0].setPoint(0.46,0.85);
//    prm.eta[1].setPoint(0.24,0.24);
//    prm.xi[0].setPoint(0.63,0.52);
//    prm.xi[1].setPoint(0.84,0.68);
//    jfunc.setParameter(prm);

    DoubleVector r; r << 0.100 << 0.200 << 0.400 << 0.800 << 1.000 << 2.000 << 5.000 << 10.00 << 20.00 << 50.00 << 100.0;
    DoubleVector e; e << 0.000 << 0.000 << 0.000 << 0.000 << 0.000 << 0.000 << 0.000 << 0.000 << 0.000 << 0.000 << 0.000;
    DoubleVector s; s << 10.00 << 5.000 << 2.000 << 1.000 << 1.000 << 1.000 << 1.000 << 1.000 << 1.000 << 1.000 << 1.000;

    DoubleVector x;
    jfunc.toVector(prm, x);

    for (unsigned int i=0; i<r.length(); i++)
    {
        jfunc.setPenaltyCoefficient(r[i]);
        jfunc.setRegEpsilon(e[i]);

        ConjugateGradient grad;
        //SteepestDescentGradient grad;
        grad.setFunction(&jfunc);
        grad.setPrinter(&jfunc);
        grad.setGradient(&jfunc);
        grad.setProjection(&jfunc);
        grad.setEpsilon1(0.0000001);
        grad.setEpsilon2(0.0000001);
        grad.setEpsilon3(0.0000001);
        grad.setR1MinimizeEpsilon(10.0, 0.00001);
        grad.setNormalize(true);
        grad.showEndMessage(false);
        grad.setResetIteration(false);
        grad.calculate(x);
    }
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
    jfunc.setRegEpsilon(0.001);
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

    DoubleVector hx; jfunc.toVector(prm0, hx);
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

    //    Parameter prm0(Lc, Lo);
    //    prm0.k[0][0] = -0.85;  prm0.k[1][0] = -0.78;
    //    prm0.k[0][1] = -0.77;  prm0.k[1][1] = -0.79;
    //    prm0.z[0][0] = +10.28; prm0.z[1][0] = +11.47;
    //    prm0.z[0][1] = +11.17; prm0.z[1][1] = +10.29;
    //    prm0.eta[0].setPoint(0.6829,0.1908);
    //    prm0.eta[1].setPoint(0.3947,0.6322);
    //    prm0.xi[0].setPoint(0.6471,0.3204);
    //    prm0.xi[1].setPoint(0.5307,0.5904);
    //    jfunc.setParameter0(prm0);

    //    DoubleVector hx; prm0.toVector(hx);
    //    IPrinter::print(hx, hx.length(), 6, 4);

    //    Parameter prm(Lc, Lo);
    //    prm.k[0][0] = -1.15; prm.k[1][0] = -1.25;
    //    prm.k[0][1] = -1.18; prm.k[1][1] = -1.11;
    //    prm.z[0][0] = +8.48; prm.z[1][0] = +9.70;
    //    prm.z[0][1] = +7.04; prm.z[1][1] = +10.45;
    //    prm.eta[0].setPoint(0.3164,0.6854);
    //    prm.eta[1].setPoint(0.5847,0.3524);
    //    prm.xi[0].setPoint(0.7341,0.8248);
    //    prm.xi[1].setPoint(0.2116,0.2329);
    //    jfunc.setParameter(prm);

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
    jfunc.toVector(prm, pv);
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
    jfunc.setRegEpsilon(0.001);
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

    DoubleVector hx; jfunc.toVector(prm0, hx);
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

    DoubleVector x; jfunc.toVector(prm, x);
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
    jfunc.setRegEpsilon(0.001);
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

    DoubleVector hx; jfunc.toVector(prm0, hx);
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

    DoubleVector x; jfunc.toVector(prm, x);
    g.calculate(x);
}

/////////////////////////////////////////////////////////////////////////////////

void Problem2Article::image1()
{
    JFunctional jfunc;
    jfunc.optimizeK = jfunc.optimizeZ = jfunc.optimizeC = jfunc.optimizeO = true;
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
    jfunc.setRegEpsilon(0.000);
    jfunc.setPenaltyCoefficient(0.1);
    jfunc.setPenaltyLimits(DoubleVector(Lc, -5.0), DoubleVector(Lc, +20.0));

    Parameter prm0(Lc, Lo);
    prm0.k[0][0] = -2.00; prm0.k[0][1] = +0.70;
    prm0.k[1][0] = +0.71; prm0.k[1][1] = -2.38;
    prm0.z[0][0] = +7.96; prm0.z[0][1] = +5.83;
    prm0.z[1][0] = +7.68; prm0.z[1][1] = +9.72;
    prm0.eta[0].setPoint(0.4149,0.7549);
    prm0.eta[1].setPoint(0.4052,0.7077);
    prm0.xi[0].setPoint(0.0501,0.0501);
    prm0.xi[1].setPoint(0.9500,0.0751);
    jfunc.setParameter0(prm0);

    DoubleVector hx; jfunc.toVector(prm0, hx);
    IPrinter::print(hx, hx.length(), 6, 4);

//    Parameter prm(Lc, Lo);
//    prm.k[0][0] = -5.85; prm.k[0][1] = -3.48;
//    prm.k[1][0] = -4.74; prm.k[1][1] = -6.15;
//    prm.z[0][0] = +14.91; prm.z[0][1] = +11.45;
//    prm.z[1][0] = +16.84; prm.z[1][1] = +12.38;
//    prm.eta[0].setPoint(0.85,0.86);
//    prm.eta[1].setPoint(0.23,0.23);
//    prm.xi[0].setPoint(0.69,0.65);
//    prm.xi[1].setPoint(0.42,0.47);
//    jfunc.setParameter(prm);

    Parameter prm(Lc, Lo);
    prm.k[0][0] = -2.0846; prm.k[0][1] = +0.2125;
    prm.k[1][0] = +2.4595; prm.k[1][1] = -3.5826;
    prm.z[0][0] = +10.1600; prm.z[0][1] = +9.2323;
    prm.z[1][0] = +15.0211; prm.z[1][1] = +13.680;
    prm.eta[0].setPoint(0.8827, 0.8579);
    prm.eta[1].setPoint(0.3493, 0.3284);
    prm.xi[0].setPoint(0.7935, 0.9400);
    prm.xi[1].setPoint(0.5650, 0.6197);
    jfunc.setParameter(prm);

    DoubleVector pv; jfunc.toVector(prm, pv);
    for (unsigned int i=2; i<=300; i+=2)
    {
        jfunc.setGridParameters(Dimension(0.005, 0, i), Dimension(0.010, 0, 100), Dimension(0.010, 0, 100));
        printf("%f %f\n", i*0.005, jfunc.fx(pv));
    }
}

#include "problem22dex5.h"

void Problem22DEx5::Main(int argc UNUSED_PARAM, char *argv[] UNUSED_PARAM)
{
    Table1Y1();
    //Table1Y2();
    //experiment2();
}

Problem22DEx5::Problem22DEx5() : AbstactProblem22D()
{   
    Dimension time = Dimension(0.005, 0, 200);
    Dimension dimX = Dimension(0.010, 0, 100);
    Dimension dimY = Dimension(0.010, 0, 100);
    setGridParameters(time, dimX, dimY);

    U.resize(dimY.sizeN()+1, dimX.sizeN()+1, 10.0);
}

void Problem22DEx5::Table1Y1()
{
    Problem22DEx5 p22d5;
    unsigned int Lc = 2;
    unsigned int Lo = 2;

    p22d5.fis.resize(3);
    p22d5.fis[0] = 0.2;
    p22d5.fis[1] = 0.3;
    p22d5.fis[2] = 0.5;

    p22d5.thetas.resize(3);
    p22d5.thetas[0] = 6.3;
    p22d5.thetas[1] = 6.4;
    p22d5.thetas[2] = 6.5;

    p22d5.setEquationParameters(1.0, 0.01, 0.01);
    //p22d5.setIntTemperature(0.1);
    //p22d5.setEnvTemperature(6.3);
    p22d5.setEpsilon(0.0);
    p22d5.setPenaltyCoefficient(1.0);
    p22d5.setPenaltyLimits(DoubleVector(Lc, -5.0), DoubleVector(Lc, +20.0));

    {
        Parameter rpm0(Lc, Lo);

        rpm0.k[0][0] = -1.100; rpm0.k[1][0] = -1.128;
        rpm0.k[0][1] = -1.110; rpm0.k[1][1] = -1.104;

        rpm0.z[0][0] = +10.50; rpm0.z[1][0] = +10.70;
        rpm0.z[0][1] = +12.40; rpm0.z[1][1] = +10.50;

        rpm0.eta[0].setPoint(0.3000, 0.6000);
        rpm0.eta[1].setPoint(0.6000, 0.2000);

        rpm0.xi[0].setPoint(0.5000, 0.8000);
        rpm0.xi[1].setPoint(0.2500, 0.3000);

        p22d5.setParameter0(rpm0);
    }

    Parameter rpm(Lc, Lo);

    rpm.k[0][0] = -1.86; rpm.k[1][0] = -1.98;
    rpm.k[0][1] = -1.93; rpm.k[1][1] = -1.82;

    rpm.z[0][0] = +10.5; rpm.z[1][0] = +11.7;
    rpm.z[0][1] = +11.4; rpm.z[1][1] = +10.5;

    rpm.eta[0].setPoint(0.3164,0.6854);
    rpm.eta[1].setPoint(0.5847,0.3524);

    rpm.xi[0].setPoint(0.7341,0.8248);
    rpm.xi[1].setPoint(0.2116,0.2329);

    p22d5.setParameter(rpm);

    DoubleMatrix u;
    vector<ExtendedSpaceNode2D> info;
    p22d5.forward->calculateMVD(u, info, false);
    QPixmap pxm;
    visualizeMatrixHeat(u, u.min(), u.max(), pxm, 101, 101);
    pxm.save("imageU.png", "PNG");
    DoubleVector v;
    rpm.toVector(v);
    printf("%f %f %f\n", u.min(), u.max(), p22d5.fx(v));
    //return;

    IPrinter::printSeperatorLine();
    DoubleVector pv;
    rpm.toVector(pv);
    DoubleVector ag(pv.length());
    IPrinter::print(pv, pv.length(), 6, 4);
    IPrinter::printSeperatorLine();

    puts("Calculating gradients....");
    p22d5.gradient(pv,ag);
    puts("Gradients are calculated.");

    double functional = p22d5.fx(pv);
    printf("Functional: %f\n", functional);

    DoubleVector ng(pv.length());
    ng.resize(pv.length(), 0.0);

    puts("Calculating numerical gradients....");
    IGradient::Gradient(&p22d5, 0.010, pv, ng, 0*rpm.Lc*rpm.Lo,          1*rpm.Lc*rpm.Lo-1);
    IGradient::Gradient(&p22d5, 0.010, pv, ng, 1*rpm.Lc*rpm.Lo,          2*rpm.Lc*rpm.Lo-1);
    IGradient::Gradient(&p22d5, 0.010, pv, ng, 2*rpm.Lc*rpm.Lo+0*Lc,     2*rpm.Lc*rpm.Lo+2*rpm.Lc-1);
    IGradient::Gradient(&p22d5, 0.010, pv, ng, 2*rpm.Lc*rpm.Lo+2*rpm.Lc, 2*rpm.Lc*rpm.Lo+2*rpm.Lc+2*rpm.Lo-1);
    puts("Numerical gradients are calculated.");

    //k------------------------------------------------------//
    DoubleVector pk = pv.mid(0, rpm.Lc*rpm.Lo-1);
    DoubleVector ak = ag.mid(0, rpm.Lc*rpm.Lo-1);
    DoubleVector nk = ng.mid(0, rpm.Lc*rpm.Lo-1);

    IPrinter::print(pk,pk.length(),14,4);
    IPrinter::print(ak,ak.length(),14,4); ak.L2Normalize();
    IPrinter::print(nk,nk.length(),14,4); nk.L2Normalize();
    IPrinter::print(ak,ak.length(),14,4);
    IPrinter::print(nk,nk.length(),14,4);
    IPrinter::printSeperatorLine();

    //z------------------------------------------------------//
    DoubleVector pz = pv.mid(rpm.Lc*rpm.Lo, 2*rpm.Lc*rpm.Lo-1);
    DoubleVector az = ag.mid(rpm.Lc*rpm.Lo, 2*rpm.Lc*rpm.Lo-1);
    DoubleVector nz = ng.mid(rpm.Lc*rpm.Lo, 2*rpm.Lc*rpm.Lo-1);

    IPrinter::print(pz,pz.length(),14,4);
    IPrinter::print(az,az.length(),14,4); az.L2Normalize();
    IPrinter::print(nz,nz.length(),14,4); nz.L2Normalize();
    IPrinter::print(az,az.length(),14,4);
    IPrinter::print(nz,nz.length(),14,4);
    IPrinter::printSeperatorLine();

    //eta------------------------------------------------------//
    DoubleVector pe = pv.mid(2*rpm.Lc*rpm.Lo, 2*rpm.Lc*rpm.Lo+2*rpm.Lc-1);
    DoubleVector ae = ag.mid(2*rpm.Lc*rpm.Lo, 2*rpm.Lc*rpm.Lo+2*rpm.Lc-1);
    DoubleVector ne = ng.mid(2*rpm.Lc*rpm.Lo, 2*rpm.Lc*rpm.Lo+2*rpm.Lc-1);

    IPrinter::print(pe,pe.length(),14,4);
    IPrinter::print(ae,ae.length(),14,4); ae.L2Normalize();
    IPrinter::print(ne,ne.length(),14,4); ne.L2Normalize();
    IPrinter::print(ae,ae.length(),14,4);
    IPrinter::print(ne,ne.length(),14,4);
    IPrinter::printSeperatorLine();

    //xi------------------------------------------------------//
    DoubleVector px = pv.mid(2*rpm.Lc*rpm.Lo+2*rpm.Lc, 2*rpm.Lc*rpm.Lo+2*rpm.Lc+2*rpm.Lo-1);
    DoubleVector ax = ag.mid(2*rpm.Lc*rpm.Lo+2*rpm.Lc, 2*rpm.Lc*rpm.Lo+2*rpm.Lc+2*rpm.Lo-1);
    DoubleVector nx = ng.mid(2*rpm.Lc*rpm.Lo+2*rpm.Lc, 2*rpm.Lc*rpm.Lo+2*rpm.Lc+2*rpm.Lo-1);

    IPrinter::print(px,px.length(),14,4);
    IPrinter::print(ax,ax.length(),14,4); ax.L2Normalize();
    IPrinter::print(nx,nx.length(),14,4); nx.L2Normalize();
    IPrinter::print(ax,ax.length(),14,4);
    IPrinter::print(nx,nx.length(),14,4);
    IPrinter::printSeperatorLine();
}

void Problem22DEx5::Table1Y2()
{
    Problem22DEx5 p22d5;
    unsigned int Lc = 2;
    unsigned int Lo = 2;

    //    Dimension time = Dimension(0.005, 0, 200);
    //    Dimension dimX = Dimension(0.01, 0, 100);
    //    Dimension dimY = Dimension(0.01, 0, 100);
    //    p22d5.setGridParameters(time, dimX, dimY);

    //    p22d5.U.resize(dimY.sizeN()+1, dimX.sizeN()+1, 10.0);

    p22d5.setEquationParameters(1.0, 0.01, 0.01);
    p22d5.setIntTemperature(0.1);
    p22d5.setEnvTemperature(6.3);
    p22d5.setEpsilon(0.0);
    p22d5.setPenaltyCoefficient(0.0);
    p22d5.setPenaltyLimits(DoubleVector(Lc, -5.0), DoubleVector(Lc, +10.0));

    {
        Parameter rpm0(Lc, Lo);

        rpm0.k[0][0] = -0.100; rpm0.k[1][0] = -0.128;
        rpm0.k[0][1] = -0.110; rpm0.k[1][1] = -0.104;

        rpm0.z[0][0] = +10.50; rpm0.z[1][0] = +10.70;
        rpm0.z[0][1] = +12.40; rpm0.z[1][1] = +10.50;

        rpm0.eta[0].setPoint(0.3000, 0.6000);
        rpm0.eta[1].setPoint(0.6000, 0.2000);

        rpm0.xi[0].setPoint(0.5000, 0.8000);
        rpm0.xi[1].setPoint(0.2500, 0.3000);

        p22d5.setParameter0(rpm0);
    }

    Parameter rpm(Lc, Lo);

    rpm.k[0][0] = -1.130; rpm.k[1][0] = -1.150;
    rpm.k[0][1] = -1.650; rpm.k[1][1] = -1.340;

    rpm.z[0][0] = +12.54; rpm.z[1][0] = +12.58;
    rpm.z[0][1] = +14.44; rpm.z[1][1] = +13.48;

    rpm.eta[0].setPoint(0.7341,0.8248);
    rpm.eta[1].setPoint(0.2116,0.2329);

    rpm.xi[0].setPoint(0.3164,0.6854);
    rpm.xi[1].setPoint(0.5847,0.3524);

    p22d5.setParameter(rpm);

    DoubleMatrix u;
    vector<ExtendedSpaceNode2D> info;
    p22d5.forward->calculateMVD(u, info, false);
    QPixmap pxm;
    visualizeMatrixHeat(u, u.min(), u.max(), pxm, 101, 101);
    pxm.save("imageU.png", "PNG");
    DoubleVector v;
    rpm.toVector(v);
    printf("%f %f %f\n", u.min(), u.max(), p22d5.fx(v));
    //return;

    IPrinter::printSeperatorLine();
    DoubleVector pv;
    rpm.toVector(pv);
    DoubleVector ag(pv.length());
    IPrinter::print(pv, pv.length(), 6, 4);
    IPrinter::printSeperatorLine();

    puts("Calculating gradients....");
    p22d5.gradient(pv,ag);
    puts("Gradients are calculated.");

    double functional = p22d5.fx(pv);
    printf("Functional: %f\n", functional);

    DoubleVector ng(pv.length());
    ng.resize(pv.length(), 0.0);

    puts("Calculating numerical gradients....");
    IGradient::Gradient(&p22d5, 0.010, pv, ng, 0*rpm.Lc*rpm.Lo,          1*rpm.Lc*rpm.Lo-1);
    IGradient::Gradient(&p22d5, 0.010, pv, ng, 1*rpm.Lc*rpm.Lo,          2*rpm.Lc*rpm.Lo-1);
    IGradient::Gradient(&p22d5, 0.010, pv, ng, 2*rpm.Lc*rpm.Lo+0*Lc,     2*rpm.Lc*rpm.Lo+2*rpm.Lc-1);
    IGradient::Gradient(&p22d5, 0.010, pv, ng, 2*rpm.Lc*rpm.Lo+2*rpm.Lc, 2*rpm.Lc*rpm.Lo+2*rpm.Lc+2*rpm.Lo-1);
    puts("Numerical gradients are calculated.");

    //k------------------------------------------------------//
    DoubleVector pk = pv.mid(0, rpm.Lc*rpm.Lo-1);
    DoubleVector ak = ag.mid(0, rpm.Lc*rpm.Lo-1);
    DoubleVector nk = ng.mid(0, rpm.Lc*rpm.Lo-1);

    IPrinter::print(pk,pk.length(),14,4);
    IPrinter::print(ak,ak.length(),14,4); ak.L2Normalize();
    IPrinter::print(nk,nk.length(),14,4); nk.L2Normalize();
    IPrinter::print(ak,ak.length(),14,4);
    IPrinter::print(nk,nk.length(),14,4);
    IPrinter::printSeperatorLine();

    //z------------------------------------------------------//
    DoubleVector pz = pv.mid(rpm.Lc*rpm.Lo, 2*rpm.Lc*rpm.Lo-1);
    DoubleVector az = ag.mid(rpm.Lc*rpm.Lo, 2*rpm.Lc*rpm.Lo-1);
    DoubleVector nz = ng.mid(rpm.Lc*rpm.Lo, 2*rpm.Lc*rpm.Lo-1);

    IPrinter::print(pz,pz.length(),14,4);
    IPrinter::print(az,az.length(),14,4); az.L2Normalize();
    IPrinter::print(nz,nz.length(),14,4); nz.L2Normalize();
    IPrinter::print(az,az.length(),14,4);
    IPrinter::print(nz,nz.length(),14,4);
    IPrinter::printSeperatorLine();

    //eta------------------------------------------------------//
    DoubleVector pe = pv.mid(2*rpm.Lc*rpm.Lo, 2*rpm.Lc*rpm.Lo+2*rpm.Lc-1);
    DoubleVector ae = ag.mid(2*rpm.Lc*rpm.Lo, 2*rpm.Lc*rpm.Lo+2*rpm.Lc-1);
    DoubleVector ne = ng.mid(2*rpm.Lc*rpm.Lo, 2*rpm.Lc*rpm.Lo+2*rpm.Lc-1);

    IPrinter::print(pe,pe.length(),14,4);
    IPrinter::print(ae,ae.length(),14,4); ae.L2Normalize();
    IPrinter::print(ne,ne.length(),14,4); ne.L2Normalize();
    IPrinter::print(ae,ae.length(),14,4);
    IPrinter::print(ne,ne.length(),14,4);
    IPrinter::printSeperatorLine();

    //xi------------------------------------------------------//
    DoubleVector px = pv.mid(2*rpm.Lc*rpm.Lo+2*rpm.Lc, 2*rpm.Lc*rpm.Lo+2*rpm.Lc+2*rpm.Lo-1);
    DoubleVector ax = ag.mid(2*rpm.Lc*rpm.Lo+2*rpm.Lc, 2*rpm.Lc*rpm.Lo+2*rpm.Lc+2*rpm.Lo-1);
    DoubleVector nx = ng.mid(2*rpm.Lc*rpm.Lo+2*rpm.Lc, 2*rpm.Lc*rpm.Lo+2*rpm.Lc+2*rpm.Lo-1);

    IPrinter::print(px,px.length(),14,4);
    IPrinter::print(ax,ax.length(),14,4); ax.L2Normalize();
    IPrinter::print(nx,nx.length(),14,4); nx.L2Normalize();
    IPrinter::print(ax,ax.length(),14,4);
    IPrinter::print(nx,nx.length(),14,4);
    IPrinter::printSeperatorLine();
}

void Problem22DEx5::experiment1()
{
    Problem22DEx5 p22d5;
    unsigned int Lc = 1;
    unsigned int Lo = 2;

    p22d5.setEquationParameters(1.0, 0.01, 0.01);
    p22d5.setIntTemperature(0.2);
    p22d5.setEnvTemperature(6.3);
    p22d5.setEpsilon(0.0);
    p22d5.setPenaltyCoefficient(0.0);
    p22d5.setPenaltyLimits(DoubleVector(Lc, -10.0), DoubleVector(Lc, +10.0));

    {
        Parameter rpm0(Lc, Lo);

        rpm0.k[0][0] = -0.100; //rpm0.k[1][0] = -0.128;
        rpm0.k[0][1] = -0.110; //rpm0.k[1][1] = -0.104;
        //rpm0.k[0][2] = -12.0; //rpm0.k[1][2] = -10.7;

        rpm0.z[0][0] = +10.5; //rpm0.z[1][0] = +10.7;
        rpm0.z[0][1] = +12.4; //rpm0.z[1][1] = +10.5;
        //rpm0.z[0][2] = +13.4; //rpm0.z[1][1] = +12.5;

        rpm0.eta[0].setPoint(0.3000, 0.6000);
        //rpm0.eta[1].setPoint(0.6000, 0.2000);

        rpm0.xi[0].setPoint(0.5000, 0.8000);
        rpm0.xi[1].setPoint(0.2500, 0.3000);
        //rpm0.xi[2].setPoint(0.8000, 0.4000);

        p22d5.setParameter0(rpm0);
        //p22d5.calculateU();
    }

    Parameter rpm(Lc, Lo);

    rpm.k[0][0] = -0.86; //rpm.k[1][0] = -0.98;
    rpm.k[0][1] = -0.83; //rpm.k[1][1] = -0.82;
    //rpm.k[0][2] = -11.2; rpm.k[1][2] = -10.7;

    rpm.z[0][0] = +10.5; //rpm.z[1][0] = +11.7;
    rpm.z[0][1] = +11.4; //rpm.z[1][1] = +12.5;
    //rpm.z[0][2] = +10.5; rpm.z[1][2] = +12.5;

    //rpm.eta[0].setPoint(0.3964, 0.5654);
    rpm.eta[0].setPoint(0.4000, 0.5000);
    //rpm.eta[0].setPoint(0.3964, 0.8010);
    //rpm.eta[1].setPoint(0.5847, 0.2124);

    rpm.xi[0].setPoint(0.7341, 0.8248);
    //rpm.xi[0].setPoint(0.7341, 0.8010);
    rpm.xi[1].setPoint(0.1516, 0.3329);
    //rpm.xi[2].setPoint(0.8157, 0.4359);

    p22d5.setParameter(rpm);

    DoubleMatrix u;
    vector<ExtendedSpaceNode2D> info;
    p22d5.forward->calculateMVD(u, info, false);
    QPixmap pxm;
    visualizeMatrixHeat(u, u.min(), u.max(), pxm, 101, 101);
    pxm.save("imageU.png", "PNG");
    DoubleVector v;
    rpm.toVector(v);
    printf("%f %f %f\n", u.min(), u.max(), p22d5.fx(v));
    //    return;

    //    DoubleVector v;
    //    rpm.toVector(v);
    //    printf("%f %f %f\n", u.min(), u.max(), p22d5.fx(v));
    //    double aa = v.at(2*Lc*Lo);
    //    v.at(2*Lc*Lo) = aa+0.01;
    //    IPrinter::print(v, v.length(), 6, 4);
    //    printf("%f %f %f\n", u.min(), u.max(), p22d5.fx(v));
    //    v.at(2*Lc*Lo) = aa-0.01;
    //    IPrinter::print(v, v.length(), 6, 4);
    //    printf("%f %f %f\n", u.min(), u.max(), p22d5.fx(v));
    //    v.at(2*Lc*Lo) = aa;

    IPrinter::printSeperatorLine();
    DoubleVector pv;
    rpm.toVector(pv);
    DoubleVector ag(pv.length());
    IPrinter::print(pv, pv.length(), 6, 4);
    IPrinter::printSeperatorLine();

    puts("Calculating gradients....");
    p22d5.gradient(pv,ag);
    puts("Gradients are calculated.");

    double functional = p22d5.fx(pv);
    printf("Functional: %f\n", functional);

    DoubleVector ng(pv.length());
    ng.resize(pv.length(), 0.0);

    puts("Calculating numerical gradients....");
    IGradient::Gradient(&p22d5, 0.001, pv, ng, 0*rpm.Lc*rpm.Lo,          1*rpm.Lc*rpm.Lo-1);
    IGradient::Gradient(&p22d5, 0.001, pv, ng, 1*rpm.Lc*rpm.Lo,          2*rpm.Lc*rpm.Lo-1);
    IGradient::Gradient(&p22d5, 0.010, pv, ng, 2*rpm.Lc*rpm.Lo+0*Lc,     2*rpm.Lc*rpm.Lo+2*rpm.Lc-1);
    IGradient::Gradient(&p22d5, 0.010, pv, ng, 2*rpm.Lc*rpm.Lo+2*rpm.Lc, 2*rpm.Lc*rpm.Lo+2*rpm.Lc+2*rpm.Lo-1);
    puts("Numerical gradients are calculated.");

    //k------------------------------------------------------//
    DoubleVector pk = pv.mid(0, rpm.Lc*rpm.Lo-1);
    DoubleVector ak = ag.mid(0, rpm.Lc*rpm.Lo-1);
    DoubleVector nk = ng.mid(0, rpm.Lc*rpm.Lo-1);

    IPrinter::print(pk,pk.length(),14,6);
    IPrinter::print(ak,ak.length(),14,6); ak.L2Normalize();
    IPrinter::print(nk,nk.length(),14,6); nk.L2Normalize();
    IPrinter::print(ak,ak.length(),14,6);
    IPrinter::print(nk,nk.length(),14,6);
    IPrinter::printSeperatorLine();

    //z------------------------------------------------------//
    DoubleVector pz = pv.mid(rpm.Lc*rpm.Lo, 2*rpm.Lc*rpm.Lo-1);
    DoubleVector az = ag.mid(rpm.Lc*rpm.Lo, 2*rpm.Lc*rpm.Lo-1);
    DoubleVector nz = ng.mid(rpm.Lc*rpm.Lo, 2*rpm.Lc*rpm.Lo-1);

    IPrinter::print(pz,pz.length(),14,6);
    IPrinter::print(az,az.length(),14,6); az.L2Normalize();
    IPrinter::print(nz,nz.length(),14,6); nz.L2Normalize();
    IPrinter::print(az,az.length(),14,6);
    IPrinter::print(nz,nz.length(),14,6);
    IPrinter::printSeperatorLine();

    //eta------------------------------------------------------//
    DoubleVector pe = pv.mid(2*rpm.Lc*rpm.Lo, 2*rpm.Lc*rpm.Lo+2*rpm.Lc-1);
    DoubleVector ae = ag.mid(2*rpm.Lc*rpm.Lo, 2*rpm.Lc*rpm.Lo+2*rpm.Lc-1);
    DoubleVector ne = ng.mid(2*rpm.Lc*rpm.Lo, 2*rpm.Lc*rpm.Lo+2*rpm.Lc-1);

    IPrinter::print(pe,pe.length(),14,6);
    IPrinter::print(ae,ae.length(),14,6); ae.L2Normalize();
    IPrinter::print(ne,ne.length(),14,6); ne.L2Normalize();
    IPrinter::print(ae,ae.length(),14,6);
    IPrinter::print(ne,ne.length(),14,6);
    IPrinter::printSeperatorLine();

    //xi------------------------------------------------------//
    DoubleVector px = pv.mid(2*rpm.Lc*rpm.Lo+2*rpm.Lc, 2*rpm.Lc*rpm.Lo+2*rpm.Lc+2*rpm.Lo-1);
    DoubleVector ax = ag.mid(2*rpm.Lc*rpm.Lo+2*rpm.Lc, 2*rpm.Lc*rpm.Lo+2*rpm.Lc+2*rpm.Lo-1);
    DoubleVector nx = ng.mid(2*rpm.Lc*rpm.Lo+2*rpm.Lc, 2*rpm.Lc*rpm.Lo+2*rpm.Lc+2*rpm.Lo-1);

    IPrinter::print(px,px.length(),14,6);
    IPrinter::print(ax,ax.length(),14,6); ax.L2Normalize();
    IPrinter::print(nx,nx.length(),14,6); nx.L2Normalize();
    IPrinter::print(ax,ax.length(),14,6);
    IPrinter::print(nx,nx.length(),14,6);
    IPrinter::printSeperatorLine();
}

void Problem22DEx5::experiment2()
{
    Problem22DEx5 p22d5;
    unsigned int Lc = 2;
    unsigned int Lo = 2;

    p22d5.setEquationParameters(1.0, 0.01, 0.01);
    p22d5.setIntTemperature(0.1);
    p22d5.setEnvTemperature(6.3);
    p22d5.setEpsilon(0.0);
    p22d5.setPenaltyCoefficient(0.0);
    p22d5.setPenaltyLimits(DoubleVector(Lc, -10.0), DoubleVector(Lc, +10.0));

    {
        Parameter rpm0(Lc, Lo);

        rpm0.k[0][0] = -1.2;  rpm0.k[1][0] = -12.8;
        rpm0.k[0][1] = -11.0; rpm0.k[1][1] = -10.4;
//        rpm0.k[0][2] = -12.0; rpm0.k[1][2] = -10.7;

        rpm0.z[0][0] = +10.5; rpm0.z[1][0] = +10.7;
        rpm0.z[0][1] = +12.4; rpm0.z[1][1] = +10.5;
        //rpm0.z[0][2] = +13.4; //rpm0.z[1][1] = +12.5;

        rpm0.eta[0].setPoint(0.3000, 0.6000);
        rpm0.eta[1].setPoint(0.6000, 0.2000);

        rpm0.xi[0].setPoint(0.5000, 0.8000);
        rpm0.xi[1].setPoint(0.2500, 0.3000);
        //rpm0.xi[2].setPoint(0.8000, 0.4000);

        p22d5.setParameter0(rpm0);
        //p22d5.calculateU();
    }

    Parameter rpm(Lc, Lo);

    rpm.k[0][0] = -0.86; rpm.k[1][0] = -0.98;
    rpm.k[0][1] = -0.83; rpm.k[1][1] = -0.82;
    //rpm.k[0][2] = -11.2; rpm.k[1][2] = -10.7;

    rpm.z[0][0] = +10.5; rpm.z[1][0] = +11.7;
    rpm.z[0][1] = +11.4; rpm.z[1][1] = +10.5;
    //rpm.z[0][2] = +10.5; rpm.z[1][2] = +12.5;

    rpm.eta[0].setPoint(0.3164, 0.5854);
    rpm.eta[1].setPoint(0.5847, 0.2124);

    rpm.xi[0].setPoint(0.7341, 0.8248);
    rpm.xi[1].setPoint(0.1516, 0.3329);
    //rpm.xi[2].setPoint(0.8157, 0.4359);

    p22d5.setParameter(rpm);

    DoubleVector x;
    rpm.toVector(x);
    p22d5.optimization(x);
}

double Problem22DEx5::fx(const DoubleVector &prms) const
{
    Problem22DEx5* ap22d = const_cast<Problem22DEx5*>(this);

    double p_fi = 1.0/3.0;
    double p_th = 1.0/3.0;

    double sum = 0.0;
    for (unsigned int i=0; i<fis.size(); i++)
    {
        ap22d->setIntTemperature(fis.at(i));
        for (unsigned j=0; j<thetas.size(); j++)
        {
            ap22d->setEnvTemperature(thetas.at(j));
            //printf("Init temp: %f Env temp: %f\n", fis.at(i), thetas.at(j));
            sum += AbstactProblem22D::fx(prms)*p_fi*p_th;
        }
    }
    return sum;
}

void Problem22DEx5::gradient(const DoubleVector &prms, DoubleVector &g)
{
    double p_fi = 1.0/3.0;
    double p_th = 1.0/3.0;

    g.resize(prms.length(), 0.0);

    for (unsigned int i=0; i<fis.size(); i++)
    {
        setIntTemperature(fis[i]);
        for (unsigned j=0; j<thetas.size(); j++)
        {
            setEnvTemperature(thetas[j]);

            printf("Init temp: %f Env temp: %f\n", fis.at(i), thetas.at(j));

            DoubleVector g0(prms.length(), 0.0);
            AbstactProblem22D::gradient(prms, g0);
            for (unsigned int k=0; k<prms.length(); k++)
            {
                g[k] += g0[k]*p_fi*p_th;
            }
            g0.clear();
        }
    }
}

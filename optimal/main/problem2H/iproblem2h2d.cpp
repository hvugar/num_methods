#include "iproblem2h2d.h"

#include "iproblem2hforward2d.h"
#include "iproblem2hbackward2d.h"
#include "iproblem2h2d_ifunctional.h"

using namespace IProblem2H;

void initParameters(EquationParameter &e_prm, OptimizeParameter &o_prm, OptimizeParameter &o_prm0)
{
    e_prm.a = 1.0;
    e_prm.lambda = +0.01;

    e_prm.Ns = 2;
    e_prm.q.resize(e_prm.Ns);
    e_prm.theta.resize(e_prm.Ns);

    //    e_prm.q[0] = -0.2; e_prm.theta[0].x = 0.5000; e_prm.theta[0].y = 0.5000;
    e_prm.q[0] = -0.2; e_prm.theta[0].x = 0.200; e_prm.theta[0].y = 0.200;
    e_prm.q[1] = -0.5; e_prm.theta[1].x = 0.800; e_prm.theta[1].y = 0.800;

    {
        e_prm.Nc = 2;
        e_prm.No = 2;

//        o_prm.k.resize(e_prm.Nc, e_prm.No, 0.0);
//        o_prm.z.resize(e_prm.Nc, e_prm.No, 0.0);
//        o_prm.k[0][0] = -0.9186; o_prm.k[0][1] = -0.9051; o_prm.k[1][0] = -1.0169; o_prm.k[1][1] = -1.0221;
//        o_prm.z[0][0] = +0.3232; o_prm.z[0][1] = +0.2617; o_prm.z[1][0] = +0.2469; o_prm.z[1][1] = +0.0996;
//        o_prm.xi.resize(e_prm.No);
//        o_prm.xi[0].x = 0.2078; o_prm.xi[0].y = 0.2080; o_prm.xi[1].x = 0.8455; o_prm.xi[1].y = 0.8491;
//        o_prm.eta.resize(e_prm.Nc);
//        o_prm.eta[0].x = 0.7623; o_prm.eta[0].y = 0.2708; o_prm.eta[1].x = 0.2636; o_prm.eta[1].y = 0.7636;

//        o_prm.k.resize(e_prm.Nc, e_prm.No, 0.0);
//        o_prm.z.resize(e_prm.Nc, e_prm.No, 0.0);
//        o_prm.xi.resize(e_prm.No);
//        o_prm.eta.resize(e_prm.Nc);
//        o_prm.k[0][0] = -1.7000; o_prm.k[0][1] = -1.3000; o_prm.k[1][0] = -1.6000; o_prm.k[1][1] = -1.5000;
//        o_prm.z[0][0] = +0.0000; o_prm.z[0][1] = +0.0000; o_prm.z[1][0] = +0.0000; o_prm.z[1][1] = +0.0000;
//        o_prm.xi[0].x = +0.4000; o_prm.xi[0].y = +0.6000; o_prm.xi[1].x = 0.9000; o_prm.xi[1].y = 0.2000;
//        o_prm.eta[0].x = +0.2000; o_prm.eta[0].y = +0.4500; o_prm.eta[1].x = 0.6500; o_prm.eta[1].y = 0.7500;

    }

    {
                e_prm.No = 3;
                e_prm.Nc = 2;

                o_prm.k.resize(e_prm.Nc, e_prm.No, 0.0);
                o_prm.z.resize(e_prm.Nc, e_prm.No, 0.0);
                o_prm.xi.resize(e_prm.No);
                o_prm.eta.resize(e_prm.Nc);
                o_prm.k[0][0] = -0.1200; o_prm.k[0][1] = -0.2400; o_prm.k[0][2] = -2.2400; o_prm.k[1][0] = -0.4500; o_prm.k[1][1] = -0.1800; o_prm.k[1][2] = -2.1800;
                o_prm.z[0][0] = +0.5000; o_prm.z[0][1] = +0.4000; o_prm.z[0][2] = +0.4000; o_prm.z[1][0] = +0.7000; o_prm.z[1][1] = +0.5000; o_prm.z[1][2] = +0.5000;
                o_prm.xi[0].x = 0.5000; o_prm.xi[0].y = 0.6000; o_prm.xi[1].x = 0.7000; o_prm.xi[1].y = 0.2000; o_prm.xi[2].x = 0.5000; o_prm.xi[2].y = 0.5000;
                o_prm.eta[0].x = 0.2000; o_prm.eta[0].y = 0.7000; o_prm.eta[1].x = 0.8000; o_prm.eta[1].y = 0.3000;
    }

    {
        //    e_prm.No = 1;
        //    e_prm.Nc = 2;

        //    OptimizeParameter o_prm;

        //    o_prm.xi.resize(e_prm.No);
        //    o_prm.xi[0].x = 0.3000; o_prm.xi[0].y = 0.6000;

        //    o_prm.eta.resize(e_prm.Nc);
        //    o_prm.eta[0].x = 0.5000; o_prm.eta[0].y = 0.6000;
        //    o_prm.eta[1].x = 0.7000; o_prm.eta[1].y = 0.3000;

        //    o_prm.k.resize(e_prm.Nc, e_prm.No, 0.0);
        //    o_prm.z.resize(e_prm.Nc, e_prm.No, 0.0);

        //    o_prm.k[0][0] = -2.12; //o_prm.k[0][1] = -2.24; //o_prm.k[0][2] = -2.24;
        //    o_prm.k[1][0] = -2.45; //o_prm.k[1][1] = -2.18; //o_prm.k[1][2] = -2.18;
        //    o_prm.z[0][0] = +0.50; //o_prm.z[0][1] = +0.40; //o_prm.z[0][2] = +0.40;
        //    o_prm.z[1][0] = +0.70; //o_prm.z[1][1] = +0.50; //o_prm.z[1][2] = +0.50;
    }

    {
        //    o_prm0.xi.resize(e_prm.No);
        //    o_prm0.xi[0].x = 0.400; o_prm0.xi[0].y = 0.400;
        //    o_prm0.xi[1].x = 0.600; o_prm0.xi[1].y = 0.600;
        //    o_prm0.xi[2].x = 0.500; o_prm0.xi[2].y = 0.500;

        //    o_prm0.eta.resize(e_prm.Nc);
        //    o_prm0.eta[0].x = 0.300; o_prm0.eta[0].y = 0.400;
        //    o_prm0.eta[1].x = 0.700; o_prm0.eta[1].y = 0.300;

        //    o_prm0.k.resize(e_prm.Nc, e_prm.No, 0.0);
        //    o_prm0.z.resize(e_prm.Nc, e_prm.No, 0.0);

        //    o_prm0.k[0][0] = -0.12; o_prm0.k[0][1] = -0.24; o_prm0.k[0][2] = -0.24;
        //    o_prm0.k[1][0] = -0.45; o_prm0.k[1][1] = -0.18; o_prm0.k[1][2] = -0.18;
        //    o_prm0.z[0][0] = +5.50; o_prm0.z[0][1] = +4.40; o_prm0.z[0][2] = +4.40;
        //    o_prm0.z[1][0] = +4.70; o_prm0.z[1][1] = +5.50; o_prm0.z[1][2] = +5.50;
    }
    o_prm0 = o_prm;

    //    I[ 40]: 0.024768 0.300001 0.324769 R:0.00 e:0.000
    //    k: -0.9095  -0.8781  -0.9951  -1.0015 z:  0.3378   0.2736   0.2531   0.1034   o:0.1995 0.2099 0.8449 0.8503   c:0.7553 0.2647 0.2676 0.7584
    //    k:  0.0332   0.0427   0.0295   0.0409 z:  0.0167   0.0161   0.0137   0.0138   o:-0.1201 0.1606 -0.0308 0.0131   c:-0.1623 -0.1512 0.1166 -0.0317
    //    0.004568 1

//    OptimizeParameter o_prm0;
//    o_prm0.k.resize(e_prm.Nc, e_prm.No, 0.0);
//    o_prm0.z.resize(e_prm.Nc, e_prm.No, 0.0);

//    o_prm0.k[0][0] = -2.6388; o_prm0.k[0][1] = -2.4754; o_prm0.k[1][0] = -2.8992; o_prm0.k[1][1] = -2.7987;
//    o_prm0.z[0][0] = +0.0991; o_prm0.z[0][1] = +0.1211; o_prm0.z[1][0] = -0.0860; o_prm0.z[1][1] = +0.3928;

//    o_prm0.xi.resize(e_prm.No);
//    o_prm0.xi[0].x = 0.2402; o_prm0.xi[0].y = 0.1762; o_prm0.xi[1].x = 0.5736; o_prm0.xi[1].y = 0.3214;

//    o_prm0.eta.resize(e_prm.Nc);
//    o_prm0.eta[0].x = 0.2743; o_prm0.eta[0].y = 0.7485; o_prm0.eta[1].x = 0.7657; o_prm0.eta[1].y = 0.2613;
}

void IProblem2H2D::Main(int argc UNUSED_PARAM, char *argv[] UNUSED_PARAM)
{
    //    forward();
    //    forwardS();
    checkGradient();
    IPrinter::printSeperatorLine();
    optimization1();
    //    optimization2();
}

void IProblem2H2D::forward()
{
    IProblem2HForward2D frw;
    frw.addSpaceDimension(Dimension(0.01, 0, 100));
    frw.addSpaceDimension(Dimension(0.01, 0, 100));
    frw.setTimeDimension(Dimension(0.01, 0, 5000));

    EquationParameter e_prm;
    e_prm.a = 1.0;
    e_prm.lambda = 0.0;

    e_prm.Ns = 1;
    e_prm.q.resize(e_prm.Ns);
    e_prm.theta.resize(e_prm.Ns);

    e_prm.q[0] = -1.0; e_prm.theta[0].x = 0.5000; e_prm.theta[0].y = 0.5000;
    //e_prm.q[0] = +1.0; e_prm.theta[0].x = 0.2000; e_prm.theta[0].y = 0.2000;
    //e_prm.q[1] = +5.0; e_prm.theta[1].x = 0.8000; e_prm.theta[1].y = 0.8000;

    e_prm.Nc = 0;
    e_prm.No = 0;

    OptimizeParameter o_prm;
    o_prm.xi.resize(e_prm.No);
    //o_prm.xi[0].x = 0.400; o_prm.xi[0].y = 0.400;
    //o_prm.xi[1].x = 0.600; o_prm.xi[1].y = 0.600;

    o_prm.eta.resize(e_prm.Nc);
    //o_prm.eta[0].x = 0.300; o_prm.eta[0].y = 0.700;
    //o_prm.eta[1].x = 0.700; o_prm.eta[1].y = 0.300;

    o_prm.k.resize(e_prm.Nc, e_prm.No, 0.0);
    o_prm.z.resize(e_prm.Nc, e_prm.No, 0.0);

    //    o_prm.k[0][0] = -1.12; o_prm.k[0][1] = -1.24;
    //    o_prm.k[1][0] = -1.08; o_prm.k[1][1] = -2.18;
    //o_prm.z[0][0] = +0.00; o_prm.z[0][1] = +0.00;
    //o_prm.z[1][0] = +0.00; o_prm.z[1][1] = +0.00;

    //o_prm.k[0][0] = +0.00; o_prm.k[0][1] = +0.00;
    //o_prm.k[1][0] = +0.00; o_prm.k[1][1] = +0.00;

    frw.mEquParameter = e_prm;
    frw.mOptParameter = o_prm;

    DoubleMatrix u;
    DoubleMatrix ut;
    std::vector<ExtendedSpaceNode2DH> u_info;
    frw.calculateMVD(u, ut, u_info, false);

    //    IPrinter::printMatrix(u);

    //    IProblem2HBackward2D bkw;
    //    bkw.addSpaceDimension(Dimension(0.01, 0, 100));
    //    bkw.addSpaceDimension(Dimension(0.01, 0, 100));
    //    bkw.setTimeDimension(Dimension(0.005, 0, 200));
    //    bkw.mEquParameter = e_prm;
    //    bkw.mOptParameter = o_prm;

    //    IFunctional ifunc;
    //    ifunc.alpha0 = 1.0; ifunc.V0.resize(101, 101, 0.0);
    //    ifunc.alpha1 = 1.0; ifunc.V1.resize(101, 101, 0.0);

    //    bkw.UT = u;
    //    bkw.UTt = ut;
    //    bkw.ifunc = &ifunc;

    //    DoubleMatrix p;
    //    std::vector<ExtendedSpaceNode2DH> p_info;
    //    bkw.calculateMVD(p, p_info, false, u_info);
}

void IProblem2H2D::checkGradient()
{
    IFunctional ifunc;
    ifunc.optimizeK = true;
    ifunc.optimizeZ = true;
    ifunc.optimizeC = true;
    ifunc.optimizeO = true;

    double hx, hy; hx = hy = 0.01;
    unsigned Nx, Ny; Nx = Ny = 100;
    ifunc.mSpaceDimensionX = Dimension(hx, 0, Nx);
    ifunc.mSpaceDimensionY = Dimension(hy, 0, Ny);
    ifunc.mTimeDimension = Dimension(0.01, 0, 100);

    ifunc.alpha0 = 1.0; ifunc.V0.resize(Ny+1, Nx+1, 0.0);
    ifunc.alpha1 = 1.0; ifunc.V1.resize(Ny+1, Nx+1, 0.0);

    EquationParameter e_prm;
    OptimizeParameter o_prm;
    OptimizeParameter o_prm0;
    initParameters(e_prm, o_prm, o_prm0);

    ifunc.mEquParameter = e_prm;
    ifunc.mOptParameter0 = o_prm0;
    ifunc.regEpsilon = 0.0;

    ifunc.r = 0.0;
    ifunc.vmin.resize(e_prm.Nc, -5.0);
    ifunc.vmax.resize(e_prm.Nc, +5.0);

    {
        IPrinter::printSeperatorLine();
        DoubleVector pv;
        ifunc.toVector(o_prm, pv);
        IPrinter::print(pv, pv.length(), 6, 4);
        IPrinter::printSeperatorLine();
        DoubleVector pv0;
        ifunc.toVector(o_prm0, pv0);
        IPrinter::print(pv0, pv0.length(), 6, 4);
        IPrinter::printSeperatorLine();
        DoubleVector ag(pv.length());

        double functional = ifunc.fx(pv);
        printf("Functional: %f\n", functional);
        //    return;
        puts("Calculating gradients....");
        ifunc.gradient(pv, ag);
        puts("Gradients are calculated.");

        DoubleVector ng1(pv.length(), 0.0);
        DoubleVector ng2(pv.length(), 0.0);

        puts("Calculating numerical gradients.... dh=0.01");
        IGradient::Gradient(&ifunc, 0.01, pv, ng1, 0*e_prm.Nc*e_prm.No,            1*e_prm.Nc*e_prm.No-1);
        IGradient::Gradient(&ifunc, 0.01, pv, ng1, 1*e_prm.Nc*e_prm.No,            2*e_prm.Nc*e_prm.No-1);
        IGradient::Gradient(&ifunc, 0.01, pv, ng1, 2*e_prm.Nc*e_prm.No+0*e_prm.No, 2*e_prm.Nc*e_prm.No+2*e_prm.No-1);
        IGradient::Gradient(&ifunc, 0.01, pv, ng1, 2*e_prm.Nc*e_prm.No+2*e_prm.No, 2*e_prm.Nc*e_prm.No+2*e_prm.No+2*e_prm.Nc-1);
        puts("Numerical gradients are calculated.");

        //    puts("Calculating numerical gradients.... hx=0.001");
        //    IGradient::Gradient(&ifunc, 0.001, pv, ng2, 0*e_prm.Nc*e_prm.No,            1*e_prm.Nc*e_prm.No-1);
        //    IGradient::Gradient(&ifunc, 0.001, pv, ng2, 1*e_prm.Nc*e_prm.No,            2*e_prm.Nc*e_prm.No-1);
        //    IGradient::Gradient(&ifunc, 0.001, pv, ng2, 2*e_prm.Nc*e_prm.No+0*e_prm.No, 2*e_prm.Nc*e_prm.No+2*e_prm.No-1);
        //    IGradient::Gradient(&ifunc, 0.001, pv, ng2, 2*e_prm.Nc*e_prm.No+2*e_prm.No, 2*e_prm.Nc*e_prm.No+2*e_prm.No+2*e_prm.Nc-1);
        //    puts("Numerical gradients are calculated.");

        //k------------------------------------------------------//
        IPrinter::printSeperatorLine("k");
        DoubleVector pk0 = pv.mid(0, e_prm.Nc*e_prm.No-1);
        DoubleVector ak0 = ag.mid(0, e_prm.Nc*e_prm.No-1);
        DoubleVector nk1 = ng1.mid(0, e_prm.Nc*e_prm.No-1);
        DoubleVector nk2 = ng2.mid(0, e_prm.Nc*e_prm.No-1);

        IPrinter::print(pk0,pk0.length(),14,4);
        IPrinter::print(ak0,ak0.length(),14,4); ak0.L2Normalize();
        IPrinter::print(nk1,nk1.length(),14,4); nk1.L2Normalize();
        IPrinter::print(nk2,nk2.length(),14,4); nk2.L2Normalize();
        IPrinter::print(ak0,ak0.length(),14,4);
        IPrinter::print(nk1,nk1.length(),14,4);
        IPrinter::print(nk2,nk2.length(),14,4);

        //z------------------------------------------------------//
        IPrinter::printSeperatorLine("z");
        DoubleVector pz0 = pv.mid(e_prm.Nc*e_prm.No, 2*e_prm.Nc*e_prm.No-1);
        DoubleVector az0 = ag.mid(e_prm.Nc*e_prm.No, 2*e_prm.Nc*e_prm.No-1);
        DoubleVector nz1 = ng1.mid(e_prm.Nc*e_prm.No, 2*e_prm.Nc*e_prm.No-1);
        DoubleVector nz2 = ng2.mid(e_prm.Nc*e_prm.No, 2*e_prm.Nc*e_prm.No-1);

        IPrinter::print(pz0,pz0.length(),14,4);
        IPrinter::print(az0,az0.length(),14,4); az0.L2Normalize();
        IPrinter::print(nz1,nz1.length(),14,4); nz1.L2Normalize();
        IPrinter::print(nz2,nz2.length(),14,4); nz2.L2Normalize();
        IPrinter::print(az0,az0.length(),14,4);
        IPrinter::print(nz1,nz1.length(),14,4);
        IPrinter::print(nz2,nz2.length(),14,4);

        //xi------------------------------------------------------//
        IPrinter::printSeperatorLine("xi");
        DoubleVector pe0 = pv.mid(2*e_prm.Nc*e_prm.No, 2*e_prm.Nc*e_prm.No+2*e_prm.No-1);
        DoubleVector ae0 = ag.mid(2*e_prm.Nc*e_prm.No, 2*e_prm.Nc*e_prm.No+2*e_prm.No-1);
        DoubleVector ne1 = ng1.mid(2*e_prm.Nc*e_prm.No, 2*e_prm.Nc*e_prm.No+2*e_prm.No-1);
        DoubleVector ne2 = ng2.mid(2*e_prm.Nc*e_prm.No, 2*e_prm.Nc*e_prm.No+2*e_prm.No-1);

        IPrinter::print(pe0,pe0.length(),14,4);
        IPrinter::print(ae0,ae0.length(),14,4); ae0.L2Normalize();
        IPrinter::print(ne1,ne1.length(),14,4); ne1.L2Normalize();
        IPrinter::print(ne2,ne2.length(),14,4); ne2.L2Normalize();
        IPrinter::print(ae0,ae0.length(),14,4);
        IPrinter::print(ne1,ne1.length(),14,4);
        IPrinter::print(ne2,ne2.length(),14,4);

        //eta------------------------------------------------------//
        IPrinter::printSeperatorLine("eta");
        DoubleVector px0 = pv.mid(2*e_prm.Nc*e_prm.No+2*e_prm.No, 2*e_prm.Nc*e_prm.No+2*e_prm.No+2*e_prm.Nc-1);
        DoubleVector ax0 = ag.mid(2*e_prm.Nc*e_prm.No+2*e_prm.No, 2*e_prm.Nc*e_prm.No+2*e_prm.No+2*e_prm.Nc-1);
        DoubleVector nx1 = ng1.mid(2*e_prm.Nc*e_prm.No+2*e_prm.No, 2*e_prm.Nc*e_prm.No+2*e_prm.No+2*e_prm.Nc-1);
        DoubleVector nx2 = ng2.mid(2*e_prm.Nc*e_prm.No+2*e_prm.No, 2*e_prm.Nc*e_prm.No+2*e_prm.No+2*e_prm.Nc-1);

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

void IProblem2H2D::optimization1()
{
    IFunctional ifunc;
    ifunc.optimizeK = true;
    ifunc.optimizeZ = true;
    ifunc.optimizeO = true;
    ifunc.optimizeC = true;

    double hx, hy; hx = hy = 0.01;
    unsigned Nx, Ny; Nx = Ny = 100;
    ifunc.mSpaceDimensionX = Dimension(hx, 0, Nx);
    ifunc.mSpaceDimensionY = (Dimension(hy, 0, Ny));
    ifunc.mTimeDimension = Dimension(0.01, 0, 100);

    ifunc.alpha0 = 1.00; ifunc.V0.resize(Ny+1, Nx+1, 0.0);
    ifunc.alpha1 = 1.00; ifunc.V1.resize(Ny+1, Nx+1, 0.0);

    EquationParameter e_prm;
    OptimizeParameter o_prm;
    OptimizeParameter o_prm0;
    initParameters(e_prm, o_prm, o_prm0);

    ifunc.mEquParameter = e_prm;
    ifunc.mOptParameter0 = o_prm0 = o_prm;
    ifunc.regEpsilon = 0.0;

    ifunc.r = 0.0;
    ifunc.vmin.resize(e_prm.Nc, -5.0);
    ifunc.vmax.resize(e_prm.Nc, +5.0);

    ConjugateGradient g;
    //SteepestDescentGradient g;
    g.setFunction(&ifunc);
    g.setGradient(&ifunc);
    g.setPrinter(&ifunc);
    g.setProjection(&ifunc);
    g.setEpsilon1(0.000000001);
    g.setEpsilon2(0.000000001);
    g.setEpsilon3(0.000000001);
    g.setR1MinimizeEpsilon(0.1, 0.001);
    g.setNormalize(true);
    g.showEndMessage(true);
    //    g.setResetIteration(false);

    DoubleVector x; ifunc.toVector(o_prm, x);
    g.calculate(x);
}

void IProblem2H2D::optimization2()
{
    IFunctional ifunc;
    ifunc.optimizeK = true;
    ifunc.optimizeZ = true;
    ifunc.optimizeO = true;
    ifunc.optimizeC = true;

    ifunc.mSpaceDimensionX = Dimension(0.01, 0, 100);
    ifunc.mSpaceDimensionY = (Dimension(0.01, 0, 100));
    ifunc.mTimeDimension = Dimension(0.01, 0, 500);

    ifunc.alpha0 = 1.00; ifunc.V0.resize(101, 101, 0.0);
    ifunc.alpha1 = 1.00; ifunc.V1.resize(101, 101, 0.0);

    EquationParameter e_prm;
    e_prm.a = +1.0;
    e_prm.lambda = +0.01;

    e_prm.Ns = 3;
    e_prm.q.resize(e_prm.Ns);
    e_prm.theta.resize(e_prm.Ns);

    e_prm.q[0] = -0.02; e_prm.theta[0].x = 0.300; e_prm.theta[0].y = 0.700;
    e_prm.q[1] = -0.05; e_prm.theta[1].x = 0.450; e_prm.theta[1].y = 0.300;
    e_prm.q[2] = -0.03; e_prm.theta[1].x = 0.800; e_prm.theta[1].y = 0.850;

    e_prm.Nc = 2;
    e_prm.No = 2;

    OptimizeParameter o_prm;
    o_prm.k.resize(e_prm.Nc, e_prm.No, 0.0);
    o_prm.z.resize(e_prm.Nc, e_prm.No, 0.0);

    o_prm.k[0][0] = -1.7000; o_prm.k[0][1] = -1.3000; o_prm.k[1][0] = -1.6000; o_prm.k[1][1] = -1.5000;

    o_prm.z[0][0] = +0.0000; o_prm.z[0][1] = +0.0000; o_prm.z[1][0] = +0.0000; o_prm.z[1][1] = +0.0000;

    o_prm.xi.resize(e_prm.No);
    o_prm.xi[0].x = 0.4000; o_prm.xi[0].y = 0.6000; o_prm.xi[1].x = 0.9000; o_prm.xi[1].y = 0.2000;

    o_prm.eta.resize(e_prm.Nc);
    o_prm.eta[0].x = 0.2000; o_prm.eta[0].y = 0.4500; o_prm.eta[1].x = 0.6500; o_prm.eta[1].y = 0.7500;

    //    OptimizeParameter o_prm0;
    //    o_prm0.xi.resize(e_prm.No);
    //    o_prm0.xi[0].x = 0.400; o_prm0.xi[0].y = 0.400;
    //    o_prm0.xi[1].x = 0.600; o_prm0.xi[1].y = 0.600;

    //    o_prm0.eta.resize(e_prm.Nc);
    //    o_prm0.eta[0].x = 0.300; o_prm0.eta[0].y = 0.700;
    //    o_prm0.eta[1].x = 0.700; o_prm0.eta[1].y = 0.300;

    //    o_prm0.k.resize(e_prm.Nc, e_prm.No, 0.0);
    //    o_prm0.z.resize(e_prm.Nc, e_prm.No, 0.0);

    //    o_prm0.k[0][0] = -0.12; o_prm0.k[0][1] = -0.24;
    //    o_prm0.k[1][0] = -0.45; o_prm0.k[1][1] = -0.18;
    //    o_prm0.z[0][0] = +5.50; o_prm0.z[0][1] = +4.40;
    //    o_prm0.z[1][0] = +4.70; o_prm0.z[1][1] = +5.50;

    ifunc.mEquParameter = e_prm;
    ifunc.mOptParameter0 = o_prm;
    ifunc.regEpsilon = 0.0;

    ifunc.r = 0.0;
    ifunc.vmin.resize(e_prm.Nc, -5.0);
    ifunc.vmax.resize(e_prm.Nc, +5.0);

    ConjugateGradient g;
    //SteepestDescentGradient g;
    g.setFunction(&ifunc);
    g.setGradient(&ifunc);
    g.setPrinter(&ifunc);
    g.setProjection(&ifunc);
    g.setEpsilon1(0.000000001);
    g.setEpsilon2(0.000000001);
    g.setEpsilon3(0.000000001);
    g.setR1MinimizeEpsilon(10.0, 0.01);
    g.setNormalize(true);
    g.showEndMessage(true);
    //    g.setResetIteration(false);

    DoubleVector x;
    ifunc.toVector(o_prm, x);
    g.calculate(x);

    //    IPrinter::printSeperatorLine();
    //    DoubleVector pv;
    //    ifunc.toVector(o_prm, pv);
    //    IPrinter::print(pv, pv.length(), 6, 4);
    //    IPrinter::printSeperatorLine();
    //    DoubleVector pv0;
    //    ifunc.toVector(o_prm0, pv0);
    //    IPrinter::print(pv0, pv0.length(), 6, 4);
    //    IPrinter::printSeperatorLine();
    //    DoubleVector ag(pv.length());
}

void IProblem2H2D::forwardS()
{
    IProblem2HForward2D frw;
    frw.addSpaceDimension(Dimension(0.01, 0, 100));
    frw.addSpaceDimension(Dimension(0.01, 0, 100));
    frw.setTimeDimension(Dimension(0.01, 0, 100000));

    EquationParameter e_prm;
    e_prm.a = 1.0;
    e_prm.lambda = +2.0;
    frw.mEquParameter = e_prm;

    DoubleMatrix u;
    DoubleMatrix ut;
    frw.calculateMVD_N(u, ut);

    IPrinter::printMatrix(u);
}

void IProblem2H2D::distributeDelta(const SpacePoint &pt, std::vector<ExtendedSpacePointNode> &nodes, unsigned int id,
                                   const Dimension &xd, const Dimension &yd, unsigned int k)
{
    double hx = xd.step();
    double hy = yd.step();
    unsigned int Nx = xd.sizeN();
    unsigned int Ny = yd.sizeN();

    double sigmaX = hx;
    double sigmaY = hy;

    unsigned int rx = (unsigned int)(round(pt.x*Nx));
    unsigned int ry = (unsigned int)(round(pt.y*Ny));

    double sumX = 0.0;
    for (unsigned int n=rx-k; n<=rx+k; n++)
    {
        sumX += exp(-((n*hx-pt.x)*(n*hx-pt.x))/(2.0*sigmaX*sigmaX));
    }
    sumX *= hx;

    double sumY = 0.0;
    for (unsigned int m=ry-k; m<=ry+k; m++)
    {
        sumY += exp(-((m*hy-pt.y)*(m*hy-pt.y))/(2.0*sigmaY*sigmaY));
    }
    sumY *= hy;

    double sigma = (sumX*sumY) / (2.0*M_PI);
    double factor = 1.0/((2.0*M_PI)*sigma);

    for (unsigned int n=rx-k; n<=rx+k; n++)
    {
        for (unsigned int m=ry-k; m<=ry+k; m++)
        {
            ExtendedSpacePointNode node;
            node.i = n; node.x = n*hx;
            node.j = m; node.y = m*hy;
            node.pt = pt; node.id = id;
            node.w = factor*exp(-0.5*(((node.x-pt.x)*(node.x-pt.x))/(sigmaX*sigmaX)+((node.y-pt.y)*(node.y-pt.y))/(sigmaY*sigmaY)));
            nodes.push_back(node);
        }
    }
}

ExtendedSpaceNode2DH::ExtendedSpaceNode2DH()
{
    x = y = 0.0;
    i = j = 0;
    wi = NULL;
    layerNumber = 0;
}

ExtendedSpaceNode2DH::~ExtendedSpaceNode2DH()
{}

void ExtendedSpaceNode2DH::setSpaceNode(const SpacePoint &sn)
{
    this->x = sn.x;
    this->y = sn.y;
}

void ExtendedSpaceNode2DH::extendWeights(const Dimension &dimX, const Dimension &dimY, unsigned int layerNumber, unsigned int rows, unsigned int cols)
{
    if (rows==3 && cols==3)
    {
        this->rows = rows;
        this->cols = cols;
        this->layerNumber = layerNumber;
        this->wi = new WISpaceNodePDE*[rows];
        for (unsigned int rw=0; rw<rows; rw++)
        {
            wi[rw] = new WISpaceNodePDE[cols];
            for (unsigned int cl=0; cl<cols; cl++)
            {
                wi[rw][cl].u = new double[layerNumber];
            }
        }

        unsigned int Nx = dimX.sizeN();
        unsigned int Ny = dimY.sizeN();

        double hx = dimX.step();
        double hy = dimY.step();

        unsigned int rx = (unsigned int)(round(x*Nx));
        unsigned int ry = (unsigned int)(round(y*Ny));

        double dx = 0.0;
        double dy = 0.0;

        wi[0][0].i = rx - 1; wi[0][0].x = wi[0][0].i*hx; wi[0][0].j = ry - 1; wi[0][0].y = wi[0][0].j*hy; dx = fabs(wi[0][0].x-x); dy = fabs(wi[0][0].y-y); wi[0][0].w = 0.0;
        wi[0][1].i = rx + 0; wi[0][1].x = wi[0][1].i*hx; wi[0][1].j = ry - 1; wi[0][1].y = wi[0][1].j*hy; dx = fabs(wi[0][1].x-x); dy = fabs(wi[0][1].y-y); wi[0][1].w = 0.0;
        wi[0][2].i = rx + 1; wi[0][2].x = wi[0][2].i*hx; wi[0][2].j = ry - 1; wi[0][2].y = wi[0][2].j*hy; dx = fabs(wi[0][2].x-x); dy = fabs(wi[0][2].y-y); wi[0][2].w = 0.0;
        wi[1][0].i = rx - 1; wi[1][0].x = wi[1][0].i*hx; wi[1][0].j = ry + 0; wi[1][0].y = wi[1][0].j*hy; dx = fabs(wi[1][0].x-x); dy = fabs(wi[1][0].y-y); wi[1][0].w = 0.0;
        wi[1][1].i = rx + 0; wi[1][1].x = wi[1][1].i*hx; wi[1][1].j = ry + 0; wi[1][1].y = wi[1][1].j*hy; dx = fabs(wi[1][1].x-x); dy = fabs(wi[1][1].y-y); wi[1][1].w = 0.0;
        wi[1][2].i = rx + 1; wi[1][2].x = wi[1][2].i*hx; wi[1][2].j = ry + 0; wi[1][2].y = wi[1][2].j*hy; dx = fabs(wi[1][2].x-x); dy = fabs(wi[1][2].y-y); wi[1][2].w = 0.0;
        wi[2][0].i = rx - 1; wi[2][0].x = wi[2][0].i*hx; wi[2][0].j = ry + 1; wi[2][0].y = wi[2][0].j*hy; dx = fabs(wi[2][0].x-x); dy = fabs(wi[2][0].y-y); wi[2][0].w = 0.0;
        wi[2][1].i = rx + 0; wi[2][1].x = wi[2][1].i*hx; wi[2][1].j = ry + 1; wi[2][1].y = wi[2][1].j*hy; dx = fabs(wi[2][1].x-x); dy = fabs(wi[2][1].y-y); wi[2][1].w = 0.0;
        wi[2][2].i = rx + 1; wi[2][2].x = wi[2][2].i*hx; wi[2][2].j = ry + 1; wi[2][2].y = wi[2][2].j*hy; dx = fabs(wi[2][2].x-x); dy = fabs(wi[2][2].y-y); wi[2][2].w = 0.0;
    }

    if (rows==4 && cols==4)
    {
        this->rows = rows;
        this->cols = cols;
        this->layerNumber = layerNumber;
        this->wi = new WISpaceNodePDE*[rows];
        for (unsigned int rw=0; rw<rows; rw++)
        {
            wi[rw] = new WISpaceNodePDE[cols];
            for (unsigned int cl=0; cl<cols; cl++)
            {
                wi[rw][cl].u = new double[layerNumber];
            }
        }

        unsigned int Nx = dimX.sizeN();
        unsigned int Ny = dimY.sizeN();

        double hx = dimX.step();
        double hy = dimY.step();

        unsigned int rx = (unsigned int)(floor(x*Nx));
        unsigned int ry = (unsigned int)(floor(y*Ny));

        double hx3 = hx*hx*hx;
        double hx32 = (1.0/(2.0*hx3));
        double hx36 = (1.0/(6.0*hx3));

        double hy3 = hy*hy*hy;
        double hy32 = (1.0/(2.0*hy3));
        double hy36 = (1.0/(6.0*hy3));

        double dx = 0.0;
        double dy = 0.0;

        wi[1][1].i = rx + 0; wi[1][1].x = wi[1][1].i*hx;
        wi[1][1].j = ry + 0; wi[1][1].y = wi[1][1].j*hy;
        dx = fabs(wi[1][1].x-x);
        dy = fabs(wi[1][1].y-y);
        wi[1][1].w = ((2.0*hx-dx)*(hx-dx)*(hx+dx)*hx32) * ((2.0*hy-dy)*(hy-dy)*(hy+dy)*hy32);

        wi[2][1].i = rx + 0; wi[2][1].x = wi[2][1].i*hx;
        wi[2][1].j = ry + 1; wi[2][1].y = wi[2][1].j*hy;
        dx = fabs(wi[2][1].x-x);
        dy = fabs(wi[2][1].y-y);
        wi[2][1].w = ((2.0*hx-dx)*(hx-dx)*(hx+dx)*hx32) * ((2.0*hy-dy)*(hy-dy)*(hy+dy)*hy32);

        wi[2][2].i = rx + 1; wi[2][2].x = wi[2][2].i*hx;
        wi[2][2].j = ry + 1; wi[2][2].y = wi[2][2].j*hy;
        dx = fabs(wi[2][2].x-x);
        dy = fabs(wi[2][2].y-y);
        wi[2][2].w = ((2.0*hx-dx)*(hx-dx)*(hx+dx)*hx32) * ((2.0*hy-dy)*(hy-dy)*(hy+dy)*hy32);

        wi[1][2].i = rx + 1; wi[1][2].x = wi[1][2].i*hx;
        wi[1][2].j = ry + 0; wi[1][2].y = wi[1][2].j*hy;
        dx = fabs(wi[1][2].x-x);
        dy = fabs(wi[1][2].y-y);
        wi[1][2].w = ((2.0*hx-dx)*(hx-dx)*(hx+dx)*hx32) * ((2.0*hy-dy)*(hy-dy)*(hy+dy)*hy32);

        wi[0][0].i = rx - 1; wi[0][0].x = wi[0][0].i*hx;
        wi[0][0].j = ry - 1; wi[0][0].y = wi[0][0].j*hy;
        dx = fabs(wi[0][0].x-x);
        dy = fabs(wi[0][0].y-y);
        wi[0][0].w = ((2.0*hx-dx)*(hx-dx)*(3.0*hx-dx)*hx36) * ((2.0*hy-dy)*(hy-dy)*(3.0*hy-dy)*hy36);

        wi[1][0].i = rx - 1; wi[1][0].x = wi[1][0].i*hx;
        wi[1][0].j = ry + 0; wi[1][0].y = wi[1][0].j*hy;
        dx = fabs(wi[1][0].x-x);
        dy = fabs(wi[1][0].y-y);
        wi[1][0].w = ((2.0*hx-dx)*(hx-dx)*(3.0*hx-dx)*hx36) * ((2.0*hy-dy)*(hy-dy)*(hy+dy)*hy32);

        wi[2][0].i = rx - 1; wi[2][0].x = wi[2][0].i*hx;
        wi[2][0].j = ry + 1; wi[2][0].y = wi[2][0].j*hy;
        dx = fabs(wi[2][0].x-x);
        dy = fabs(wi[2][0].y-y);
        wi[2][0].w = ((2.0*hx-dx)*(hx-dx)*(3.0*hx-dx)*hx36)*((2.0*hy-dy)*(hy-dy)*(hy+dy)*hy32);

        wi[3][0].i = rx - 1; wi[3][0].x = wi[3][0].i*hx;
        wi[3][0].j = ry + 2; wi[3][0].y = wi[3][0].j*hy;
        dx = fabs(wi[3][0].x-x);
        dy = fabs(wi[3][0].y-y);
        wi[3][0].w = ((2.0*hx-dx)*(hx-dx)*(3.0*hx-dx)*hx36)*((2.0*hy-dy)*(hy-dy)*(3.0*hy-dy)*hy36);

        wi[3][1].i = rx + 0; wi[3][1].x = wi[3][1].i*hx;
        wi[3][1].j = ry + 2; wi[3][1].y = wi[3][1].j*hy;
        dx = fabs(wi[3][1].x-x);
        dy = fabs(wi[3][1].y-y);
        wi[3][1].w = ((2.0*hx-dx)*(hx-dx)*(hx+dx)*hx32)*((2.0*hy-dy)*(hy-dy)*(3.0*hy-dy)*hy36);

        wi[3][2].i = rx + 1; wi[3][2].x = wi[3][2].i*hx;
        wi[3][2].j = ry + 2; wi[3][2].y = wi[3][2].j*hy;
        dx = fabs(wi[3][2].x-x);
        dy = fabs(wi[3][2].y-y);
        wi[3][2].w = ((2.0*hx-dx)*(hx-dx)*(hx+dx)*hx32)*((2.0*hy-dy)*(hy-dy)*(3.0*hy-dy)*hy36);

        wi[3][3].i = rx + 2; wi[3][3].x = wi[3][3].i*hx;
        wi[3][3].j = ry + 2; wi[3][3].y = wi[3][3].j*hy;
        dx = fabs(wi[3][3].x-x);
        dy = fabs(wi[3][3].y-y);
        wi[3][3].w = ((2.0*hx-dx)*(hx-dx)*(3.0*hx-dx)*hx36)*((2.0*hy-dy)*(hy-dy)*(3.0*hy-dy)*hy36);

        wi[2][3].i = rx + 2; wi[2][3].x = wi[2][3].i*hx;
        wi[2][3].j = ry + 1; wi[2][3].y = wi[2][3].j*hy;
        dx = fabs(wi[2][3].x-x);
        dy = fabs(wi[2][3].y-y);
        wi[2][3].w = ((2.0*hx-dx)*(hx-dx)*(3.0*hx-dx)*hx36)*((2.0*hy-dy)*(hy-dy)*(hy+dy)*hy32);

        wi[1][3].i = rx + 2; wi[1][3].x = wi[1][3].i*hx;
        wi[1][3].j = ry + 0; wi[1][3].y = wi[1][3].j*hy;
        dx = fabs(wi[1][3].x-x);
        dy = fabs(wi[1][3].y-y);
        wi[1][3].w = ((2.0*hx-dx)*(hx-dx)*(3.0*hx-dx)*hx36)*((2.0*hy-dy)*(hy-dy)*(hy+dy)*hy32);

        wi[0][3].i = rx + 2; wi[0][3].x = wi[0][3].i*hx;
        wi[0][3].j = ry - 1; wi[0][3].y = wi[0][3].j*hy;
        dx = fabs(wi[0][3].x-x);
        dy = fabs(wi[0][3].y-y);
        wi[0][3].w = ((2.0*hx-dx)*(hx-dx)*(3.0*hx-dx)*hx36)*((2.0*hy-dy)*(hy-dy)*(3.0*hy-dy)*hy36);

        wi[0][2].i = rx + 1; wi[0][2].x = wi[0][2].i*hx;
        wi[0][2].j = ry - 1; wi[0][2].y = wi[0][2].j*hy;
        dx = fabs(wi[0][2].x-x);
        dy = fabs(wi[0][2].y-y);
        wi[0][2].w = ((2.0*hx-dx)*(hx-dx)*(hx+dx)*hx32)*((2.0*hy-dy)*(hy-dy)*(3.0*hy-dy)*hy36);

        wi[0][1].i = rx + 0; wi[0][1].x = wi[0][1].i*hx;
        wi[0][1].j = ry - 1; wi[0][1].y = wi[0][1].j*hy;
        dx = fabs(wi[0][1].x-x);
        dy = fabs(wi[0][1].y-y);
        wi[0][1].w = ((2.0*hx-dx)*(hx-dx)*(hx+dx)*hx32)*((2.0*hy-dy)*(hy-dy)*(3.0*hy-dy)*hy36);
    }
}

void ExtendedSpaceNode2DH::clearWeights()
{
    for (unsigned int rw=0; rw<rows; rw++)
    {
        for (unsigned int cl=0; cl<cols; cl++)
        {
            delete [] (wi[rw][cl].u);
        }
        delete [] wi[rw];
    }
    delete [] wi;
}

double ExtendedSpaceNode2DH::value(unsigned int layer) const
{
    if (_INFO_ROWS_ == 4 && _INFO_COLS_ == 4) return value1(layer);
    if (_INFO_ROWS_ == 3 && _INFO_COLS_ == 3) return value2(layer);
    return NAN;
}

double ExtendedSpaceNode2DH::value1(unsigned int layer) const
{
    double Lx[] = {0.0, 0.0, 0.0, 0.0};
    double Ly[] = {0.0, 0.0, 0.0, 0.0};

    double y0 = wi[0][0].y;
    double y1 = wi[1][0].y;
    double y2 = wi[2][0].y;
    double y3 = wi[3][0].y;
    Ly[0] = ((y-y1)*(y-y2)*(y-y3))/((y0-y1)*(y0-y2)*(y0-y3));
    Ly[1] = ((y-y0)*(y-y2)*(y-y3))/((y1-y0)*(y1-y2)*(y1-y3));
    Ly[2] = ((y-y0)*(y-y1)*(y-y3))/((y2-y0)*(y2-y1)*(y2-y3));
    Ly[3] = ((y-y0)*(y-y1)*(y-y2))/((y3-y0)*(y3-y1)*(y3-y2));

    double x0 = wi[0][0].x;
    double x1 = wi[0][1].x;
    double x2 = wi[0][2].x;
    double x3 = wi[0][3].x;
    Lx[0] = ((x-x1)*(x-x2)*(x-x3))/((x0-x1)*(x0-x2)*(x0-x3));
    Lx[1] = ((x-x0)*(x-x2)*(x-x3))/((x1-x0)*(x1-x2)*(x1-x3));
    Lx[2] = ((x-x0)*(x-x1)*(x-x3))/((x2-x0)*(x2-x1)*(x2-x3));
    Lx[3] = ((x-x0)*(x-x1)*(x-x2))/((x3-x0)*(x3-x1)*(x3-x2));

    double P = 0.0;
    for (unsigned int j=0; j<rows; j++)
    {
        for (unsigned int i=0; i<cols; i++)
        {
            double u = wi[j][i].u[layer];
            P += Ly[j]*Lx[i]*u;
        }
    }
    return P;
}

double ExtendedSpaceNode2DH::value2(unsigned int layer) const
{
    return wi[1][1].u[layer];
}

double ExtendedSpaceNode2DH::valueDx(unsigned int layer) const
{
    if (_INFO_ROWS_ == 4 && _INFO_COLS_ == 4) return valueDx1(layer);
    if (_INFO_ROWS_ == 3 && _INFO_COLS_ == 3) return valueDx2(layer);
    return NAN;
}

double ExtendedSpaceNode2DH::valueDy(unsigned int layer) const
{
    if (_INFO_ROWS_ == 4 && _INFO_COLS_ == 4) return valueDy1(layer);
    if (_INFO_ROWS_ == 3 && _INFO_COLS_ == 3) return valueDy2(layer);
    return NAN;
}

double ExtendedSpaceNode2DH::valueDx1(unsigned int layer) const
{
    double Lx[] = {0.0, 0.0, 0.0, 0.0};
    double Ly[] = {0.0, 0.0, 0.0, 0.0};

    double y0 = wi[0][0].y;
    double y1 = wi[1][0].y;
    double y2 = wi[2][0].y;
    double y3 = wi[3][0].y;
    Ly[0] = ((y-y1)*(y-y2)*(y-y3))/((y0-y1)*(y0-y2)*(y0-y3));
    Ly[1] = ((y-y0)*(y-y2)*(y-y3))/((y1-y0)*(y1-y2)*(y1-y3));
    Ly[2] = ((y-y0)*(y-y1)*(y-y3))/((y2-y0)*(y2-y1)*(y2-y3));
    Ly[3] = ((y-y0)*(y-y1)*(y-y2))/((y3-y0)*(y3-y1)*(y3-y2));

    double x0 = wi[0][0].x;
    double x1 = wi[0][1].x;
    double x2 = wi[0][2].x;
    double x3 = wi[0][3].x;
    Lx[0] = ((x-x1)*(x-x2)+(x-x2)*(x-x3)+(x-x1)*(x-x3))/((x0-x1)*(x0-x2)*(x0-x3));
    Lx[1] = ((x-x0)*(x-x2)+(x-x2)*(x-x3)+(x-x0)*(x-x3))/((x1-x0)*(x1-x2)*(x1-x3));
    Lx[2] = ((x-x0)*(x-x1)+(x-x1)*(x-x3)+(x-x0)*(x-x3))/((x2-x0)*(x2-x1)*(x2-x3));
    Lx[3] = ((x-x0)*(x-x1)+(x-x1)*(x-x2)+(x-x0)*(x-x2))/((x3-x0)*(x3-x1)*(x3-x2));

    double Px = 0.0;
    for (unsigned int j=0; j<rows; j++)
    {
        for (unsigned int i=0; i<cols; i++)
        {
            double u = wi[j][i].u[layer];
            Px += Ly[j]*Lx[i]*u;
        }
    }
    return Px;


}

double ExtendedSpaceNode2DH::valueDy1(unsigned int layer) const
{
    double Lx[] = {0.0, 0.0, 0.0, 0.0};
    double Ly[] = {0.0, 0.0, 0.0, 0.0};

    double y0 = wi[0][0].y;
    double y1 = wi[1][0].y;
    double y2 = wi[2][0].y;
    double y3 = wi[3][0].y;
    Ly[0] = ((y-y1)*(y-y2)+(y-y2)*(y-y3)+(y-y1)*(y-y3))/((y0-y1)*(y0-y2)*(y0-y3));
    Ly[1] = ((y-y0)*(y-y2)+(y-y2)*(y-y3)+(y-y0)*(y-y3))/((y1-y0)*(y1-y2)*(y1-y3));
    Ly[2] = ((y-y0)*(y-y1)+(y-y1)*(y-y3)+(y-y0)*(y-y3))/((y2-y0)*(y2-y1)*(y2-y3));
    Ly[3] = ((y-y0)*(y-y1)+(y-y1)*(y-y2)+(y-y0)*(y-y2))/((y3-y0)*(y3-y1)*(y3-y2));

    double x0 = wi[0][0].x;
    double x1 = wi[0][1].x;
    double x2 = wi[0][2].x;
    double x3 = wi[0][3].x;
    Lx[0] = ((x-x1)*(x-x2)*(x-x3))/((x0-x1)*(x0-x2)*(x0-x3));
    Lx[1] = ((x-x0)*(x-x2)*(x-x3))/((x1-x0)*(x1-x2)*(x1-x3));
    Lx[2] = ((x-x0)*(x-x1)*(x-x3))/((x2-x0)*(x2-x1)*(x2-x3));
    Lx[3] = ((x-x0)*(x-x1)*(x-x2))/((x3-x0)*(x3-x1)*(x3-x2));

    double Py = 0.0;
    for (unsigned int j=0; j<rows; j++)
    {
        for (unsigned int i=0; i<cols; i++)
        {
            double u = wi[j][i].u[layer];
            Py += Ly[j]*Lx[i]*u;
        }
    }
    return Py;
}

double ExtendedSpaceNode2DH::valueDx2(unsigned int layer) const
{
    return (wi[1][2].u[layer]-wi[1][0].u[layer])/(2.0*0.01);
}

double ExtendedSpaceNode2DH::valueDy2(unsigned int layer) const
{
    return (wi[2][1].u[layer]-wi[0][1].u[layer])/(2.0*0.01);
}

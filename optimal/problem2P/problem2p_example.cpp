#include "problem2p_example.h"

void Problem2PNeumann::Main(int argc UNUSED_PARAM, char** argv UNUSED_PARAM)
{
    //example1();
    //example2();
    //example3();
    example4();
}

auto Problem2PNeumann::checkGradient1(const Problem2PNeumann &prob, const OptimizeParameterP &o_prm) -> void
{
    EquationParameterP e_prm = prob.mEquParameter;
    OptimizeParameterP r_prm = prob.mRegParameter;

    IPrinter::printSeperatorLine();
    DoubleVector pv;
    prob.PrmToVector(o_prm, pv);
    IPrinter::print(pv, pv.length(), 8, 4);
    IPrinter::printSeperatorLine();
    DoubleVector r_pv;
    prob.PrmToVector(r_prm, r_pv);
    IPrinter::print(r_pv, r_pv.length(), 8, 4);
    IPrinter::printSeperatorLine();
    double functional = prob.fx(pv);
    printf("Functional: %f\n", functional);
    puts("Calculating gradients....");
    DoubleVector ag(pv.length());
    prob.gradient(pv, ag);
    puts("Gradients are calculated.");

    DoubleVector ng1(pv.length(), 0.0);
    DoubleVector ng2(pv.length(), 0.0);

    const unsigned int sk = 0*e_prm.Nc*e_prm.No;
    const unsigned int fk = 1*e_prm.Nc*e_prm.No-1;
    const unsigned int sz = 1*e_prm.Nc*e_prm.No;
    const unsigned int fz = 2*e_prm.Nc*e_prm.No-1;
    const unsigned int so = 2*e_prm.Nc*e_prm.No+0*e_prm.No;
    const unsigned int fo = 2*e_prm.Nc*e_prm.No+2*e_prm.No-1;
    const unsigned int sc = 2*e_prm.Nc*e_prm.No+2*e_prm.No;
    const unsigned int fc = 2*e_prm.Nc*e_prm.No+2*e_prm.No+2*e_prm.Nc-1;
#ifdef TIME_DISCRETE
    unsigned int st = 2*e_prm.Nc*e_prm.No + 2*e_prm.No + 2*e_prm.Nc;
    unsigned int ft = 2*e_prm.Nc*e_prm.No + 2*e_prm.No + 2*e_prm.Nc + e_prm.Nt-1;
#endif
    puts("Calculating numerical gradients.... dh=0.01");
    puts("*** Calculating numerical gradients for k...... dh=0.01");
    IGradient::Gradient(&prob, 0.01, pv, ng1, sk, fk);
    puts("*** Calculating numerical gradients for z...... dh=0.01");
    IGradient::Gradient(&prob, 0.01, pv, ng1, sz, fz);
    puts("*** Calculating numerical gradients for xi..... dh=0.01");
    IGradient::Gradient(&prob, 0.01, pv, ng1, so, fo);
    puts("*** Calculating numerical gradients for eta.... dh=0.01");
    IGradient::Gradient(&prob, 0.01, pv, ng1, sc, fc);
#ifdef TIME_DISCRETE
    puts("*** Calculating numerical gradients for tau.... dh=0.01");
    IGradient::Gradient(&prob, 0.01, pv, ng1, st, ft);
#endif
    puts("Numerical gradients are calculated.");

    puts("Calculating numerical gradients.... hx=0.001");
    puts("*** Calculating numerical gradients for k...... dh=0.001");
    IGradient::Gradient(&prob, 0.001, pv, ng2, sk, fk);
    puts("*** Calculating numerical gradients for z...... dh=0.001");
    IGradient::Gradient(&prob, 0.001, pv, ng2, sz, fz);
    puts("*** Calculating numerical gradients for xi..... dh=0.001");
    IGradient::Gradient(&prob, 0.001, pv, ng2, so, fo);
    puts("*** Calculating numerical gradients for eta.... dh=0.001");
    IGradient::Gradient(&prob, 0.001, pv, ng2, sc, fc);
#ifdef TIME_DISCRETE
    puts("*** Calculating numerical gradients for tau.... dh=0.005");
    IGradient::Gradient(&prob, 0.005, pv, ng2, st, ft);
#endif
    puts("Numerical gradients are calculated.");

    //k------------------------------------------------------//
    IPrinter::printSeperatorLine("k");
    DoubleVector pk0 = pv.mid(sk, fk);
    DoubleVector ak0 = ag.mid(sk, fk);
    DoubleVector nk1 = ng1.mid(sk, fk);
    DoubleVector nk2 = ng2.mid(sk, fk);

    IPrinter::print(pk0,pk0.length(),14,4);
    IPrinter::print(ak0,ak0.length(),14,4); ak0.L2Normalize();
    IPrinter::print(nk1,nk1.length(),14,4); nk1.L2Normalize();
    IPrinter::print(nk2,nk2.length(),14,4); nk2.L2Normalize();
    IPrinter::print(ak0,ak0.length(),14,4);
    IPrinter::print(nk1,nk1.length(),14,4);
    IPrinter::print(nk2,nk2.length(),14,4);

    //z------------------------------------------------------//
    IPrinter::printSeperatorLine("z");
    DoubleVector pz0 = pv.mid(sz, fz);
    DoubleVector az0 = ag.mid(sz, fz);
    DoubleVector nz1 = ng1.mid(sz, fz);
    DoubleVector nz2 = ng2.mid(sz, fz);

    IPrinter::print(pz0,pz0.length(),14,4);
    IPrinter::print(az0,az0.length(),14,4); az0.L2Normalize();
    IPrinter::print(nz1,nz1.length(),14,4); nz1.L2Normalize();
    IPrinter::print(nz2,nz2.length(),14,4); nz2.L2Normalize();
    IPrinter::print(az0,az0.length(),14,4);
    IPrinter::print(nz1,nz1.length(),14,4);
    IPrinter::print(nz2,nz2.length(),14,4);

    //xi------------------------------------------------------//
    IPrinter::printSeperatorLine("xi");
    DoubleVector pe0 = pv.mid(so, fo);
    DoubleVector ae0 = ag.mid(so, fo);
    DoubleVector ne1 = ng1.mid(so, fo);
    DoubleVector ne2 = ng2.mid(so, fo);

    IPrinter::print(pe0,pe0.length(),14,4);
    IPrinter::print(ae0,ae0.length(),14,4); ae0.L2Normalize();
    IPrinter::print(ne1,ne1.length(),14,4); ne1.L2Normalize();
    IPrinter::print(ne2,ne2.length(),14,4); ne2.L2Normalize();
    IPrinter::print(ae0,ae0.length(),14,4);
    IPrinter::print(ne1,ne1.length(),14,4);
    IPrinter::print(ne2,ne2.length(),14,4);

    //eta------------------------------------------------------//
    IPrinter::printSeperatorLine("eta");
    DoubleVector px0 = pv.mid(sc, fc);
    DoubleVector ax0 = ag.mid(sc, fc);
    DoubleVector nx1 = ng1.mid(sc, fc);
    DoubleVector nx2 = ng2.mid(sc, fc);

    IPrinter::print(px0,px0.length(),14,4);
    IPrinter::print(ax0,ax0.length(),14,4); ax0.L2Normalize();
    IPrinter::print(nx1,nx1.length(),14,4); nx1.L2Normalize();
    IPrinter::print(nx2,nx2.length(),14,4); nx2.L2Normalize();
    IPrinter::print(ax0,ax0.length(),14,4);
    IPrinter::print(nx1,nx1.length(),14,4);
    IPrinter::print(nx2,nx2.length(),14,4);

#ifdef TIME_DISCRETE
    //tau------------------------------------------------------//
    IPrinter::printSeperatorLine("tau");

    DoubleVector pt0 = pv.mid(st, ft);
    DoubleVector at0 = ag.mid(st, ft);
    DoubleVector nt1 = ng1.mid(st, ft);
    DoubleVector nt2 = ng2.mid(st, ft);

    IPrinter::print(pt0,pt0.length(),14,4);
    IPrinter::print(at0,at0.length(),14,4); at0.L2Normalize();
    IPrinter::print(nt1,nt1.length(),14,4); nt1.L2Normalize();
    IPrinter::print(nt2,nt2.length(),14,4); nt2.L2Normalize();
    IPrinter::print(at0,at0.length(),14,4);
    IPrinter::print(nt1,nt1.length(),14,4);
    IPrinter::print(nt2,nt2.length(),14,4);
#endif
    IPrinter::printSeperatorLine();
}

auto Problem2PNeumann::checkGradient2(const Problem2PNeumann &prob, const OptimizeParameterP &o_prm) -> void
{
    EquationParameterP e_prm = prob.mEquParameter;
    OptimizeParameterP r_prm = prob.mRegParameter;

    IPrinter::printSeperatorLine();
    DoubleVector pv;
    prob.PrmToVector(o_prm, pv);
    IPrinter::print(pv, pv.length(), 8, 4);
    IPrinter::printSeperatorLine();
    DoubleVector r_pv;
    prob.PrmToVector(r_prm, r_pv);
    IPrinter::print(r_pv, r_pv.length(), 8, 4);
    IPrinter::printSeperatorLine();
    DoubleVector ag(pv.length());
    double functional = prob.fx(pv);
    printf("Functional: %f\n", functional);
    puts("Calculating gradients....");
    prob.gradient(pv, ag);
    puts("Gradients are calculated.");

    DoubleVector ng1(pv.length(), 0.0);
    DoubleVector ng2(pv.length(), 0.0);

    const unsigned int sk = 0*e_prm.Nc*e_prm.No;
    const unsigned int fk = 1*e_prm.Nc*e_prm.No-1;
    const unsigned int sz = 1*e_prm.Nc*e_prm.No;
    const unsigned int fz = 2*e_prm.Nc*e_prm.No-1;
    const unsigned int so = 2*e_prm.Nc*e_prm.No+0*e_prm.No;
    const unsigned int fo = 2*e_prm.Nc*e_prm.No+2*e_prm.No-1;
    const unsigned int sc = 2*e_prm.Nc*e_prm.No+2*e_prm.No;
    const unsigned int fc = 2*e_prm.Nc*e_prm.No+2*e_prm.No+2*e_prm.Nc-1;
#ifdef TIME_DISCRETE
    unsigned int st = 2*e_prm.Nc*e_prm.No + 2*e_prm.No + 2*e_prm.Nc;
    unsigned int ft = 2*e_prm.Nc*e_prm.No + 2*e_prm.No + 2*e_prm.Nc + e_prm.Nt-1;
#endif
    puts("Calculating numerical gradients.... dh=0.01");
    puts("*** Calculating numerical gradients for k...... dh=0.01");
    IGradient::Gradient(&prob, 0.01, pv, ng1, sk, fk);
    puts("*** Calculating numerical gradients for z...... dh=0.01");
    IGradient::Gradient(&prob, 0.01, pv, ng1, sz, fz);
    puts("*** Calculating numerical gradients for xi..... dh=0.01");
    IGradient::Gradient(&prob, 0.01, pv, ng1, so, fo);
    puts("*** Calculating numerical gradients for eta.... dh=0.01");
    IGradient::Gradient(&prob, 0.01, pv, ng1, sc, fc);
#ifdef TIME_DISCRETE
    puts("*** Calculating numerical gradients for tau.... dh=0.01");
    IGradient::Gradient(&prob, 0.01, pv, ng1, st, ft);
#endif
    puts("Numerical gradients are calculated.");

    puts("Calculating numerical gradients.... hx=0.001");
    puts("*** Calculating numerical gradients for k...... dh=0.001");
    IGradient::Gradient(&prob, 0.001, pv, ng2, sk, fk);
    puts("*** Calculating numerical gradients for z...... dh=0.001");
    IGradient::Gradient(&prob, 0.001, pv, ng2, sz, fz);
    puts("*** Calculating numerical gradients for xi..... dh=0.001");
    IGradient::Gradient(&prob, 0.001, pv, ng2, so, fo);
    puts("*** Calculating numerical gradients for eta.... dh=0.001");
    IGradient::Gradient(&prob, 0.001, pv, ng2, sc, fc);
#ifdef TIME_DISCRETE
    puts("*** Calculating numerical gradients for tau.... dh=0.005");
    IGradient::Gradient(&prob, 0.005, pv, ng2, st, ft);
#endif
    puts("Numerical gradients are calculated.");

    DoubleVector nag0  = ag; nag0.EuclideanNormalize();
    DoubleVector nng1 = ng1; nng1.EuclideanNormalize();
    DoubleVector nng2 = ng2; nng2.EuclideanNormalize();

    //k------------------------------------------------------//
    IPrinter::printSeperatorLine("k");
    DoubleVector pk0 = pv.mid(sk, fk);
    IPrinter::print(pk0,pk0.length(),14,4);

    DoubleVector ak0 = ag.mid(sk, fk);
    DoubleVector nk1 = ng1.mid(sk, fk);
    DoubleVector nk2 = ng2.mid(sk, fk);
    IPrinter::print(ak0,ak0.length(),14,4);
    IPrinter::print(nk1,nk1.length(),14,4);
    IPrinter::print(nk2,nk2.length(),14,4);

    DoubleVector nak0 = nag0.mid(sk, fk);
    DoubleVector nnk1 = nng1.mid(sk, fk);
    DoubleVector nnk2 = nng2.mid(sk, fk);

    IPrinter::print(nak0,nak0.length(),14,4);
    IPrinter::print(nnk1,nnk1.length(),14,4);
    IPrinter::print(nnk2,nnk2.length(),14,4);

    //z------------------------------------------------------//
    IPrinter::printSeperatorLine("z");
    DoubleVector pz0 = pv.mid(sz, fz);
    IPrinter::print(pz0,pz0.length(),14,4);

    DoubleVector az0 = ag.mid(sz, fz);
    DoubleVector nz1 = ng1.mid(sz, fz);
    DoubleVector nz2 = ng2.mid(sz, fz);
    IPrinter::print(az0,az0.length(),14,4);
    IPrinter::print(nz1,nz1.length(),14,4);
    IPrinter::print(nz2,nz2.length(),14,4);

    DoubleVector naz0 = nag0.mid(sz, fz);
    DoubleVector nnz1 = nng1.mid(sz, fz);
    DoubleVector nnz2 = nng2.mid(sz, fz);
    IPrinter::print(naz0,naz0.length(),14,4);
    IPrinter::print(nnz1,nnz1.length(),14,4);
    IPrinter::print(nnz2,nnz2.length(),14,4);

    //xi------------------------------------------------------//
    IPrinter::printSeperatorLine("xi");
    DoubleVector pe0 = pv.mid(so, fo);
    IPrinter::print(pe0,pe0.length(),14,4);

    DoubleVector ae0 = ag.mid(so, fo);
    DoubleVector ne1 = ng1.mid(so, fo);
    DoubleVector ne2 = ng2.mid(so, fo);
    IPrinter::print(ae0,ae0.length(),14,4);
    IPrinter::print(ne1,ne1.length(),14,4);
    IPrinter::print(ne2,ne2.length(),14,4);

    DoubleVector nae0 = nag0.mid(so, fo);
    DoubleVector nne1 = nng1.mid(so, fo);
    DoubleVector nne2 = nng2.mid(so, fo);
    IPrinter::print(nae0,nae0.length(),14,4);
    IPrinter::print(nne1,nne1.length(),14,4);
    IPrinter::print(nne2,nne2.length(),14,4);

    //eta------------------------------------------------------//
    IPrinter::printSeperatorLine("eta");
    DoubleVector px0 = pv.mid(sc, fc);
    IPrinter::print(px0,px0.length(),14,4);

    DoubleVector ax0 = ag.mid(sc, fc);
    DoubleVector nx1 = ng1.mid(sc, fc);
    DoubleVector nx2 = ng2.mid(sc, fc);
    IPrinter::print(ax0,ax0.length(),14,4);
    IPrinter::print(nx1,nx1.length(),14,4);
    IPrinter::print(nx2,nx2.length(),14,4);

    DoubleVector nax0 = nag0.mid(sc, fc);
    DoubleVector nnx1 = nng1.mid(sc, fc);
    DoubleVector nnx2 = nng2.mid(sc, fc);
    IPrinter::print(nax0,nax0.length(),14,4);
    IPrinter::print(nnx1,nnx1.length(),14,4);
    IPrinter::print(nnx2,nnx2.length(),14,4);

#ifdef TIME_DISCRETE
    //tau------------------------------------------------------//
    IPrinter::printSeperatorLine("tau");
    DoubleVector pt0 = pv.mid(st, ft);
    IPrinter::print(pt0,pt0.length(),14,4);

    DoubleVector at0 = ag.mid(st, ft);
    DoubleVector nt1 = ng1.mid(st, ft);
    DoubleVector nt2 = ng2.mid(st, ft);
    IPrinter::print(at0,at0.length(),14,4);
    IPrinter::print(nt1,nt1.length(),14,4);
    IPrinter::print(nt2,nt2.length(),14,4);

    DoubleVector nat0 = nag0.mid(st, ft);
    DoubleVector nnt1 = nng1.mid(st, ft);
    DoubleVector nnt2 = nng2.mid(st, ft);
    IPrinter::print(nat0,nat0.length(),14,4);
    IPrinter::print(nnt1,nnt1.length(),14,4);
    IPrinter::print(nnt2,nnt2.length(),14,4);
#endif

    IPrinter::printSeperatorLine();
}

auto example1() -> void
{
    EquationParameterP e_prm;
    e_prm.a = 1.0;
    e_prm.alpha = 0.01;
    e_prm.lambda = +0.01;
    e_prm.theta = +6.3;
    e_prm.phi = +0.2;

    e_prm.Nc = 2;
    e_prm.No = 2;

    /*********************************************************************************************************
     * Optimization parameters Table 1
     *********************************************************************************************************/
    OptimizeParameterP o_prm;
    o_prm.k.resize(e_prm.Nc, e_prm.No, 0.0);
    o_prm.z.resize(e_prm.Nc, e_prm.No, 0.0);
    o_prm.xi.resize(e_prm.No);
    o_prm.eta.resize(e_prm.Nc);
    /**********************************************************************************************************
     * Example 1
     **********************************************************************************************************/
    o_prm.k[0][0]  = -5.8500;  o_prm.k[0][1]  = -3.4800;  o_prm.k[1][0]  = -4.7400;  o_prm.k[1][1]  = -6.1500;
    o_prm.z[0][0]  = 14.9100;  o_prm.z[0][1]  = 11.4500;  o_prm.z[1][0]  = 16.8400;  o_prm.z[1][1]  = 12.3800;
    o_prm.xi[0].x  =  0.6900;  o_prm.xi[0].y  =  0.6500;  o_prm.xi[1].x  =  0.4200;  o_prm.xi[1].y  =  0.4700;
    o_prm.eta[0].x =  0.8500;  o_prm.eta[0].y =  0.8600;  o_prm.eta[1].x =  0.2300;  o_prm.eta[1].y =  0.2300;
    /**********************************************************************************************************
     * Table 1, Example 2
     ********************************************************************************************************/
    //o_prm.k[0][0]  = -2.1200;  o_prm.k[0][1]  = +1.2400;  o_prm.k[1][0]  = -2.3800;  o_prm.k[1][1]  = +2.5800;
    //o_prm.z[0][0]  = +8.5000;  o_prm.z[0][1]  = +7.4000;  o_prm.z[1][0]  = +7.7000;  o_prm.z[1][1]  = +9.5000;
    //o_prm.xi[0].x  =  0.6300;  o_prm.xi[0].y  =  0.5200;  o_prm.xi[1].x  =  0.8400;  o_prm.xi[1].y  =  0.6800;
    //o_prm.eta[0].x =  0.4600;  o_prm.eta[0].y =  0.8500;  o_prm.eta[1].x =  0.2400;  o_prm.eta[1].y =  0.2400;
    /********************************************************************************************************/

    /***************************************************************************************************/
    //o_prm.k[0][0] = -0.12; o_prm.k[0][1] = -0.24; o_prm.k[1][0] = -0.38; o_prm.k[1][1] = -0.58;
    //o_prm.z[0][0] = +8.50; o_prm.z[0][1] = +7.40; o_prm.z[1][0] = +7.70; o_prm.z[1][1] = +9.50;
    //o_prm.xi[0].x  = 0.25; o_prm.xi[0].y  = 0.25; o_prm.xi[1].x  = 0.75; o_prm.xi[1].y  = 0.75;
    //o_prm.eta[0].x = 0.25; o_prm.eta[0].y = 0.75; o_prm.eta[1].x = 0.75; o_prm.eta[1].y = 0.25;
    /***************************************************************************************************/

    /*********************************************************************************************************
     * Regularization parameters
     * Table 1, Example 1,2
     *********************************************************************************************************/
    OptimizeParameterP r_prm;
    r_prm.k.resize(e_prm.Nc, e_prm.No, 0.0);
    r_prm.z.resize(e_prm.Nc, e_prm.No, 0.0);
    r_prm.xi.resize(e_prm.No);
    r_prm.eta.resize(e_prm.Nc);
    r_prm.k[0][0]  = -2.0000;  r_prm.k[0][1]  = +0.7000;  r_prm.k[1][0]  = +0.7100;  r_prm.k[1][1]  = -2.3800;
    r_prm.z[0][0]  = +7.9600;  r_prm.z[0][1]  = +5.8300;  r_prm.z[1][0]  = +7.6800;  r_prm.z[1][1]  = +9.7200;
    r_prm.xi[0].x  =  0.0501;  r_prm.xi[0].y  =  0.0501;  r_prm.xi[1].x  =  0.9500;  r_prm.xi[1].y  =  0.0751;
    r_prm.eta[0].x =  0.4149;  r_prm.eta[0].y =  0.7549;  r_prm.eta[1].x =  0.4052;  r_prm.eta[1].y =  0.7077;
    /*********************************************************************************************************/

    // Grid parameters
    double hx = 0.010; int Nx = 100;
    double hy = 0.010; int Ny = 100;
    double ht = 0.005; int Nt = 200;

    Dimension time(ht, 0, Nt);
    Dimension dimx(hx, 0, Nx);
    Dimension dimy(hy, 0, Ny);

    // Penalty paramteres
    DoubleVector r; r << 1.0000;

    // Regularization coefficients
    DoubleVector e; e << 0.0000;

    Problem2PNeumann prob;
    prob.setTimeDimension(time);
    prob.addSpaceDimension(dimx);
    prob.addSpaceDimension(dimy);
    prob.mEquParameter = e_prm;
    prob.mRegParameter = r_prm;
    prob.optimizeK = true;
    prob.optimizeZ = true;
    prob.optimizeC = true;
    prob.optimizeO = true;
    prob.vmin.resize(e_prm.Nc, -5.00);
    prob.vmax.resize(e_prm.Nc, +20.0);
    prob.U.resize(static_cast<const unsigned int>(Ny+1), static_cast<const unsigned int>(Nx+1), 10.0);
    prob.regEpsilon = 0.0;
    prob.r = 0.1;
    prob.checkGradient1(prob, o_prm);
}

auto example2() -> void
{
    EquationParameterP e_prm;
    e_prm.a = 1.0;
    e_prm.alpha = 0.001;
    e_prm.lambda = +0.001;
    e_prm.theta = +0.3;
    e_prm.phi = +0.2;

    e_prm.Nc = 2;
    e_prm.No = 2;

    /*********************************************************************************************************
     * Optimization parameters Table 1
     *********************************************************************************************************/
    OptimizeParameterP o_prm;
    o_prm.k.resize(e_prm.Nc, e_prm.No, -0.5);
    o_prm.z.resize(e_prm.Nc, e_prm.No, 10.0);
    o_prm.xi.resize(e_prm.No);
    o_prm.eta.resize(e_prm.Nc);
    /**********************************************************************************************************
     * Example 1
     **********************************************************************************************************/
    //o_prm.k[0][0]  = -5.8500;  o_prm.k[0][1]  = -3.4800;  o_prm.k[1][0]  = -4.7400;  o_prm.k[1][1]  = -6.1500;
    //o_prm.z[0][0]  = 14.9100;  o_prm.z[0][1]  = 11.4500;  o_prm.z[1][0]  = 16.8400;  o_prm.z[1][1]  = 12.3800;
    o_prm.xi[0].x  =  0.6900;  o_prm.xi[0].y  =  0.6500;  o_prm.xi[1].x  =  0.4200;  o_prm.xi[1].y  =  0.4700;
    o_prm.eta[0].x =  0.8500;  o_prm.eta[0].y =  0.8600;  o_prm.eta[1].x =  0.2300;  o_prm.eta[1].y =  0.2300;
    /**********************************************************************************************************
     * Table 1, Example 2
     ********************************************************************************************************/
    //o_prm.k[0][0]  = -2.1200;  o_prm.k[0][1]  = +1.2400;  o_prm.k[1][0]  = -2.3800;  o_prm.k[1][1]  = +2.5800;
    //o_prm.z[0][0]  = +8.5000;  o_prm.z[0][1]  = +7.4000;  o_prm.z[1][0]  = +7.7000;  o_prm.z[1][1]  = +9.5000;
    //o_prm.xi[0].x  =  0.6300;  o_prm.xi[0].y  =  0.5200;  o_prm.xi[1].x  =  0.8400;  o_prm.xi[1].y  =  0.6800;
    //o_prm.eta[0].x =  0.4600;  o_prm.eta[0].y =  0.8500;  o_prm.eta[1].x =  0.2400;  o_prm.eta[1].y =  0.2400;
    /********************************************************************************************************/

    /***************************************************************************************************/
    //o_prm.k[0][0] = -0.12; o_prm.k[0][1] = -0.24; o_prm.k[1][0] = -0.38; o_prm.k[1][1] = -0.58;
    //o_prm.z[0][0] = +8.50; o_prm.z[0][1] = +7.40; o_prm.z[1][0] = +7.70; o_prm.z[1][1] = +9.50;
    //o_prm.xi[0].x  = 0.25; o_prm.xi[0].y  = 0.25; o_prm.xi[1].x  = 0.75; o_prm.xi[1].y  = 0.75;
    //o_prm.eta[0].x = 0.25; o_prm.eta[0].y = 0.75; o_prm.eta[1].x = 0.75; o_prm.eta[1].y = 0.25;
    /***************************************************************************************************/

    /*********************************************************************************************************
     * Regularization parameters
     * Table 1, Example 1,2
     *********************************************************************************************************/
    OptimizeParameterP r_prm;
    r_prm.k.resize(e_prm.Nc, e_prm.No, 0.0);
    r_prm.z.resize(e_prm.Nc, e_prm.No, 0.0);
    r_prm.xi.resize(e_prm.No);
    r_prm.eta.resize(e_prm.Nc);
    r_prm.k[0][0]  = -2.0000;  r_prm.k[0][1]  = +0.7000;  r_prm.k[1][0]  = +0.7100;  r_prm.k[1][1]  = -2.3800;
    r_prm.z[0][0]  = +7.9600;  r_prm.z[0][1]  = +5.8300;  r_prm.z[1][0]  = +7.6800;  r_prm.z[1][1]  = +9.7200;
    r_prm.xi[0].x  =  0.0501;  r_prm.xi[0].y  =  0.0501;  r_prm.xi[1].x  =  0.9500;  r_prm.xi[1].y  =  0.0751;
    r_prm.eta[0].x =  0.4149;  r_prm.eta[0].y =  0.7549;  r_prm.eta[1].x =  0.4052;  r_prm.eta[1].y =  0.7077;
    /*********************************************************************************************************/

    // Grid parameters
    double hx = 0.010; int Nx = 100;
    double hy = 0.010; int Ny = 100;
    double ht = 0.005; int Nt = 200;

    Dimension time(ht, 0, Nt);
    Dimension dimx(hx, 0, Nx);
    Dimension dimy(hy, 0, Ny);

    // Penalty paramteres
    DoubleVector r; r << 0.0000 << 20.000 << 50.0000 << 100.00;

    // Regularization coefficients
    DoubleVector e; e << 0.0000 << 0.0000 << 0.00000 << 0.0000;

    DoubleVector x;
    for (unsigned int i=0; i<r.length(); i++)
    {
        Problem2PNeumann prob;
        prob.setTimeDimension(time);
        prob.addSpaceDimension(dimx);
        prob.addSpaceDimension(dimy);
        prob.mEquParameter = e_prm;
        prob.mRegParameter = r_prm;
        prob.optimizeK = true;
        prob.optimizeZ = true;
        prob.optimizeC = true;
        prob.optimizeO = true;
        prob.vmin.resize(e_prm.Nc,  -0.00);
        prob.vmax.resize(e_prm.Nc, +25.00);
        prob.U.resize(static_cast<const unsigned int>(Ny+1), static_cast<const unsigned int>(Nx+1), 10.0);

        prob.regEpsilon = e[i];
        prob.r = r[i];
        if (i==0)
        {
            prob.PrmToVector(o_prm, x);
            //prob.checkGradient1(prob, o_prm);
            IPrinter::printSeperatorLine();
        }

        //ConjugateGradient g;
        SteepestDescentGradient g;
        g.setFunction(&prob);
        g.setGradient(&prob);
        g.setPrinter(&prob);
        g.setProjection(&prob);
        //g.setGradientNormalizer(&prob);
        g.setOptimalityTolerance(0.00001);
        g.setStepTolerance(0.00001);
        g.setFunctionTolerance(0.00001);
        g.setR1MinimizeEpsilon(0.1, 0.0001);
        g.setNormalize(true);
        g.showExitMessage(true);

        g.calculate(x);

        IPrinter::printSeperatorLine(NULL, '=');
    }
}

auto example3() -> void
{
    EquationParameterP e_prm;
    e_prm.a = 1.0;
    e_prm.alpha = 0.001;
    e_prm.lambda = +0.001;
    e_prm.theta = +0.3;
    e_prm.phi = +0.2;

    e_prm.Nc = 2;
    e_prm.No = 2;

    /*********************************************************************************************************
     * Optimization parameters Table 1
     *********************************************************************************************************/
    OptimizeParameterP o_prm;
    o_prm.k.resize(e_prm.Nc, e_prm.No, -0.5);
    o_prm.z.resize(e_prm.Nc, e_prm.No, 10.0);
    o_prm.xi.resize(e_prm.No);
    o_prm.eta.resize(e_prm.Nc);
    /**********************************************************************************************************
     * Example 1
     **********************************************************************************************************/
    //o_prm.k[0][0]  = -5.8500;  o_prm.k[0][1]  = -3.4800;  o_prm.k[1][0]  = -4.7400;  o_prm.k[1][1]  = -6.1500;
    //o_prm.z[0][0]  = 14.9100;  o_prm.z[0][1]  = 11.4500;  o_prm.z[1][0]  = 16.8400;  o_prm.z[1][1]  = 12.3800;
    o_prm.xi[0].x  =  0.6900;  o_prm.xi[0].y  =  0.6500;  o_prm.xi[1].x  =  0.4200;  o_prm.xi[1].y  =  0.4700;
    o_prm.eta[0].x =  0.8500;  o_prm.eta[0].y =  0.8600;  o_prm.eta[1].x =  0.2300;  o_prm.eta[1].y =  0.2300;
    /**********************************************************************************************************
     * Table 1, Example 2
     ********************************************************************************************************/
    //o_prm.k[0][0]  = -2.1200;  o_prm.k[0][1]  = +1.2400;  o_prm.k[1][0]  = -2.3800;  o_prm.k[1][1]  = +2.5800;
    //o_prm.z[0][0]  = +8.5000;  o_prm.z[0][1]  = +7.4000;  o_prm.z[1][0]  = +7.7000;  o_prm.z[1][1]  = +9.5000;
    //o_prm.xi[0].x  =  0.6300;  o_prm.xi[0].y  =  0.5200;  o_prm.xi[1].x  =  0.8400;  o_prm.xi[1].y  =  0.6800;
    //o_prm.eta[0].x =  0.4600;  o_prm.eta[0].y =  0.8500;  o_prm.eta[1].x =  0.2400;  o_prm.eta[1].y =  0.2400;
    /********************************************************************************************************/

    /***************************************************************************************************/
    //o_prm.k[0][0] = -0.12; o_prm.k[0][1] = -0.24; o_prm.k[1][0] = -0.38; o_prm.k[1][1] = -0.58;
    //o_prm.z[0][0] = +8.50; o_prm.z[0][1] = +7.40; o_prm.z[1][0] = +7.70; o_prm.z[1][1] = +9.50;
    //o_prm.xi[0].x  = 0.25; o_prm.xi[0].y  = 0.25; o_prm.xi[1].x  = 0.75; o_prm.xi[1].y  = 0.75;
    //o_prm.eta[0].x = 0.25; o_prm.eta[0].y = 0.75; o_prm.eta[1].x = 0.75; o_prm.eta[1].y = 0.25;
    /***************************************************************************************************/

    /*********************************************************************************************************
     * Regularization parameters
     * Table 1, Example 1,2
     *********************************************************************************************************/
    OptimizeParameterP r_prm;
    r_prm.k.resize(e_prm.Nc, e_prm.No, 0.0);
    r_prm.z.resize(e_prm.Nc, e_prm.No, 0.0);
    r_prm.xi.resize(e_prm.No);
    r_prm.eta.resize(e_prm.Nc);
    r_prm.k[0][0]  = -2.0000;  r_prm.k[0][1]  = +0.7000;  r_prm.k[1][0]  = +0.7100;  r_prm.k[1][1]  = -2.3800;
    r_prm.z[0][0]  = +7.9600;  r_prm.z[0][1]  = +5.8300;  r_prm.z[1][0]  = +7.6800;  r_prm.z[1][1]  = +9.7200;
    r_prm.xi[0].x  =  0.0501;  r_prm.xi[0].y  =  0.0501;  r_prm.xi[1].x  =  0.9500;  r_prm.xi[1].y  =  0.0751;
    r_prm.eta[0].x =  0.4149;  r_prm.eta[0].y =  0.7549;  r_prm.eta[1].x =  0.4052;  r_prm.eta[1].y =  0.7077;

    DoubleVector rv;

    //I[ 13]: F:  0.000001
    rv << -1.1964 << -1.2618 << -1.2777 << -1.3481
       << 10.0071 << 10.0033 << 10.0071 << 10.0026
       <<  0.9500 <<  0.3487 <<  0.2838 <<  0.6895
        <<  0.8935 <<  0.9469 <<  0.0500 <<  0.0721;
    /*********************************************************************************************************/

    // Grid parameters
    double hx = 0.010; int Nx = 100;
    double hy = 0.010; int Ny = 100;
    double ht = 0.005; int Nt = 200;

    Dimension time(ht, 0, Nt);
    Dimension dimx(hx, 0, Nx);
    Dimension dimy(hy, 0, Ny);

    // Penalty paramteres
    DoubleVector r; r << 0.0000;// << 20.000 << 50.0000 << 100.00;

    // Regularization coefficients
    DoubleVector e; e << 0.0000;// << 0.0000 << 0.00000 << 0.0000;

    DoubleVector x;
    for (unsigned int i=0; i<r.length(); i++)
    {
        Problem2PNeumann prob;
        prob.setTimeDimension(time);
        prob.addSpaceDimension(dimx);
        prob.addSpaceDimension(dimy);
        prob.mEquParameter = e_prm;
        prob.mRegParameter = r_prm;
        prob.optimizeK = true;
        prob.optimizeZ = true;
        prob.optimizeC = true;
        prob.optimizeO = true;
        prob.vmin.resize(e_prm.Nc, -5.00);
        prob.vmax.resize(e_prm.Nc, +5.00);
        prob.U.resize(static_cast<const unsigned int>(Ny+1), static_cast<const unsigned int>(Nx+1), 10.0);

        prob.VectorToPrm(rv, r_prm);
        prob.mRegParameter = r_prm;

        prob.regEpsilon = e[i];
        prob.r = r[i];
        if (i==0)
        {
            prob.PrmToVector(o_prm, x);
            //prob.checkGradient1(prob, o_prm);
            IPrinter::printSeperatorLine();
        }

        //ConjugateGradient g;
        SteepestDescentGradient g;
        g.setFunction(&prob);
        g.setGradient(&prob);
        g.setPrinter(&prob);
        g.setProjection(&prob);
        //g.setGradientNormalizer(&prob);
        g.setOptimalityTolerance(0.00001);
        g.setStepTolerance(0.00001);
        g.setFunctionTolerance(0.00001);
        g.setR1MinimizeEpsilon(0.1, 0.0001);
        g.setNormalize(true);
        g.showExitMessage(true);
        g.setMaxIterations(100);

        g.calculate(x);

        IPrinter::printSeperatorLine(NULL, '=');
    }
}

auto example4() -> void
{
    EquationParameterP e_prm;
    e_prm.a = 1.0;
    e_prm.alpha = 0.01;
    e_prm.lambda = +0.01;
    e_prm.theta = +6.3;
    e_prm.phi = +0.2;

    e_prm.Nc = 2;
    e_prm.No = 2;

    /*********************************************************************************************************
     * Optimization parameters Table 1
     *********************************************************************************************************/
    OptimizeParameterP o_prm;
    o_prm.k.resize(e_prm.Nc, e_prm.No, 0.0);
    o_prm.z.resize(e_prm.Nc, e_prm.No, 0.0);
    o_prm.xi.resize(e_prm.No);
    o_prm.eta.resize(e_prm.Nc);
    /**********************************************************************************************************
     * Example 1
     **********************************************************************************************************/
    //o_prm.k[0][0]  = -0.0585;  o_prm.k[0][1]  = -0.0348;  o_prm.k[1][0]  = -0.0474;  o_prm.k[1][1]  = -0.0615;
    o_prm.k[0][0]  = -0.5850;  o_prm.k[0][1]  = -0.3480;  o_prm.k[1][0]  = -0.4740;  o_prm.k[1][1]  = -0.6150;
    //o_prm.k[0][0]  = -5.8500;  o_prm.k[0][1]  = -3.4800;  o_prm.k[1][0]  = -4.7400;  o_prm.k[1][1]  = -6.1500;
    o_prm.z[0][0]  = 14.9100;  o_prm.z[0][1]  = 11.4500;  o_prm.z[1][0]  = 16.8400;  o_prm.z[1][1]  = 12.3800;
    o_prm.xi[0].x  =  0.6900;  o_prm.xi[0].y  =  0.6500;  o_prm.xi[1].x  =  0.4200;  o_prm.xi[1].y  =  0.4700;
    o_prm.eta[0].x =  0.8500;  o_prm.eta[0].y =  0.8600;  o_prm.eta[1].x =  0.2300;  o_prm.eta[1].y =  0.2300;
    /**********************************************************************************************************
     * Table 1, Example 2
     ********************************************************************************************************/
    //o_prm.k[0][0]  = -2.1200;  o_prm.k[0][1]  = +1.2400;  o_prm.k[1][0]  = -2.3800;  o_prm.k[1][1]  = +2.5800;
    //o_prm.z[0][0]  = +8.5000;  o_prm.z[0][1]  = +7.4000;  o_prm.z[1][0]  = +7.7000;  o_prm.z[1][1]  = +9.5000;
    //o_prm.xi[0].x  =  0.6300;  o_prm.xi[0].y  =  0.5200;  o_prm.xi[1].x  =  0.8400;  o_prm.xi[1].y  =  0.6800;
    //o_prm.eta[0].x =  0.4600;  o_prm.eta[0].y =  0.8500;  o_prm.eta[1].x =  0.2400;  o_prm.eta[1].y =  0.2400;
    /********************************************************************************************************/

    /***************************************************************************************************/
    //o_prm.k[0][0] = -0.12; o_prm.k[0][1] = -0.24; o_prm.k[1][0] = -0.38; o_prm.k[1][1] = -0.58;
    //o_prm.z[0][0] = +8.50; o_prm.z[0][1] = +7.40; o_prm.z[1][0] = +7.70; o_prm.z[1][1] = +9.50;
    //o_prm.xi[0].x  = 0.25; o_prm.xi[0].y  = 0.25; o_prm.xi[1].x  = 0.75; o_prm.xi[1].y  = 0.75;
    //o_prm.eta[0].x = 0.25; o_prm.eta[0].y = 0.75; o_prm.eta[1].x = 0.75; o_prm.eta[1].y = 0.25;
    /***************************************************************************************************/

    /*********************************************************************************************************
     * Regularization parameters
     * Table 1, Example 1,2
     *********************************************************************************************************/
    OptimizeParameterP r_prm;
    r_prm.k.resize(e_prm.Nc, e_prm.No, 0.0);
    r_prm.z.resize(e_prm.Nc, e_prm.No, 0.0);
    r_prm.xi.resize(e_prm.No);
    r_prm.eta.resize(e_prm.Nc);
    r_prm.k[0][0]  = -2.0000;  r_prm.k[0][1]  = +0.7000;  r_prm.k[1][0]  = +0.7100;  r_prm.k[1][1]  = -2.3800;
    r_prm.z[0][0]  = +7.9600;  r_prm.z[0][1]  = +5.8300;  r_prm.z[1][0]  = +7.6800;  r_prm.z[1][1]  = +9.7200;
    r_prm.xi[0].x  =  0.0501;  r_prm.xi[0].y  =  0.0501;  r_prm.xi[1].x  =  0.9500;  r_prm.xi[1].y  =  0.0751;
    r_prm.eta[0].x =  0.4149;  r_prm.eta[0].y =  0.7549;  r_prm.eta[1].x =  0.4052;  r_prm.eta[1].y =  0.7077;
    /*********************************************************************************************************/

    // Grid parameters
    double hx = 0.010; int Nx = 100;
    double hy = 0.010; int Ny = 100;
    double ht = 0.005; int Nt = 200;

    Dimension time(ht, 0, Nt);
    Dimension dimx(hx, 0, Nx);
    Dimension dimy(hy, 0, Ny);

    // Penalty paramteres
    //DoubleVector r; r << 0.0000;

    // Regularization coefficients
    //DoubleVector e; e << 0.0000;

    Problem2PNeumann prob;

#ifdef TIME_DISCRETE
    e_prm.Nt = 3;

    o_prm.tau.push_back(0.10);
    o_prm.tau.push_back(0.40);
    o_prm.tau.push_back(0.80);

    r_prm.tau.push_back(0.25);
    r_prm.tau.push_back(0.50);
    r_prm.tau.push_back(0.75);
#endif

    prob.setTimeDimension(time);
    prob.addSpaceDimension(dimx);
    prob.addSpaceDimension(dimy);
    prob.mEquParameter = e_prm;
    prob.mRegParameter = r_prm;
    prob.optimizeK = true;
    prob.optimizeZ = true;
    prob.optimizeC = true;
    prob.optimizeO = true;
#ifdef TIME_DISCRETE
    prob.optimezeT = true;
#endif
    prob.vmin.resize(e_prm.Nc, -5.0);
    prob.vmax.resize(e_prm.Nc, +5.0);
    prob.U.resize(static_cast<const unsigned int>(Ny+1), static_cast<const unsigned int>(Nx+1), 10.0);
    prob.regEpsilon = 0.0;
    prob.r = 0.0;
    prob.checkGradient1(prob, o_prm);

    DoubleMatrix u;
    spif_vector u_info;
    prob.solveForwardIBVP(u, u_info, true, o_prm);
    const SpacePointInfoP &ui1 = u_info[0];
    const SpacePointInfoP &ui2 = u_info[1];
    for (unsigned int i=0; i<ui1.length; i+=2)
    {
        printf("%6f %10.6f %10.6f\n", (i/2)*0.005, ui1.vl[i], ui2.vl[i]);
    }
}



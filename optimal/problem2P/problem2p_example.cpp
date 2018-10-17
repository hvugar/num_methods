#include "problem2p_example.h"

void Problem2PNeumann::Main(int argc UNUSED_PARAM, char** argv UNUSED_PARAM)
{
    example1();
}

auto Problem2PNeumann::checkGradient1(const Problem2PNeumann &prob, const OptimizeParameter &o_prm) -> void
{
    EquationParameter e_prm = prob.mEquParameter;
    OptimizeParameter r_prm = prob.mRegParameter;

    IPrinter::printSeperatorLine();
    DoubleVector pv;
    prob.PrmToVector(o_prm, pv);
    IPrinter::print(pv, pv.length(), 6, 4);
    IPrinter::printSeperatorLine();
    DoubleVector r_pv;
    prob.PrmToVector(r_prm, r_pv);
    IPrinter::print(r_pv, r_pv.length(), 6, 4);
    IPrinter::printSeperatorLine();
    double functional = prob.fx(pv);
    printf("Functional: %f\n", functional);
    puts("Calculating gradients....");
    DoubleVector ag(pv.length());
    prob.gradient(pv, ag);
    puts("Gradients are calculated.");

    DoubleVector ng1(pv.length(), 0.0);
    DoubleVector ng2(pv.length(), 0.0);

    puts("Calculating numerical gradients.... dh=0.01");
    puts("*** Calculating numerical gradients for k...... dh=0.01");
    IGradient::Gradient(&prob, 0.01, pv, ng1, 0*e_prm.Nc*e_prm.No,            1*e_prm.Nc*e_prm.No-1);
    puts("*** Calculating numerical gradients for z...... dh=0.01");
    IGradient::Gradient(&prob, 0.01, pv, ng1, 1*e_prm.Nc*e_prm.No,            2*e_prm.Nc*e_prm.No-1);
    puts("*** Calculating numerical gradients for xi..... dh=0.01");
    IGradient::Gradient(&prob, 0.01, pv, ng1, 2*e_prm.Nc*e_prm.No+0*e_prm.No, 2*e_prm.Nc*e_prm.No+2*e_prm.No-1);
    puts("*** Calculating numerical gradients for eta.... dh=0.01");
    IGradient::Gradient(&prob, 0.01, pv, ng1, 2*e_prm.Nc*e_prm.No+2*e_prm.No, 2*e_prm.Nc*e_prm.No+2*e_prm.No+2*e_prm.Nc-1);
    puts("Numerical gradients are calculated.");

    puts("Calculating numerical gradients.... hx=0.001");
    puts("*** Calculating numerical gradients for k...... dh=0.001");
    IGradient::Gradient(&prob, 0.001, pv, ng2, 0*e_prm.Nc*e_prm.No,            1*e_prm.Nc*e_prm.No-1);
    puts("*** Calculating numerical gradients for z...... dh=0.001");
    IGradient::Gradient(&prob, 0.001, pv, ng2, 1*e_prm.Nc*e_prm.No,            2*e_prm.Nc*e_prm.No-1);
    puts("*** Calculating numerical gradients for xi..... dh=0.001");
    IGradient::Gradient(&prob, 0.001, pv, ng2, 2*e_prm.Nc*e_prm.No+0*e_prm.No, 2*e_prm.Nc*e_prm.No+2*e_prm.No-1);
    puts("*** Calculating numerical gradients for eta.... dh=0.001");
    IGradient::Gradient(&prob, 0.001, pv, ng2, 2*e_prm.Nc*e_prm.No+2*e_prm.No, 2*e_prm.Nc*e_prm.No+2*e_prm.No+2*e_prm.Nc-1);
    puts("Numerical gradients are calculated.");

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

auto Problem2PNeumann::checkGradient2(const Problem2PNeumann &prob, const OptimizeParameter &o_prm) -> void
{
    EquationParameter e_prm = prob.mEquParameter;
    //OptimizeParameter o_prm = prob.mOptParameter;
    OptimizeParameter r_prm = prob.mRegParameter;

    IPrinter::printSeperatorLine();
    DoubleVector pv;
    prob.PrmToVector(o_prm, pv);
    IPrinter::print(pv, pv.length(), 6, 4);
    IPrinter::printSeperatorLine();
    DoubleVector r_pv;
    prob.PrmToVector(r_prm, r_pv);
    IPrinter::print(r_pv, r_pv.length(), 6, 4);
    IPrinter::printSeperatorLine();
    DoubleVector ag(pv.length());
    double functional = prob.fx(pv);
    printf("Functional: %f\n", functional);
    puts("Calculating gradients....");
    prob.gradient(pv, ag);
    puts("Gradients are calculated.");

    DoubleVector ng1(pv.length(), 0.0);
    DoubleVector ng2(pv.length(), 0.0);

    puts("Calculating numerical gradients.... dh=0.01");
    puts("*** Calculating numerical gradients for k...... dh=0.01");
    IGradient::Gradient(&prob, 0.01, pv, ng1, 0*e_prm.Nc*e_prm.No,            1*e_prm.Nc*e_prm.No-1);
    puts("*** Calculating numerical gradients for z...... dh=0.01");
    IGradient::Gradient(&prob, 0.01, pv, ng1, 1*e_prm.Nc*e_prm.No,            2*e_prm.Nc*e_prm.No-1);
    puts("*** Calculating numerical gradients for xi..... dh=0.01");
    IGradient::Gradient(&prob, 0.01, pv, ng1, 2*e_prm.Nc*e_prm.No+0*e_prm.No, 2*e_prm.Nc*e_prm.No+2*e_prm.No-1);
    puts("*** Calculating numerical gradients for eta.... dh=0.01");
    IGradient::Gradient(&prob, 0.01, pv, ng1, 2*e_prm.Nc*e_prm.No+2*e_prm.No, 2*e_prm.Nc*e_prm.No+2*e_prm.No+2*e_prm.Nc-1);
    puts("Numerical gradients are calculated.");

    puts("Calculating numerical gradients.... hx=0.001");
    puts("*** Calculating numerical gradients for k...... dh=0.001");
    IGradient::Gradient(&prob, 0.001, pv, ng2, 0*e_prm.Nc*e_prm.No,            1*e_prm.Nc*e_prm.No-1);
    puts("*** Calculating numerical gradients for z...... dh=0.001");
    IGradient::Gradient(&prob, 0.001, pv, ng2, 1*e_prm.Nc*e_prm.No,            2*e_prm.Nc*e_prm.No-1);
    puts("*** Calculating numerical gradients for xi..... dh=0.001");
    IGradient::Gradient(&prob, 0.001, pv, ng2, 2*e_prm.Nc*e_prm.No+0*e_prm.No, 2*e_prm.Nc*e_prm.No+2*e_prm.No-1);
    puts("*** Calculating numerical gradients for eta.... dh=0.001");
    IGradient::Gradient(&prob, 0.001, pv, ng2, 2*e_prm.Nc*e_prm.No+2*e_prm.No, 2*e_prm.Nc*e_prm.No+2*e_prm.No+2*e_prm.Nc-1);
    puts("Numerical gradients are calculated.");

    DoubleVector nag  = ag;  nag.EuclideanNormalize();
    DoubleVector nng1 = ng1; nng1.EuclideanNormalize();
    DoubleVector nng2 = ng2; nng2.EuclideanNormalize();

    //k------------------------------------------------------//
    IPrinter::printSeperatorLine("k");
    DoubleVector pk0 = pv.mid(0, e_prm.Nc*e_prm.No-1);
    IPrinter::print(pk0,pk0.length(),14,4);

    DoubleVector ak0 = ag.mid(0, e_prm.Nc*e_prm.No-1);
    DoubleVector nk1 = ng1.mid(0, e_prm.Nc*e_prm.No-1);
    DoubleVector nk2 = ng2.mid(0, e_prm.Nc*e_prm.No-1);
    IPrinter::print(ak0,ak0.length(),14,4);
    IPrinter::print(nk1,nk1.length(),14,4);
    IPrinter::print(nk2,nk2.length(),14,4);

    DoubleVector nak0 = nag.mid(0, e_prm.Nc*e_prm.No-1);
    DoubleVector nnk1 = nng1.mid(0, e_prm.Nc*e_prm.No-1);
    DoubleVector nnk2 = nng2.mid(0, e_prm.Nc*e_prm.No-1);

    IPrinter::print(nak0,nak0.length(),14,4);
    IPrinter::print(nnk1,nnk1.length(),14,4);
    IPrinter::print(nnk2,nnk2.length(),14,4);

    //z------------------------------------------------------//
    IPrinter::printSeperatorLine("z");
    DoubleVector pz0 = pv.mid(e_prm.Nc*e_prm.No, 2*e_prm.Nc*e_prm.No-1);
    IPrinter::print(pz0,pz0.length(),14,4);

    DoubleVector az0 = ag.mid(e_prm.Nc*e_prm.No, 2*e_prm.Nc*e_prm.No-1);
    DoubleVector nz1 = ng1.mid(e_prm.Nc*e_prm.No, 2*e_prm.Nc*e_prm.No-1);
    DoubleVector nz2 = ng2.mid(e_prm.Nc*e_prm.No, 2*e_prm.Nc*e_prm.No-1);
    IPrinter::print(az0,az0.length(),14,4);
    IPrinter::print(nz1,nz1.length(),14,4);
    IPrinter::print(nz2,nz2.length(),14,4);

    DoubleVector naz0 = nag.mid(e_prm.Nc*e_prm.No, 2*e_prm.Nc*e_prm.No-1);
    DoubleVector nnz1 = nng1.mid(e_prm.Nc*e_prm.No, 2*e_prm.Nc*e_prm.No-1);
    DoubleVector nnz2 = nng2.mid(e_prm.Nc*e_prm.No, 2*e_prm.Nc*e_prm.No-1);
    IPrinter::print(naz0,naz0.length(),14,4);
    IPrinter::print(nnz1,nnz1.length(),14,4);
    IPrinter::print(nnz2,nnz2.length(),14,4);

    //xi------------------------------------------------------//
    IPrinter::printSeperatorLine("xi");
    DoubleVector pe0 = pv.mid(2*e_prm.Nc*e_prm.No, 2*e_prm.Nc*e_prm.No+2*e_prm.No-1);
    IPrinter::print(pe0,pe0.length(),14,4);

    DoubleVector ae0 = ag.mid(2*e_prm.Nc*e_prm.No, 2*e_prm.Nc*e_prm.No+2*e_prm.No-1);
    DoubleVector ne1 = ng1.mid(2*e_prm.Nc*e_prm.No, 2*e_prm.Nc*e_prm.No+2*e_prm.No-1);
    DoubleVector ne2 = ng2.mid(2*e_prm.Nc*e_prm.No, 2*e_prm.Nc*e_prm.No+2*e_prm.No-1);
    IPrinter::print(ae0,ae0.length(),14,4);
    IPrinter::print(ne1,ne1.length(),14,4);
    IPrinter::print(ne2,ne2.length(),14,4);

    DoubleVector nae0 = nag.mid(2*e_prm.Nc*e_prm.No, 2*e_prm.Nc*e_prm.No+2*e_prm.No-1);
    DoubleVector nne1 = nng1.mid(2*e_prm.Nc*e_prm.No, 2*e_prm.Nc*e_prm.No+2*e_prm.No-1);
    DoubleVector nne2 = nng2.mid(2*e_prm.Nc*e_prm.No, 2*e_prm.Nc*e_prm.No+2*e_prm.No-1);
    IPrinter::print(nae0,nae0.length(),14,4);
    IPrinter::print(nne1,nne1.length(),14,4);
    IPrinter::print(nne2,nne2.length(),14,4);

    //eta------------------------------------------------------//
    IPrinter::printSeperatorLine("eta");
    DoubleVector px0 = pv.mid(2*e_prm.Nc*e_prm.No+2*e_prm.No, 2*e_prm.Nc*e_prm.No+2*e_prm.No+2*e_prm.Nc-1);
    IPrinter::print(px0,px0.length(),14,4);

    DoubleVector ax0 = ag.mid(2*e_prm.Nc*e_prm.No+2*e_prm.No, 2*e_prm.Nc*e_prm.No+2*e_prm.No+2*e_prm.Nc-1);
    DoubleVector nx1 = ng1.mid(2*e_prm.Nc*e_prm.No+2*e_prm.No, 2*e_prm.Nc*e_prm.No+2*e_prm.No+2*e_prm.Nc-1);
    DoubleVector nx2 = ng2.mid(2*e_prm.Nc*e_prm.No+2*e_prm.No, 2*e_prm.Nc*e_prm.No+2*e_prm.No+2*e_prm.Nc-1);
    IPrinter::print(ax0,ax0.length(),14,4);
    IPrinter::print(nx1,nx1.length(),14,4);
    IPrinter::print(nx2,nx2.length(),14,4);

    DoubleVector nax0 = nag.mid(2*e_prm.Nc*e_prm.No+2*e_prm.No, 2*e_prm.Nc*e_prm.No+2*e_prm.No+2*e_prm.Nc-1);
    DoubleVector nnx1 = nng1.mid(2*e_prm.Nc*e_prm.No+2*e_prm.No, 2*e_prm.Nc*e_prm.No+2*e_prm.No+2*e_prm.Nc-1);
    DoubleVector nnx2 = nng2.mid(2*e_prm.Nc*e_prm.No+2*e_prm.No, 2*e_prm.Nc*e_prm.No+2*e_prm.No+2*e_prm.Nc-1);
    IPrinter::print(nax0,nax0.length(),14,4);
    IPrinter::print(nnx1,nnx1.length(),14,4);
    IPrinter::print(nnx2,nnx2.length(),14,4);
    IPrinter::printSeperatorLine();
}

auto example1() -> void
{
    EquationParameter e_prm;
    e_prm.a = 1.0;
    e_prm.alpha = 0.01;
    e_prm.lambda = +0.01;
    e_prm.theta = +6.3;
    e_prm.phi = +0.2;

    e_prm.Nc = 2;
    e_prm.No = 2;

    // Optimization parameters
    OptimizeParameter o_prm;
    o_prm.k.resize(e_prm.Nc, e_prm.No, 0.0);
    o_prm.z.resize(e_prm.Nc, e_prm.No, 0.0);
    o_prm.xi.resize(e_prm.No);
    o_prm.eta.resize(e_prm.Nc);

    /***************************************************************************************************/
    //o_prm.k[0][0]  =  -5.85;  o_prm.k[0][1]  =  -3.48;  o_prm.k[1][0]  =  -4.74;  o_prm.k[1][1] =  -6.15;
    //o_prm.z[0][0]  = +14.91;  o_prm.z[0][1]  = +11.45;  o_prm.z[1][0]  = +16.84;  o_prm.z[1][1] = +12.38;
    //o_prm.xi[0].x  =   0.69;  o_prm.xi[0].y  =   0.65;  o_prm.xi[1].x  =   0.42;  o_prm.xi[1].y =   0.47;
    //o_prm.eta[0].x =   0.85;  o_prm.eta[0].y =   0.86;  o_prm.eta[1].x =   0.23;  o_prm.eta[1].y =  0.23;
    /***************************************************************************************************/

    /***************************************************************************************************/
    o_prm.k[0][0] = -2.12; o_prm.k[0][1] = +1.24; o_prm.k[1][0] = -2.38; o_prm.k[1][1] = +2.58;
    o_prm.z[0][0] = +8.50; o_prm.z[0][1] = +7.40; o_prm.z[1][0] = +7.70; o_prm.z[1][1] = +9.50;
    o_prm.xi[0].x =  0.63; o_prm.xi[0].y =  0.52; o_prm.xi[1].x =  0.84; o_prm.xi[1].y =  0.68;
    o_prm.eta[0].x = 0.46; o_prm.eta[0].y = 0.85; o_prm.eta[1].x = 0.24; o_prm.eta[1].y = 0.24;
    /***************************************************************************************************/

    /***************************************************************************************************/
    //o_prm.k[0][0] = -0.12; o_prm.k[0][1] = -0.24; o_prm.k[1][0] = -0.38; o_prm.k[1][1] = -0.58;
    //o_prm.z[0][0] = +8.50; o_prm.z[0][1] = +7.40; o_prm.z[1][0] = +7.70; o_prm.z[1][1] = +9.50;
    //o_prm.xi[0].x  = 0.25; o_prm.xi[0].y  = 0.25; o_prm.xi[1].x  = 0.75; o_prm.xi[1].y  = 0.75;
    //o_prm.eta[0].x = 0.25; o_prm.eta[0].y = 0.75; o_prm.eta[1].x = 0.75; o_prm.eta[1].y = 0.25;
    /***************************************************************************************************/

    // Regularization parameters
    OptimizeParameter r_prm;
    r_prm.k.resize(e_prm.Nc, e_prm.No, 0.0);
    r_prm.z.resize(e_prm.Nc, e_prm.No, 0.0);
    r_prm.xi.resize(e_prm.No);
    r_prm.eta.resize(e_prm.Nc);
    /***************************************************************************************************/
    r_prm.k[0][0]  =  -2.00;  r_prm.k[0][1]  =  +0.70;  r_prm.k[1][0]  =  +0.71;  r_prm.k[1][1]  =  -2.38;
    r_prm.z[0][0]  =  +7.96;  r_prm.z[0][1]  =  +5.83;  r_prm.z[1][0]  =  +7.68;  r_prm.z[1][1]  =  +9.72;
    r_prm.xi[0].x  = 0.0501;  r_prm.xi[0].y  = 0.0501;  r_prm.xi[1].x  = 0.9500;  r_prm.xi[1].y  = 0.0751;
    r_prm.eta[0].x = 0.4149;  r_prm.eta[0].y = 0.7549;  r_prm.eta[1].x = 0.4052;  r_prm.eta[1].y = 0.7077;
    /***************************************************************************************************/

    // Grid parameters
    double hx = 0.010; int Nx = 100;
    double hy = 0.010; int Ny = 100;
    double ht = 0.005; int Nt = 200;

    Dimension time(ht, 0, Nt);
    Dimension dimx(hx, 0, Nx);
    Dimension dimy(hy, 0, Ny);

    // Penalty paramteres
    DoubleVector r; r << 1.0000 << 10.000 << 50.0000 << 100.00;

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
        prob.vmin.resize(e_prm.Nc, +0.00);
        prob.vmax.resize(e_prm.Nc, +11.0);
        prob.U.resize(static_cast<const unsigned int>(Ny+1),
                      static_cast<const unsigned int>(Nx+1), 10.0);

        prob.regEpsilon = e[i];
        prob.r = r[i];
        if (i==0)
        {
            prob.PrmToVector(o_prm, x);
            prob.checkGradient1(prob, o_prm);
            IPrinter::printSeperatorLine();
        }

        //ConjugateGradient g;
        SteepestDescentGradient g;
        g.setFunction(&prob);
        g.setGradient(&prob);
        g.setPrinter(&prob);
        g.setProjection(&prob);
        g.setGradientNormalizer(&prob);
        g.setOptimalityTolerance(0.0001);
        g.setStepTolerance(0.0001);
        g.setFunctionTolerance(0.0001);
        g.setR1MinimizeEpsilon(0.1, 0.01);
        g.setNormalize(true);
        g.showExitMessage(true);

        g.calculate(x);

        IPrinter::printSeperatorLine(NULL, '=');
    }
}

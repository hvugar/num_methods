#include "exporter.h"
#include "problem2h.h"

Problem2HDirichlet prob;

void init_pr()
{
    // Equation parameters ---------------------------------------------------------------------
    EquationParameter e_prm;
    e_prm.a = 1.0;
    e_prm.lambda = 0.01;

    e_prm.Ns = 3;
    e_prm.q.resize(e_prm.Ns);
    e_prm.theta.resize(e_prm.Ns);

    e_prm.q[0] = 0.2; e_prm.theta[0].x = 0.5000; e_prm.theta[0].y = 0.5000;
    e_prm.q[1] = 0.3; e_prm.theta[1].x = 0.2000; e_prm.theta[1].y = 0.2000;
    e_prm.q[2] = 0.5; e_prm.theta[2].x = 0.8000; e_prm.theta[2].y = 0.8000;

    e_prm.No = 2;
    e_prm.Nc = 2;

    // Optimization parameters ---------------------------------------------------------------------
    OptimizeParameter o_prm;
    o_prm.k.resize(e_prm.Nc, e_prm.No, 0.0);
    o_prm.z.resize(e_prm.Nc, e_prm.No, 0.0);
    o_prm.xi.resize(e_prm.No);
    o_prm.eta.resize(e_prm.Nc);

    o_prm.k[0][0]  = +2.3400; o_prm.k[0][1]  = -2.7400; o_prm.k[1][0]  = +1.5800; o_prm.k[1][1]  = +1.9500;
    o_prm.z[0][0]  = +0.5000; o_prm.z[0][1]  = -0.4000; o_prm.z[1][0]  = -0.3000; o_prm.z[1][1]  = +0.6000;
    o_prm.xi[0].x  = +0.5500; o_prm.xi[0].y  = +0.1400; o_prm.xi[1].x  = +0.7400; o_prm.xi[1].y  = +0.3700;
    o_prm.eta[0].x = +0.2800; o_prm.eta[0].y = +0.7500; o_prm.eta[1].x = +0.8500; o_prm.eta[1].y = +0.8900;

    // Regulirization parameters ---------------------------------------------------------------------
    OptimizeParameter r_prm = o_prm;

    DoubleVector r; r << 1.0 << 2.0 << 10.0 << 100.0;

    double hx, hy; hx = hy = 0.01;
    unsigned Nx, Ny; Nx = Ny = 100;

    prob.setTimeDimension(Dimension(0.01, 0, 500));
    prob.addSpaceDimension(Dimension(hx, 0, Nx));
    prob.addSpaceDimension(Dimension(hy, 0, Ny));
    prob.mEquParameter = e_prm;
    prob.mOptParameter = o_prm;
    prob.mRegParameter = r_prm;
    prob.optimizeK = true;
    prob.optimizeZ = true;
    prob.optimizeC = true;
    prob.optimizeO = true;

    prob.alpha0 = 1.0; prob.V0.resize(Ny+1, Nx+1, 0.0);
    prob.alpha1 = 1.0; prob.V1.resize(Ny+1, Nx+1, 0.0);

    prob.regEpsilon = 0.0;

    prob.r = 1.0;
    prob.vmin.resize(e_prm.Nc, -2.0);
    prob.vmax.resize(e_prm.Nc, +2.0);
}

double call_fx(double *x)
{
    DoubleVector px(x, 16);
    return prob.fx(px);
}

void call_gr(double *x, double *g, unsigned int size)
{
    DoubleVector px(size); for (unsigned int i=0; i<size; i++) px[i] = x[i];
    DoubleVector gr(size); for (unsigned int i=0; i<size; i++) gr[i] = 0.0;

    prob.gradient(px, gr);

    for (unsigned int i=0; i<size; i++) g[i] = gr[i];

    px.clear();
    gr.clear();
}


void setPenaltyR(double r)
{
    prob.r = r;
}

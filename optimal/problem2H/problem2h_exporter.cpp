#include "problem2h_exporter.h"
#include "problem2h_solver.h"
#include "problem2h_solver1.h"

#define EXAMPLE4_SAMPLE_1

//Problem2HDirichlet prob;

void init_pr2()
{
//    // Equation parameters ---------------------------------------------------------------------
//    EquationParameterH e_prm;
//    e_prm.a = 1.0;
//    e_prm.lambda = 0.01;

//    // Pulse influences
//    e_prm.Ns = 2;
//    e_prm.pulses.resize(e_prm.Ns);
//    e_prm.pulses[0] = InitialPulse(SpacePoint(0.4300, 0.7500), -5.2);
//    e_prm.pulses[0] = InitialPulse(SpacePoint(0.4300, 0.7500));


//    e_prm.q[0] = -5.2; e_prm.theta[0].x = 0.4300; e_prm.theta[0].y = 0.7500;
//    e_prm.q[1] = -5.3; e_prm.theta[1].x = 0.8700; e_prm.theta[1].y = 0.2300;

//    e_prm.No = 2;
//    e_prm.Nc = 2;

//    // Optimization parameters
//    OptimizeParameterH o_prm;
//    o_prm.k.resize(e_prm.Nc, e_prm.No, 0.0);
//    o_prm.z.resize(e_prm.Nc, e_prm.No, 0.0);
//    o_prm.xi.resize(e_prm.No);
//    o_prm.eta.resize(e_prm.Nc);

//    // Regularization parameters
//    OptimizeParameterH r_prm;
//    r_prm.k.resize(e_prm.Nc, e_prm.No, 0.0);
//    r_prm.z.resize(e_prm.Nc, e_prm.No, 0.0);
//    r_prm.xi.resize(e_prm.No);
//    r_prm.eta.resize(e_prm.Nc);

//#ifdef EXAMPLE4_SAMPLE_1
//    o_prm.k[0][0]  = +1.1200; o_prm.k[0][1]  = +1.2400; o_prm.k[1][0]  = +1.4500; o_prm.k[1][1]  = +1.1800;
//    o_prm.z[0][0]  = +0.5000; o_prm.z[0][1]  = -0.4000; o_prm.z[1][0]  = +0.7000; o_prm.z[1][1]  = +0.5000;
//    o_prm.xi[0].x  = +0.4274; o_prm.xi[0].y  = +0.6735; o_prm.xi[1].x  = +0.6710; o_prm.xi[1].y  = +0.3851;
//    o_prm.eta[0].x = +0.5174; o_prm.eta[0].y = +0.7635; o_prm.eta[1].x = +0.5570; o_prm.eta[1].y = +0.4751;

//    r_prm.k[0][0]  = +0.4639; r_prm.k[0][1]  = -0.0136; r_prm.k[1][0]  = +0.1977; r_prm.k[1][1]  = -0.5896;
//    r_prm.z[0][0]  = +0.3014; r_prm.z[0][1]  = -0.6160; r_prm.z[1][0]  = -0.1914; r_prm.z[1][1]  = -0.2933;
//    r_prm.xi[0].x  = +0.4679; r_prm.xi[0].y  = +0.5770; r_prm.xi[1].x  = +0.7140; r_prm.xi[1].y  = +0.2614;
//    r_prm.eta[0].x = +0.5579; r_prm.eta[0].y = +0.8282; r_prm.eta[1].x = +0.8040; r_prm.eta[1].y = +0.7535;

//    //    r_prm.k[0][0]  = +0.5636; r_prm.k[0][1]  = -0.2421; r_prm.k[1][0]  = +0.2505; r_prm.k[1][1]  = -0.7679;
//    //    r_prm.z[0][0]  = +0.3220; r_prm.z[0][1]  = -0.6179; r_prm.z[1][0]  = -0.1833; r_prm.z[1][1]  = -0.3160;
//    //    r_prm.xi[0].x  = +0.4269; r_prm.xi[0].y  = +0.3364; r_prm.xi[1].x  = +0.8651; r_prm.xi[1].y  = +0.2336;
//    //    r_prm.eta[0].x = +0.5436; r_prm.eta[0].y = +0.5788; r_prm.eta[1].x = +0.5269; r_prm.eta[1].y = +0.5860;

//    //o_prm = r_prm;
//#endif

//#ifdef EXAMPLE4_SAMPLE_2
//    o_prm.k[0][0]  = -2.6400; o_prm.k[0][1]  = +3.7400; o_prm.k[1][0]  = -2.1800; o_prm.k[1][1]  = -2.0700;
//    o_prm.z[0][0]  = -0.9500; o_prm.z[0][1]  = +0.8500; o_prm.z[1][0]  = -0.1400; o_prm.z[1][1]  = -0.4500;
//    o_prm.xi[0].x  = +0.1486; o_prm.xi[0].y  = +0.1284; o_prm.xi[1].x  = +0.7525; o_prm.xi[1].y  = +0.7920;
//    o_prm.eta[0].x = +0.8512; o_prm.eta[0].y = +0.3245; o_prm.eta[1].x = +0.2854; o_prm.eta[1].y = +0.6515;

//    r_prm.k[0][0]  = -0.7053; r_prm.k[0][1]  = +0.6419; r_prm.k[1][0]  = -0.8886; r_prm.k[1][1]  = -1.2510;
//    r_prm.z[0][0]  = -1.9027; r_prm.z[0][1]  = +1.2513; r_prm.z[1][0]  = -0.1182; r_prm.z[1][1]  = -0.3907;
//    r_prm.xi[0].x  = +0.0500; r_prm.xi[0].y  = +0.0500; r_prm.xi[1].x  = +0.2210; r_prm.xi[1].y  = +0.8799;
//    r_prm.eta[0].x = +0.5281; r_prm.eta[0].y = +0.6057; r_prm.eta[1].x = +0.3210; r_prm.eta[1].y = +0.2266;
//    o_prm = r_prm;
//#endif

//    // Grid parameters
//    double hx = 0.010; int Nx = 100;
//    double hy = 0.010; int Ny = 100;
//    double ht = 0.010; int Nt = 500;

//    Dimension time(ht, 0, Nt);
//    Dimension dimx(hx, 0, Nx);
//    Dimension dimy(hy, 0, Ny);

//    prob.setTimeDimension(time);
//    prob.addSpaceDimension(dimx);
//    prob.addSpaceDimension(dimy);
//    prob.mEquParameter = e_prm;
//    prob.mOptParameter = o_prm;
//    prob.mRegParameter = r_prm;
//    prob.optimizeK = true;
//    prob.optimizeZ = true;
//    prob.optimizeC = true;
//    prob.optimizeO = true;
//    prob.LD = 50;

//    prob.regEpsilon = 0.0;
//    prob.r = 0.0;

//    prob.vmin.resize(e_prm.Nc, -1.5);
//    prob.vmax.resize(e_prm.Nc, +1.5);
}

void init_pr()
{
//    // Equation parameters
//    EquationParameterH e_prm;
//    e_prm.a = 1.0;
//    e_prm.lambda = +0.00;

//    // Pulse influences
//    e_prm.Ns = 2;
//    e_prm.q.resize(e_prm.Ns);
//    e_prm.theta.resize(e_prm.Ns);

//    e_prm.q[0] = +0.245; e_prm.theta[0].x = 0.3500; e_prm.theta[0].y = 0.3500;
//    e_prm.q[1] = +0.257; e_prm.theta[1].x = 0.6500; e_prm.theta[1].y = 0.6500;

//    e_prm.No = 2;
//    e_prm.Nc = 2;

//    // Optimization parameters
//    OptimizeParameterH o_prm;
//    o_prm.k.resize(e_prm.Nc, e_prm.No, 0.0);
//    o_prm.z.resize(e_prm.Nc, e_prm.No, 0.0);
//    o_prm.xi.resize(e_prm.No);
//    o_prm.eta.resize(e_prm.Nc);

//    o_prm.k[0][0]  = +1.1200; o_prm.k[0][1]  = +1.2400; o_prm.k[1][0]  = +1.4500; o_prm.k[1][1]  = +1.1800;
//    o_prm.z[0][0]  = +0.5000; o_prm.z[0][1]  = -0.4000; o_prm.z[1][0]  = +0.7000; o_prm.z[1][1]  = +0.5000;
//    o_prm.xi[0].x  = +0.4274; o_prm.xi[0].y  = +0.6735; o_prm.xi[1].x  = +0.6710; o_prm.xi[1].y  = +0.3851;
//    o_prm.eta[0].x = +0.5174; o_prm.eta[0].y = +0.7635; o_prm.eta[1].x = +0.5570; o_prm.eta[1].y = +0.4751;

//    // Regularization parameters
//    OptimizeParameterH r_prm;
//    r_prm.k.resize(e_prm.Nc, e_prm.No, 0.0);
//    r_prm.z.resize(e_prm.Nc, e_prm.No, 0.0);
//    r_prm.xi.resize(e_prm.No);
//    r_prm.eta.resize(e_prm.Nc);

//    r_prm.k[0][0]  = -0.0345; r_prm.k[0][1]  = -0.0541; r_prm.k[1][0]  = -0.0431; r_prm.k[1][1]  = -0.0984;
//    r_prm.z[0][0]  = +0.1245; r_prm.z[0][1]  = +0.0325; r_prm.z[1][0]  = -0.1452; r_prm.z[1][1]  = +0.2154;
//    r_prm.xi[0].x  = +0.2754; r_prm.xi[0].y  = +0.2438; r_prm.xi[1].x  = +0.4271; r_prm.xi[1].y  = +0.3824;
//    r_prm.eta[0].x = +0.3865; r_prm.eta[0].y = +0.4365; r_prm.eta[1].x = +0.6294; r_prm.eta[1].y = +0.5466;

//    r_prm.k[0][0]  = -0.3839; r_prm.k[0][1]  = -0.7371; r_prm.k[1][0]  = -0.5532; r_prm.k[1][1]  = -1.3051;
//    r_prm.z[0][0]  = +0.0238; r_prm.z[0][1]  = +0.0467; r_prm.z[1][0]  = -0.0352; r_prm.z[1][1]  = +0.0366;
//    r_prm.xi[0].x  = +0.5203; r_prm.xi[0].y  = +0.2636; r_prm.xi[1].x  = +0.2521; r_prm.xi[1].y  = +0.7225;
//    r_prm.eta[0].x = +0.5911; r_prm.eta[0].y = +0.4344; r_prm.eta[1].x = +0.4637; r_prm.eta[1].y = +0.5160;

//    r_prm.k[0][0]  = -0.6010; r_prm.k[0][1]  = -0.8657; r_prm.k[1][0]  = -0.8284; r_prm.k[1][1]  = -1.4035;
//    r_prm.z[0][0]  = -0.0108; r_prm.z[0][1]  = -0.0128; r_prm.z[1][0]  = -0.0671; r_prm.z[1][1]  = -0.0328;
//    r_prm.xi[0].x  = +0.5117; r_prm.xi[0].y  = +0.2698; r_prm.xi[1].x  = +0.2562; r_prm.xi[1].y  = +0.7156;
//    r_prm.eta[0].x = +0.5353; r_prm.eta[0].y = +0.4041; r_prm.eta[1].x = +0.4411; r_prm.eta[1].y = +0.5339;

//    //k: -0.3839  -0.7371  -0.5532  -1.3051
//    //z:  0.0238   0.0467  -0.0352   0.0366
//    //o:  0.5203   0.2636   0.2521   0.7225
//    //c:  0.5911   0.4344   0.4637   0.5160

//    //r_prm = o_prm;
//    //o_prm = r_prm;

//    // Grid parameters
//    double hx = 0.010; int Nx = 100;
//    double hy = 0.010; int Ny = 100;
//    double ht = 0.010; int Nt = 500;

//    Dimension time(ht, 0, Nt);
//    Dimension dimx(hx, 0, Nx);
//    Dimension dimy(hy, 0, Ny);

//    prob.setTimeDimension(time);
//    prob.addSpaceDimension(dimx);
//    prob.addSpaceDimension(dimy);
//    prob.mEquParameter = e_prm;
//    prob.mOptParameter = o_prm;
//    prob.mRegParameter = r_prm;
//    prob.optimizeK = true;
//    prob.optimizeZ = true;
//    prob.optimizeO = true;
//    prob.optimizeC = true;
//    prob.vmin.resize(e_prm.Nc, -0.05);
//    prob.vmax.resize(e_prm.Nc, +0.05);
//    prob.LD = 50;

//    prob.regEpsilon = 0.0;
//    prob.r = 0.0;
}

double call_fx(double *x)
{
//    DoubleVector px(x, 16);
//    return prob.fx(px);
    return 0.0;
}

void call_gr(double *x, double *g, unsigned int size)
{
//    DoubleVector px(size); for (unsigned int i=0; i<size; i++) px[i] = x[i];
//    DoubleVector gr(size); for (unsigned int i=0; i<size; i++) gr[i] = 0.0;

//    prob.gradient(px, gr);

//    for (unsigned int i=0; i<size; i++) g[i] = gr[i];

//    px.clear();
//    gr.clear();
}

void setPenaltyR(double r)
{
//    prob.r = r;
}

#include "problem2h_example.h"

void example1()
{
    // Equation parameters ---------------------------------------------------------------------
    EquationParameter e_prm;
    e_prm.a = 1.0;
    e_prm.lambda = 0.01;

    e_prm.Ns = 5;
    e_prm.q.resize(e_prm.Ns);
    e_prm.theta.resize(e_prm.Ns);

    e_prm.q[0] = +10.2; e_prm.theta[0].x = 0.5000; e_prm.theta[0].y = 0.5000;
    e_prm.q[1] = +10.3; e_prm.theta[1].x = 0.2000; e_prm.theta[1].y = 0.2000;
    e_prm.q[2] = +10.5; e_prm.theta[2].x = 0.8000; e_prm.theta[2].y = 0.8000;
    e_prm.q[3] = +10.8; e_prm.theta[3].x = 0.3000; e_prm.theta[3].y = 0.7000;
    e_prm.q[4] = +10.4; e_prm.theta[4].x = 0.8000; e_prm.theta[4].y = 0.3000;

    e_prm.No = 2;
    e_prm.Nc = 2;

    // Optimization parameters ---------------------------------------------------------------------
    OptimizeParameter o_prm;
    o_prm.k.resize(e_prm.Nc, e_prm.No, 0.0);
    o_prm.z.resize(e_prm.Nc, e_prm.No, 0.0);
    o_prm.xi.resize(e_prm.No);
    o_prm.eta.resize(e_prm.Nc);

    //o_prm.k[0][0]  = +1.1200; o_prm.k[0][1]  = +1.2400; o_prm.k[1][0]  = +2.4500; o_prm.k[1][1]  = +2.1800;
    //o_prm.z[0][0]  = +1.5000; o_prm.z[0][1]  = -1.4000; o_prm.z[1][0]  = +1.7000; o_prm.z[1][1]  = +1.5000;
    //o_prm.xi[0].x  = +0.3000; o_prm.xi[0].y  = +0.8000; o_prm.xi[1].x  = +0.6000; o_prm.xi[1].y  = +0.4000;
    //o_prm.eta[0].x = +0.5000; o_prm.eta[0].y = +0.7000; o_prm.eta[1].x = +0.7000; o_prm.eta[1].y = +0.3000;

    //o_prm.k[0][0]  = +1.0000; o_prm.k[0][1]  = +1.0000; o_prm.k[1][0]  = +1.0000; o_prm.k[1][1]  = +1.0000;
    //o_prm.z[0][0]  = +0.5000; o_prm.z[0][1]  = +0.7000; o_prm.z[1][0]  = +0.4000; o_prm.z[1][1]  = +8.0000;
    o_prm.k[0][0]  = +1.1200; o_prm.k[0][1]  = +1.2400; o_prm.k[1][0]  = +1.4500; o_prm.k[1][1]  = +1.1800;
    o_prm.z[0][0]  = +0.5000; o_prm.z[0][1]  = -0.4000; o_prm.z[1][0]  = +0.7000; o_prm.z[1][1]  = +0.5000;
    //o_prm.xi[0].x  = +0.3000; o_prm.xi[0].y  = +0.8000; o_prm.xi[1].x  = +0.6000; o_prm.xi[1].y  = +0.4000;
    //o_prm.eta[0].x = +0.1000; o_prm.eta[0].y = +0.6000; o_prm.eta[1].x = +0.9000; o_prm.eta[1].y = +0.2000;
    o_prm.xi[0].x  = +0.4274; o_prm.xi[0].y  = +0.6735; o_prm.xi[1].x  = +0.6710; o_prm.xi[1].y  = +0.3851;
    o_prm.eta[0].x = +0.5174; o_prm.eta[0].y = +0.7635; o_prm.eta[1].x = +0.5570; o_prm.eta[1].y = +0.4751;

    //o_prm.xi[0].x  = +0.4074; o_prm.xi[0].y  = +0.6735; o_prm.xi[1].x  = +0.6710; o_prm.xi[1].y  = +0.2851;
    //o_prm.eta[0].x = +0.5000; o_prm.eta[0].y = +0.7635; o_prm.eta[1].x = +0.5000; o_prm.eta[1].y = +0.4751;

    // Regulirization parameters ---------------------------------------------------------------------
    OptimizeParameter r_prm;
    r_prm.k.resize(e_prm.Nc, e_prm.No, 0.0);
    r_prm.z.resize(e_prm.Nc, e_prm.No, 0.0);
    r_prm.xi.resize(e_prm.No);
    r_prm.eta.resize(e_prm.Nc);

    // r = 0.01
    r_prm.k[0][0]  = +0.1433; r_prm.k[0][1]  = +0.2481; r_prm.k[1][0]  = +1.0715; r_prm.k[1][1]  = +0.7702;
    r_prm.z[0][0]  = +1.1141; r_prm.z[0][1]  = -1.7490; r_prm.z[1][0]  = +0.6137; r_prm.z[1][1]  = +0.7212;
    r_prm.xi[0].x  = +0.1137; r_prm.xi[0].y  = +0.8865; r_prm.xi[1].x  = +0.0502; r_prm.xi[1].y  = +0.9498;
    r_prm.eta[0].x = +0.9420; r_prm.eta[0].y = +0.7940; r_prm.eta[1].x = +0.9472; r_prm.eta[1].y = +0.0528;

    // r = 100.0
    //r_prm.k[0][0]  = -0.0953; r_prm.k[0][1]  = +0.1448; r_prm.k[1][0]  = +1.0591; r_prm.k[1][1]  = -0.7964;
    //r_prm.z[0][0]  = +1.1060; r_prm.z[0][1]  = -1.7351; r_prm.z[1][0]  = +0.5721; r_prm.z[1][1]  = +0.7514;
    //r_prm.xi[0].x  = +0.2024; r_prm.xi[0].y  = +0.7993; r_prm.xi[1].x  = +0.0513; r_prm.xi[1].y  = +0.9487;
    //r_prm.eta[0].x = +0.9410; r_prm.eta[0].y = +0.7091; r_prm.eta[1].x = +0.8438; r_prm.eta[1].y = +0.1553;

    //r_prm = o_prm;
    // Regulirization parameters ---------------------------------------------------------------------

    DoubleVector r; r << 0.0000 << 2.0000 << 5.0000 << 10.000 << 100.00;
    DoubleVector e; e << 0.0000 << 0.0000 << 0.0000 << 0.0000 << 0.0000;

    double hx, hy; hx = hy = 0.01;
    int Nx, Ny; Nx = Ny = 100;
    double ht = 0.01;
    int Nt = 500;

    Dimension time(ht, 0, Nt);
    Dimension dimx(hx, 0, Nx);
    Dimension dimy(hy, 0, Ny);

    // Checking gradients ----------------------------------------------------------------------------

    DoubleVector x;
    for (unsigned int i=0; i<r.length(); i++)
    {
        Problem2HNDirichlet prob;
        prob.setTimeDimension(time);
        prob.addSpaceDimension(dimx);
        prob.addSpaceDimension(dimy);
        prob.mEquParameter = e_prm;
        prob.mOptParameter = o_prm;
        prob.mRegParameter = r_prm;
        prob.optimizeK = true;
        prob.optimizeZ = true;
        prob.optimizeC = true;
        prob.optimizeO = true;

        prob.V0.resize(Ny+1, Nx+1, 0.0);

        prob.regEpsilon = e[i];
        prob.r = r[i];
        prob.vmin.resize(e_prm.Nc, -5.0);
        prob.vmax.resize(e_prm.Nc, +5.0);
        prob.LD = 50;

        if (i==0)
        {
            prob.PrmToVector(o_prm, x);

            //checkGradient(prob);
            //IPrinter::printSeperatorLine();
        }

        //ConjugateGradient g;
        SteepestDescentGradient g;
        g.setFunction(&prob);
        g.setGradient(&prob);
        g.setPrinter(&prob);
        g.setProjection(&prob);
        g.setEpsilon1(0.00001);
        g.setEpsilon2(0.00001);
        g.setEpsilon3(0.00001);
        g.setR1MinimizeEpsilon(1.0, 0.001);
        g.setNormalize(true);
        g.showExitMessage(true);

        g.calculate(x);

        IPrinter::printSeperatorLine(NULL, '=');
    }
}

void example2()
{
    // Equation parameters ---------------------------------------------------------------------
    EquationParameter e_prm;
    e_prm.a = 1.0;
    e_prm.lambda = 0.01;

    e_prm.Ns = 3;
    e_prm.q.resize(e_prm.Ns);
    e_prm.theta.resize(e_prm.Ns);

    //    e_prm.q[0] = 0.8; e_prm.theta[0].x = 0.2500; e_prm.theta[0].y = 0.7500;
    //    e_prm.q[1] = 0.1; e_prm.theta[1].x = 0.3500; e_prm.theta[1].y = 0.9500;
    //    e_prm.q[2] = 0.5; e_prm.theta[2].x = 0.8500; e_prm.theta[2].y = 0.4500;

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

    o_prm.k.resize(e_prm.Nc, e_prm.No, 0.0);
    o_prm.z.resize(e_prm.Nc, e_prm.No, 0.0);
    o_prm.xi.resize(e_prm.No);
    o_prm.eta.resize(e_prm.Nc);

    o_prm.k[0][0]  = +2.3400; o_prm.k[0][1]  = -2.7400; o_prm.k[1][0]  = +1.5800; o_prm.k[1][1]  = +1.9500;
    o_prm.z[0][0]  = +0.5000; o_prm.z[0][1]  = -0.4000; o_prm.z[1][0]  = -0.3000; o_prm.z[1][1]  = +0.6000;
    o_prm.xi[0].x  = +0.5500; o_prm.xi[0].y  = +0.1400; o_prm.xi[1].x  = +0.7400; o_prm.xi[1].y  = +0.3700;
    o_prm.eta[0].x = +0.2800; o_prm.eta[0].y = +0.7500; o_prm.eta[1].x = +0.8500; o_prm.eta[1].y = +0.8900;
    //    +2.3400, -2.7400, +1.5800, +1.9500, +0.5000, -0.4000, -0.3000, +0.6000,
    //    +0.5500, +0.1400, +0.7400, +0.3700, +0.2800, +0.7500, +0.8500, +0.8900,

    //    o_prm.k[0][0] = +2.4229; o_prm.k[0][1] = -2.4453;
    //    o_prm.k[1][0] = +1.8412; o_prm.k[1][1] = +2.1660;
    //    o_prm.k[0][0] = +2.0440; o_prm.k[0][1] = -2.3185;
    //    o_prm.k[1][0] = +1.8133; o_prm.k[1][1] = +2.3553;

    //    o_prm.z[0][0] = +0.2000; o_prm.z[0][1] = -0.0439;
    //    o_prm.z[1][0] = -0.3203; o_prm.z[1][1] = +0.5721;
    //    o_prm.z[0][0] = +0.2198; o_prm.z[0][1] = -0.0629;
    //    o_prm.z[1][0] = -0.3247; o_prm.z[1][1] = +0.5672;

    //    o_prm.xi[0].x = 0.5030; o_prm.xi[0].y = 0.2125;
    //    o_prm.xi[1].x = 0.8903; o_prm.xi[1].y = 0.4200;
    //    o_prm.xi[0].x = +0.4951; o_prm.xi[0].y = +0.2267;
    //    o_prm.xi[1].x = +0.9068; o_prm.xi[1].y = +0.4109;

    //    o_prm.eta[0].x = 0.1804; o_prm.eta[0].y = 0.7732;
    //    o_prm.eta[1].x = 0.8029; o_prm.eta[1].y = 0.6982;
    //    o_prm.eta[0].x = +0.1745; o_prm.eta[0].y = +0.7807;
    //    o_prm.eta[1].x = +0.7901; o_prm.eta[1].y = +0.7071;


    //    +2.0440,-2.3185,+1.8133,+2.3553,
    //    +0.2198,-0.0629,-0.3247,+0.5672,
    //    +0.4951,+0.2267,+0.9068,+0.4109,
    //    +0.1745,+0.7807,+0.7901,+0.7071;

    //    o_prm.k[0][0]  = +2.3199; o_prm.k[0][1]  = -2.6193; o_prm.k[1][0]  = +1.5956; o_prm.k[1][1]  = +1.9759;
    //    o_prm.z[0][0]  = +0.2105; o_prm.z[0][1]  = -0.0597; o_prm.z[1][0]  = -0.2872; o_prm.z[1][1]  = +0.6144;
    //    o_prm.xi[0].x  = +0.3810; o_prm.xi[0].y  = +0.2625; o_prm.xi[1].x  = +0.8788; o_prm.xi[1].y  = +0.4078;
    //    o_prm.eta[0].x = +0.2957; o_prm.eta[0].y = +0.6892; o_prm.eta[1].x = +0.7464; o_prm.eta[1].y = +0.7672;
    //    +2.3199, -2.6193, +1.5956, +1.9759, +0.2105  -0.0597, -0.2872, +0.6144,
    //    +0.3810  +0.2625, +0.8788, +0.4078, +0.2957, +0.6892, +0.7464, +0.7672; 0.0037;

    //+2.4045,  -2.4473,  +1.9092,   +2.1860,   +0.2031,  -0.0511,  -0.3136,   +0.5790,
    //+0.5043,  +0.2105,  +0.8899,   +0.4208,   +0.1791,  +0.7716,  +0.8049,   +0.6970, 0.004803;


    // Regulirization parameters ---------------------------------------------------------------------
    OptimizeParameter r_prm = o_prm;

    DoubleVector r; r << 1.0 << 2.0 << 10.0 << 100.0;

    double hx, hy; hx = hy = 0.01;
    unsigned Nx, Ny; Nx = Ny = 100;

    Dimension time(0.01, 0, 500);
    Dimension dimx(hx, 0, Nx);
    Dimension dimy(hy, 0, Ny);

    // Checking gradients ---------------------------------------------------------------------

    //checkGradient(prob);
    //IPrinter::printSeperatorLine();

    // ----------------------------------------------------------------------------------------

    DoubleVector x;
    for (unsigned int i=0; i<r.length(); i++)
    {
        Problem2HNDirichlet prob;
        prob.setTimeDimension(time);
        prob.addSpaceDimension(dimx);
        prob.addSpaceDimension(dimy);
        prob.mEquParameter = e_prm;
        prob.mOptParameter = o_prm;
        prob.mRegParameter = r_prm;
        prob.optimizeK = true;
        prob.optimizeZ = true;
        prob.optimizeC = true;
        prob.optimizeO = true;

        prob.V0.resize(Ny+1, Nx+1, 0.0);

        prob.regEpsilon = 0.0;

        prob.r = r[i];
        prob.vmin.resize(e_prm.Nc, -2.0);
        prob.vmax.resize(e_prm.Nc, +2.0);
        prob.LD = 5;

        if (i==0)
        {
            prob.PrmToVector(o_prm, x);

            //checkGradient(prob);
            //IPrinter::printSeperatorLine();
        }

        ConjugateGradient g;
        g.setFunction(&prob);
        g.setGradient(&prob);
        g.setPrinter(&prob);
        g.setProjection(&prob);
        g.setEpsilon1(0.01);
        g.setEpsilon2(0.01);
        g.setEpsilon3(0.01);
        g.setR1MinimizeEpsilon(0.1, 0.001);
        g.setNormalize(true);
        g.showExitMessage(true);
        g.setResetIteration(true);

        g.calculate(x);

        IPrinter::printSeperatorLine(NULL, '=');
    }
}

void example3()
{
    Problem2HNDirichlet prob;
    unsigned int Nx = 100;
    unsigned int Ny = 100;
    unsigned int Nt = 100;
    double hx = 0.01;
    double hy = 0.01;
    double ht = 0.01;

    prob.setTimeDimension(Dimension(ht, 0, Nt));
    prob.addSpaceDimension(Dimension(hx, 0, Nx));
    prob.addSpaceDimension(Dimension(hy, 0, Ny));

    EquationParameter eprm;
    OptimizeParameter oprm;
    //OptimizeParameter oprm0;
    //prob.initParameters(eprm, oprm, oprm0);

    eprm.a = +1.0;
    eprm.lambda = +0.01;
    eprm.Ns = 1;
    eprm.q.resize(eprm.Ns);
    eprm.theta.resize(eprm.Ns);

    eprm.q[0] = 1.0; eprm.theta[0].x = 0.5000; eprm.theta[0].y = 0.5000;
    //eprm.q[1] = 1.0; eprm.theta[1].x = 0.2500; eprm.theta[1].y = 0.2500;
    //eprm.q[2] = 1.0; eprm.theta[2].x = 0.7500; eprm.theta[2].y = 0.7500;
    //eprm.q[3] = 1.0; eprm.theta[3].x = 0.2500; eprm.theta[3].y = 0.7500;
    //eprm.q[4] = 1.0; eprm.theta[4].x = 0.7500; eprm.theta[4].y = 0.2500;

    eprm.No = 2;
    eprm.Nc = 2;

    oprm.k.resize(eprm.Nc, eprm.No, 0.0);
    oprm.z.resize(eprm.Nc, eprm.No, 0.0);
    oprm.xi.resize(eprm.No);
    oprm.eta.resize(eprm.Nc);

    oprm.k[0][0]  = -1.0000; oprm.k[0][1]  = -1.0000; oprm.k[1][0]  = -1.0000; oprm.k[1][1]  = -1.0000;
    oprm.z[0][0]  = +1.0000; oprm.z[0][1]  = +1.0000; oprm.z[1][0]  = +1.0000; oprm.z[1][1]  = +1.0000;
    oprm.xi[0].x  = +0.2500; oprm.xi[0].y  = +0.2500; oprm.xi[1].x  = +0.7500; oprm.xi[1].y  = +0.7500;
    oprm.eta[0].x = +0.2500; oprm.eta[0].y = +0.2500; oprm.eta[1].x = +0.7500; oprm.eta[1].y = +0.7500;

    prob.mEquParameter = eprm;
    prob.mOptParameter = oprm;
    prob.mRegParameter = oprm;
    prob.regEpsilon = 0.0;
    prob.r = 0.0;
    prob.LD = 1;

    prob.optimizeK = true;
    prob.optimizeZ = true;
    prob.optimizeC = true;
    prob.optimizeO = true;
    prob.V0.resize(Ny+1, Nx+1, 0.0);

    std::vector<DoubleMatrix> u;
    spif_vector info;
    prob.solveForwardIBVP(u, info, false);

    //    DoubleVector x;
    //    prob.PrmToVector(oprm, x);
    //    double res = prob.fx(x);
    //    printf("%f\n", res);

    //    printf("%f\n", prob.integralU(u.at(0)));
}

void example4()
{
    // Equation parameters
    EquationParameter e_prm;
    e_prm.a = 1.0;
    e_prm.lambda = +0.01;

    // Pulse influences
    e_prm.Ns = 2;
    e_prm.q.resize(e_prm.Ns);
    e_prm.theta.resize(e_prm.Ns);

    e_prm.q[0] = -5.2; e_prm.theta[0].x = 0.4300; e_prm.theta[0].y = 0.7500;
    e_prm.q[1] = -5.3; e_prm.theta[1].x = 0.8700; e_prm.theta[1].y = 0.2300;

    e_prm.No = 2;
    e_prm.Nc = 2;

    // Optimization parameters
    OptimizeParameter o_prm;

    o_prm.k.resize(e_prm.Nc, e_prm.No, 0.0);
    o_prm.z.resize(e_prm.Nc, e_prm.No, 0.0);
    o_prm.xi.resize(e_prm.No);
    o_prm.eta.resize(e_prm.Nc);

    // Regularization parameters
    OptimizeParameter r_prm;
    r_prm.k.resize(e_prm.Nc, e_prm.No, 0.0);
    r_prm.z.resize(e_prm.Nc, e_prm.No, 0.0);
    r_prm.xi.resize(e_prm.No);
    r_prm.eta.resize(e_prm.Nc);

#ifdef EXAMPLE4_SAMPLE_1
    o_prm.k[0][0]  = +1.1200; o_prm.k[0][1]  = +1.2400; o_prm.k[1][0]  = +1.4500; o_prm.k[1][1]  = +1.1800;
    o_prm.z[0][0]  = +0.5000; o_prm.z[0][1]  = -0.4000; o_prm.z[1][0]  = +0.7000; o_prm.z[1][1]  = +0.5000;
    o_prm.xi[0].x  = +0.4274; o_prm.xi[0].y  = +0.6735; o_prm.xi[1].x  = +0.6710; o_prm.xi[1].y  = +0.3851;
    o_prm.eta[0].x = +0.5174; o_prm.eta[0].y = +0.7635; o_prm.eta[1].x = +0.5570; o_prm.eta[1].y = +0.4751;

    //r_prm.k[0][0]  = +0.4639; r_prm.k[0][1]  = -0.0136; r_prm.k[1][0]  = +0.1977; r_prm.k[1][1]  = -0.5896;
    //r_prm.z[0][0]  = +0.3014; r_prm.z[0][1]  = -0.6160; r_prm.z[1][0]  = -0.1914; r_prm.z[1][1]  = -0.2933;
    //r_prm.xi[0].x  = +0.4679; r_prm.xi[0].y  = +0.5770; r_prm.xi[1].x  = +0.7140; r_prm.xi[1].y  = +0.2614;
    //r_prm.eta[0].x = +0.5579; r_prm.eta[0].y = +0.8282; r_prm.eta[1].x = +0.8040; r_prm.eta[1].y = +0.7535;

    r_prm.k[0][0]  = +0.5636; r_prm.k[0][1]  = -0.2421; r_prm.k[1][0]  = +0.2505; r_prm.k[1][1]  = -0.7679;
    r_prm.z[0][0]  = +0.3220; r_prm.z[0][1]  = -0.6179; r_prm.z[1][0]  = -0.1833; r_prm.z[1][1]  = -0.3160;
    r_prm.xi[0].x  = +0.4269; r_prm.xi[0].y  = +0.3364; r_prm.xi[1].x  = +0.8651; r_prm.xi[1].y  = +0.2336;
    r_prm.eta[0].x = +0.5436; r_prm.eta[0].y = +0.5788; r_prm.eta[1].x = +0.5269; r_prm.eta[1].y = +0.5860;

    r_prm = o_prm;

    //Measure k:-1.1212 -1.2418 -1.4580 -1.1793 z: 0.4961 -0.4043  0.6879  0.4902 o: 0.3976 0.8678 0.8008 0.3278 c: 0.5370 0.7321 0.5276 0.4524 // 0.060
    //Control k:-1.1227 -1.2402 -1.4438 -1.1749 z: 0.4976 -0.4026  0.6870  0.4894 o: 0.4336 0.7616 0.6664 0.4268 c: 0.5336 0.8616 0.5664 0.5268 // 0.046
#endif

#ifdef EXAMPLE4_SAMPLE_2
    o_prm.k[0][0]  = -2.6400; o_prm.k[0][1]  = +3.7400; o_prm.k[1][0]  = -2.1800; o_prm.k[1][1]  = -2.0700;
    o_prm.z[0][0]  = -0.9500; o_prm.z[0][1]  = +0.8500; o_prm.z[1][0]  = -0.1400; o_prm.z[1][1]  = -0.4500;
    o_prm.xi[0].x  = +0.1486; o_prm.xi[0].y  = +0.1284; o_prm.xi[1].x  = +0.7525; o_prm.xi[1].y  = +0.7920;
    o_prm.eta[0].x = +0.8512; o_prm.eta[0].y = +0.3245; o_prm.eta[1].x = +0.2854; o_prm.eta[1].y = +0.6515;

    r_prm.k[0][0]  = -0.7053; r_prm.k[0][1]  = +0.6419; r_prm.k[1][0]  = -0.8886; r_prm.k[1][1]  = -1.2510;
    r_prm.z[0][0]  = -1.9027; r_prm.z[0][1]  = +1.2513; r_prm.z[1][0]  = -0.1182; r_prm.z[1][1]  = -0.3907;
    r_prm.xi[0].x  = +0.0500; r_prm.xi[0].y  = +0.0500; r_prm.xi[1].x  = +0.2210; r_prm.xi[1].y  = +0.8799;
    r_prm.eta[0].x = +0.5281; r_prm.eta[0].y = +0.6057; r_prm.eta[1].x = +0.3210; r_prm.eta[1].y = +0.2266;
    //o_prm = r_prm;
#endif

    // Grid parameters
    double hx = 0.010; int Nx = 100;
    double hy = 0.010; int Ny = 100;
    double ht = 0.010; int Nt = 500;

    Dimension time(ht, 0, Nt);
    Dimension dimx(hx, 0, Nx);
    Dimension dimy(hy, 0, Ny);

    // Penalty paramteres
    DoubleVector r; r << 0.0000;// << 10.000 << 100.000 << 1000.0;// << 20.000 << 50.000 << 100.00;
    // Regularization coefficients
    DoubleVector e; e << 0.0000;// << 0.0000 << 0.00000 << 0.0000;// << 0.0000 << 0.0000 << 0.0000;

    DoubleVector x;
    for (unsigned int i=0; i<r.length(); i++)
    {
        Problem2HNDirichlet prob;
        prob.setTimeDimension(time);
        prob.addSpaceDimension(dimx);
        prob.addSpaceDimension(dimy);
        prob.mEquParameter = e_prm;
        prob.mOptParameter = o_prm;
        prob.mRegParameter = r_prm;
        prob.optimizeK = true;
        prob.optimizeZ = true;
        prob.optimizeC = true;
        prob.optimizeO = true;

        prob.V0.resize(Ny+1, Nx+1, 0.0);

        prob.regEpsilon = e[i];
        prob.r = r[i];

        prob.vmin.resize(e_prm.Nc, -1.5);
        prob.vmax.resize(e_prm.Nc, +1.5);
        prob.LD = 50;

        if (i==0)
        {
            prob.PrmToVector(o_prm, x);
            //prob.checkGradient(prob);
            //IPrinter::printSeperatorLine();
        }

        //ConjugateGradient g;
        SteepestDescentGradient g;
        //ConstStepGradient g;
        g.setFunction(&prob);
        g.setGradient(&prob);
        g.setPrinter(&prob);
        g.setProjection(&prob);
        g.setEpsilon1(0.00001);
        g.setEpsilon2(0.00001);
        g.setEpsilon3(0.00001);
        g.setR1MinimizeEpsilon(0.1, 0.01);
        g.setNormalize(true);
        g.showExitMessage(true);
        prob.gm = &g;

//        DoubleVector gr;
//        prob.gradient(x, gr);
//        gr[0] = gr[1] = gr[2] = gr[3] = gr[4] = gr[5] = gr[6] = gr[7] = 0.0;
//        gr.L2Normalize();
//        g.mg = &gr;
//        g.mx = &x;
//        DoubleVector cx(gr.length());
//        IPrinter::print(gr.mid(8, 15),10,6,4);
//        for (int i = 0; i <= 200;  i++)
//        {
//            double a = 0.01*i;
//            for (unsigned int i=0; i<gr.length(); i++)
//            {
//                cx[i] = x[i] - a * gr[i];
//                //if (m_projection != NULL) m_projection->project(cx, i);
//            }
//            prob.project(cx);
//            double f = prob.fx(cx);
//            printf("%.3f %f ", a, f);
//            printf("ko: %6.4f %6.4f %6.4f %6.4f c: %6.4f %6.4f %6.4f %6.4f\n", cx[8], cx[9], cx[10], cx[11], cx[12], cx[13], cx[14], cx[15]);

//        }
//        exit(-1);

        g.calculate(x);

        IPrinter::printSeperatorLine(NULL, '=');
    }
}

void example5()
{
    // Equation parameters
    EquationParameter e_prm;
    e_prm.a = 1.0;
    e_prm.lambda = 0.01;

    // Pulse influences
    e_prm.Ns = 2;
    e_prm.q.resize(e_prm.Ns);
    e_prm.theta.resize(e_prm.Ns);

    e_prm.q[0] = -5.2; e_prm.theta[0].x = 0.4300; e_prm.theta[0].y = 0.7500;
    e_prm.q[1] = -5.3; e_prm.theta[1].x = 0.8700; e_prm.theta[1].y = 0.2300;

    e_prm.No = 2;
    e_prm.Nc = 3;

    // Optimization parameters
    OptimizeParameter o_prm;

    o_prm.k.resize(e_prm.Nc, e_prm.No, 0.0);
    o_prm.z.resize(e_prm.Nc, e_prm.No, 0.0);
    o_prm.xi.resize(e_prm.No);
    o_prm.eta.resize(e_prm.Nc);

    // Regularization parameters
    OptimizeParameter r_prm;
    r_prm.k.resize(e_prm.Nc, e_prm.No, 0.0);
    r_prm.z.resize(e_prm.Nc, e_prm.No, 0.0);
    r_prm.xi.resize(e_prm.No);
    r_prm.eta.resize(e_prm.Nc);

#ifdef EXAMPLE4_SAMPLE_1
    o_prm.k[0][0]  = +1.1200; o_prm.k[0][1]  = +1.2400; o_prm.k[0][2]  = +1.2400; o_prm.k[1][0]  = +1.4500; o_prm.k[1][1]  = +1.1800; o_prm.k[1][2]  = +1.1800;
    o_prm.z[0][0]  = +0.5000; o_prm.z[0][1]  = -0.4000; o_prm.z[0][2]  = -0.8000; o_prm.z[1][0]  = +0.7000; o_prm.z[1][1]  = +0.5000; o_prm.z[1][2]  = +0.5000;
    o_prm.xi[0].x  = +0.4274; o_prm.xi[0].y  = +0.6735; o_prm.xi[1].x  = +0.6710; o_prm.xi[1].y  = +0.3851;
    o_prm.eta[0].x = +0.5174; o_prm.eta[0].y = +0.7635; o_prm.eta[1].x = +0.5570; o_prm.eta[1].y = +0.4751; o_prm.eta[2].x = +0.2570; o_prm.eta[2].y = +0.1751;

    //    r_prm.k[0][0]  = +0.4639; r_prm.k[0][1]  = -0.0136; r_prm.k[1][0]  = +0.1977; r_prm.k[1][1]  = -0.5896;
    //    r_prm.z[0][0]  = +0.3014; r_prm.z[0][1]  = -0.6160; r_prm.z[1][0]  = -0.1914; r_prm.z[1][1]  = -0.2933;
    //    r_prm.xi[0].x  = +0.4679; r_prm.xi[0].y  = +0.5770; r_prm.xi[1].x  = +0.7140; r_prm.xi[1].y  = +0.2614;
    //    r_prm.eta[0].x = +0.5579; r_prm.eta[0].y = +0.8282; r_prm.eta[1].x = +0.8040; r_prm.eta[1].y = +0.7535;

    r_prm.k[0][0]  = +0.5636; r_prm.k[0][1]  = -0.2421; r_prm. k[0][2]  = -0.2421; r_prm.k[1][0]  = +0.2505; r_prm.k[1][1]  = -0.7679; r_prm.k[1][2]  = -0.7679;
    r_prm.z[0][0]  = +0.3220; r_prm.z[0][1]  = -0.6179; r_prm.z[0][2]  = -0.6179; r_prm.z[1][0]  = -0.1833; r_prm.z[1][1]  = -0.3160; r_prm.z[1][2]  = -0.3160;
    r_prm.xi[0].x  = +0.4269; r_prm.xi[0].y  = +0.3364; r_prm.xi[1].x  = +0.8651; r_prm.xi[1].y  = +0.2336;
    r_prm.eta[0].x = +0.5436; r_prm.eta[0].y = +0.5788; r_prm.eta[1].x = +0.5269; r_prm.eta[1].y = +0.5860; r_prm.eta[2].x = +0.3269; r_prm.eta[2].y = +0.8860;

    //    o_prm = r_prm;
#endif

#ifdef EXAMPLE4_SAMPLE_2
    o_prm.k[0][0]  = -2.6400; o_prm.k[0][1]  = +3.7400; o_prm.k[1][0]  = -2.1800; o_prm.k[1][1]  = -2.0700;
    o_prm.z[0][0]  = -0.9500; o_prm.z[0][1]  = +0.8500; o_prm.z[1][0]  = -0.1400; o_prm.z[1][1]  = -0.4500;
    o_prm.xi[0].x  = +0.1486; o_prm.xi[0].y  = +0.1284; o_prm.xi[1].x  = +0.7525; o_prm.xi[1].y  = +0.7920;
    o_prm.eta[0].x = +0.8512; o_prm.eta[0].y = +0.3245; o_prm.eta[1].x = +0.2854; o_prm.eta[1].y = +0.6515;

    r_prm.k[0][0]  = -0.7053; r_prm.k[0][1]  = +0.6419; r_prm.k[1][0]  = -0.8886; r_prm.k[1][1]  = -1.2510;
    r_prm.z[0][0]  = -1.9027; r_prm.z[0][1]  = +1.2513; r_prm.z[1][0]  = -0.1182; r_prm.z[1][1]  = -0.3907;
    r_prm.xi[0].x  = +0.0500; r_prm.xi[0].y  = +0.0500; r_prm.xi[1].x  = +0.2210; r_prm.xi[1].y  = +0.8799;
    r_prm.eta[0].x = +0.5281; r_prm.eta[0].y = +0.6057; r_prm.eta[1].x = +0.3210; r_prm.eta[1].y = +0.2266;
    o_prm = r_prm;
#endif

    // Grid parameters
    double hx = 0.005; int Nx = 200;
    double hy = 0.005; int Ny = 200;
    double ht = 0.005; int Nt = 200;

    Dimension time(ht, 0, Nt);
    Dimension dimx(hx, 0, Nx);
    Dimension dimy(hy, 0, Ny);

    // Penalty paramteres
    DoubleVector r; r << 0.0000 << 10.000 << 100.000 << 0.0000;// << 20.000 << 50.000 << 100.00;
    // Regularization coefficients
    DoubleVector e; e << 0.0000 << 0.0000 << 0.00000 << 0.0000;// << 0.0000 << 0.0000 << 0.0000;

    DoubleVector x;
    for (unsigned int i=0; i<r.length(); i++)
    {
        Problem2HNDirichlet prob;
        prob.setTimeDimension(time);
        prob.addSpaceDimension(dimx);
        prob.addSpaceDimension(dimy);
        prob.mEquParameter = e_prm;
        prob.mOptParameter = o_prm;
        prob.mRegParameter = r_prm;
        prob.optimizeK = true;
        prob.optimizeZ = true;
        prob.optimizeC = true;
        prob.optimizeO = true;

        prob.V0.resize(Ny+1, Nx+1, 0.0);

        prob.regEpsilon = e[i];
        prob.r = r[i];

        prob.vmin.resize(e_prm.Nc, -1.5);
        prob.vmax.resize(e_prm.Nc, +1.5);
        prob.LD = 50;

        if (i==0)
        {
            prob.PrmToVector(o_prm, x);
            //prob.checkGradient(prob);
            //IPrinter::printSeperatorLine();
        }

        //        ConjugateGradient g;
        //        g.setResetIteration(true);
        SteepestDescentGradient g;
        //ConstStepGradient g;
        //g.R1Minimizer().setCallback(new Problem2HNDirichletR1MinimizeCallback);
        g.setFunction(&prob);
        g.setGradient(&prob);
        g.setPrinter(&prob);
        g.setProjection(&prob);
        g.setEpsilon1(0.0001);
        g.setEpsilon2(0.0001);
        g.setEpsilon3(0.0001);
        g.setR1MinimizeEpsilon(0.1, 0.001);
        g.setNormalize(false);
        g.showExitMessage(true);

        g.calculate(x);

        IPrinter::printSeperatorLine(NULL, '=');
    }
}

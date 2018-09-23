#include "problem2h_example.h"

void example1()
{
    // Equation parameters
    EquationParameter e_prm;
    e_prm.a = 1.0;
    e_prm.lambda = +0.01;

    // Pulse influences
    e_prm.Ns = 2;
    e_prm.q.resize(e_prm.Ns);
    e_prm.theta.resize(e_prm.Ns);

    e_prm.q[0] = +0.012; e_prm.theta[0].x = 0.4300; e_prm.theta[0].y = 0.7500;
    e_prm.q[1] = +0.013; e_prm.theta[1].x = 0.8700; e_prm.theta[1].y = 0.2300;

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
    double ht = 0.005; int Nt = 200;

    Dimension time(ht, 0, Nt);
    Dimension dimx(hx, 0, Nx);
    Dimension dimy(hy, 0, Ny);

    // Penalty paramteres
    DoubleVector r; r << 0.0000 << 0.0100 << 0.10000 << 1.0000;
    // Regularization coefficients
    DoubleVector e; e << 0.0000 << 0.0000 << 0.00000 << 0.0000;

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
        prob.vmin.resize(e_prm.Nc, +0.004);
        prob.vmax.resize(e_prm.Nc, +0.008);
        prob.LD = 50;

        prob.regEpsilon = e[i];
        prob.r = r[i];
        if (i==0)
        {
            prob.PrmToVector(o_prm, x);
            prob.checkGradient(prob);
            IPrinter::printSeperatorLine();
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
        g.setR1MinimizeEpsilon(1.0, 0.01);
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
        //            double a = 0.001*i;
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

void example2()
{
    // Equation parameters
    EquationParameter e_prm;
    e_prm.a = +1.0;
    e_prm.lambda = +0.01;

    // Pulse influences
    e_prm.Ns = 5;
    e_prm.q.resize(e_prm.Ns);
    e_prm.theta.resize(e_prm.Ns);

    //e_prm.q[0] = +0.12; e_prm.theta[0].x = 0.4300; e_prm.theta[0].y = 0.7500;
    //e_prm.q[1] = +0.13; e_prm.theta[1].x = 0.8700; e_prm.theta[1].y = 0.2300;
    e_prm.q[0] = +0.12; e_prm.theta[0].x = 0.5000; e_prm.theta[0].y = 0.5000;
    e_prm.q[1] = +0.12; e_prm.theta[1].x = 0.2500; e_prm.theta[1].y = 0.2500;
    e_prm.q[2] = +0.12; e_prm.theta[2].x = 0.7500; e_prm.theta[2].y = 0.7500;
    e_prm.q[3] = +0.12; e_prm.theta[3].x = 0.2500; e_prm.theta[3].y = 0.7500;
    e_prm.q[4] = +0.12; e_prm.theta[4].x = 0.7500; e_prm.theta[4].y = 0.2500;

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

    o_prm.k[0][0]  = +1.1200; o_prm.k[0][1]  = +1.2400; o_prm.k[1][0]  = +1.4500; o_prm.k[1][1]  = +1.1800;
    o_prm.z[0][0]  = +0.5000; o_prm.z[0][1]  = +0.4000; o_prm.z[1][0]  = +0.7000; o_prm.z[1][1]  = +0.5000;
    o_prm.xi[0].x  = +0.4274; o_prm.xi[0].y  = +0.6735; o_prm.xi[1].x  = +0.6710; o_prm.xi[1].y  = +0.3851;
    o_prm.eta[0].x = +0.5174; o_prm.eta[0].y = +0.7635; o_prm.eta[1].x = +0.5570; o_prm.eta[1].y = +0.4751;

    r_prm = o_prm;

    // Grid parameters
    double hx = 0.010; int Nx = 100;
    double hy = 0.010; int Ny = 100;
    double ht = 0.010; int Nt = 100;

    Dimension time(ht, 0, Nt);
    Dimension dimx(hx, 0, Nx);
    Dimension dimy(hy, 0, Ny);

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
    prob.vmin.resize(e_prm.Nc, +0.004);
    prob.vmax.resize(e_prm.Nc, +0.008);
    prob.LD = 1;
    prob.regEpsilon = 0.0;

    std::vector<DoubleMatrix> u;
    spif_vector info;
    prob.solveForwardIBVP(u, info, false);

//    IPrinter::printSeperatorLine();
//    DoubleMatrix &m = u.at(0);
//    IPrinter::printMatrix(m);
//    IPrinter::printSeperatorLine();

//    printf("%f\n", prob.integralU(m));
}


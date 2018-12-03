#include "problem2h_example.h"
#include "vectornormalizer.h"
#ifdef USE_IMAGING
#include <imaging.h>
#include <QtGui/QGuiApplication>
#endif

void Problem2HNDirichlet::Main(int argc UNUSED_PARAM, char *argv[] UNUSED_PARAM)
{
#ifdef USE_IMAGING
    QGuiApplication app(argc, argv);
#endif
    //example1();
    //example2();
    example3();
}

void example1()
{
    // Equation parameters
    EquationParameterH e_prm;
    e_prm.a = 1.0;
    e_prm.lambda = +0.00;

    // Pulse influences
    e_prm.Ns = 2;
    e_prm.q.resize(e_prm.Ns);
    e_prm.theta.resize(e_prm.Ns);

    e_prm.q[0] = +0.145; e_prm.theta[0].x = 0.2500; e_prm.theta[0].y = 0.7200;
    e_prm.q[1] = +0.157; e_prm.theta[1].x = 0.5400; e_prm.theta[1].y = 0.2700;
    //e_prm.q[2] = +0.148; e_prm.theta[2].x = 0.7400; e_prm.theta[2].y = 0.6300;

    e_prm.No = 2;
    e_prm.Nc = 2;

    // Optimization parameters
    OptimizeParameterH o_prm;
    o_prm.k.resize(e_prm.Nc, e_prm.No, 0.0);
    o_prm.z.resize(e_prm.Nc, e_prm.No, 0.0);
    o_prm.xi.resize(e_prm.No);
    o_prm.eta.resize(e_prm.Nc);

    o_prm.k[0][0]  = -0.1262; o_prm.k[0][1]  = -0.4038; o_prm.k[1][0]  = -1.7431; o_prm.k[1][1]  = -2.8052;
    o_prm.z[0][0]  = -2.0245; o_prm.z[0][1]  = +3.0784; o_prm.z[1][0]  = -5.0587; o_prm.z[1][1]  = +8.0641;
    o_prm.xi[0].x  = +0.6527; o_prm.xi[0].y  = +0.8412; o_prm.xi[1].x  = +0.7412; o_prm.xi[1].y  = +0.2483;
    o_prm.eta[0].x = +0.3254; o_prm.eta[0].y = +0.3654; o_prm.eta[1].x = +0.9462; o_prm.eta[1].y = +0.4966;


    // Regularization parameters
    OptimizeParameterH r_prm;
    r_prm.k.resize(e_prm.Nc, e_prm.No, 0.0);
    r_prm.z.resize(e_prm.Nc, e_prm.No, 0.0);
    r_prm.xi.resize(e_prm.No);
    r_prm.eta.resize(e_prm.Nc);

    r_prm.k[0][0]  = -0.0182; r_prm.k[0][1]  = -0.0125; r_prm.k[1][0]  = -0.0155; r_prm.k[1][1]  = -0.0131;
    r_prm.z[0][0]  = -0.0262; r_prm.z[0][1]  = -0.0773; r_prm.z[1][0]  = -0.0570; r_prm.z[1][1]  = +0.0653;
    r_prm.xi[0].x  = +0.3849; r_prm.xi[0].y  = +0.5442; r_prm.xi[1].x  = +0.7661; r_prm.xi[1].y  = +0.6785;
    r_prm.eta[0].x = +0.6656; r_prm.eta[0].y = +0.7909; r_prm.eta[1].x = +0.4856; r_prm.eta[1].y = +0.3810;

    //r_prm = o_prm;

    // Grid parameters
    double hx = 0.010; int Nx = 100;
    double hy = 0.010; int Ny = 100;
    double ht = 0.010; int Nt = 200;

    Dimension time(ht, 0, Nt);
    Dimension dimx(hx, 0, Nx);
    Dimension dimy(hy, 0, Ny);

    // Penalty paramteres
    DoubleVector r; r << 0.10 << 1.0 << 10.0 << 100.00;
    // Regularization coefficients
    DoubleVector e; e << 1.00 << 0.10 << 0.010 << 0.00100;

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
        prob.optimizeO = true;
        prob.optimizeC = true;
        prob.vmin.resize(e_prm.Nc, -0.005);
        prob.vmax.resize(e_prm.Nc, +0.005);
        prob.LD = 20;

        prob.regEpsilon = e[i];
        prob.r = r[i];
        if (i==0)
        {
            prob.PrmToVector(o_prm, x);
            //prob.checkGradient1(prob);
            IPrinter::printSeperatorLine();

            //            std::vector<DoubleMatrix> u;
            //            spif_vectorH u_info;
            //            prob.solveForwardIBVP(u, u_info, false);
            //            return;
        }

        //ConjugateGradient g;
        SteepestDescentGradient g;
        g.setFunction(&prob);
        g.setGradient(&prob);
        g.setPrinter(&prob);
        g.setProjection(&prob);
        //g.setGradientNormalizer(&prob);
        g.setOptimalityTolerance(0.00001);
        g.setFunctionTolerance(0.00001);
        g.setStepTolerance(0.00001);
        g.setR1MinimizeEpsilon(0.1, 0.01);
        g.setMaxIterations(50);
        g.setNormalize(true);
        g.showExitMessage(true);
        prob.gm = &g;

        g.calculate(x);

        IPrinter::printSeperatorLine(nullptr, '=');
    }
}

void example2()
{
    // Equation parameters
    EquationParameterH e_prm;
    e_prm.a = 1.0;
    e_prm.lambda = +0.00;

    // Pulse influences
    e_prm.Ns = 2;
    e_prm.q.resize(e_prm.Ns);
    e_prm.theta.resize(e_prm.Ns);

    e_prm.q[0] = +0.145; e_prm.theta[0].x = 0.2500; e_prm.theta[0].y = 0.7200;
    e_prm.q[1] = +0.157; e_prm.theta[1].x = 0.5400; e_prm.theta[1].y = 0.2700;

    e_prm.No = 2;
    e_prm.Nc = 2;

    // Optimization parameters
    OptimizeParameterH o_prm;
    o_prm.k.resize(e_prm.Nc, e_prm.No, 0.0);
    o_prm.z.resize(e_prm.Nc, e_prm.No, 0.0);
    o_prm.xi.resize(e_prm.No);
    o_prm.eta.resize(e_prm.Nc);

    o_prm.k[0][0]  = -1.0345; o_prm.k[0][1]  = -2.5401; o_prm.k[1][0]  = -1.1431; o_prm.k[1][1]  = -1.0984;
    o_prm.z[0][0]  = +3.1245; o_prm.z[0][1]  = +1.8532; o_prm.z[1][0]  = -2.4512; o_prm.z[1][1]  = +2.5421;
    o_prm.xi[0].x  = +0.7524; o_prm.xi[0].y  = +0.4828; o_prm.xi[1].x  = +0.1274; o_prm.xi[1].y  = +0.8234;
    o_prm.eta[0].x = +0.8635; o_prm.eta[0].y = +0.3654; o_prm.eta[1].x = +0.2496; o_prm.eta[1].y = +0.6536;


    // Regularization parameters
    OptimizeParameterH r_prm;
    r_prm.k.resize(e_prm.Nc, e_prm.No, 0.0);
    r_prm.z.resize(e_prm.Nc, e_prm.No, 0.0);
    r_prm.xi.resize(e_prm.No);
    r_prm.eta.resize(e_prm.Nc);

    r_prm.k[0][0]  = -0.0345; r_prm.k[0][1]  = -0.0541; r_prm.k[1][0]  = -0.0431; r_prm.k[1][1]  = -0.0984;
    r_prm.z[0][0]  = +0.1245; r_prm.z[0][1]  = +0.0325; r_prm.z[1][0]  = -0.1452; r_prm.z[1][1]  = +0.2154;
    r_prm.xi[0].x  = +0.2754; r_prm.xi[0].y  = +0.2438; r_prm.xi[1].x  = +0.4271; r_prm.xi[1].y  = +0.3824;
    r_prm.eta[0].x = +0.3865; r_prm.eta[0].y = +0.4365; r_prm.eta[1].x = +0.6294; r_prm.eta[1].y = +0.5466;

    //r_prm = o_prm;

    // Grid parameters
    double hx = 0.010; int Nx = 100;
    double hy = 0.010; int Ny = 100;
    double ht = 0.010; int Nt = 200;

    Dimension time(ht, 0, Nt);
    Dimension dimx(hx, 0, Nx);
    Dimension dimy(hy, 0, Ny);

    // Penalty paramteres
    DoubleVector r; r << 0.10 << 1.0 << 10.0 << 100.00;
    // Regularization coefficients
    DoubleVector e; e << 1.00 << 0.10 << 0.010 << 0.00100;

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
        prob.optimizeO = true;
        prob.optimizeC = true;
        prob.vmin.resize(e_prm.Nc, -0.005);
        prob.vmax.resize(e_prm.Nc, +0.005);
        prob.LD = 20;

        prob.regEpsilon = e[i];
        prob.r = r[i];
        if (i==0)
        {
            prob.PrmToVector(o_prm, x);
            //prob.checkGradient1(prob);
            IPrinter::printSeperatorLine();

            //            std::vector<DoubleMatrix> u;
            //            spif_vectorH u_info;
            //            prob.solveForwardIBVP(u, u_info, false);
            //            return;
        }

        //ConjugateGradient g;
        SteepestDescentGradient g;
        g.setFunction(&prob);
        g.setGradient(&prob);
        g.setPrinter(&prob);
        g.setProjection(&prob);
        //g.setGradientNormalizer(&prob);
        g.setOptimalityTolerance(0.00001);
        g.setFunctionTolerance(0.00001);
        g.setStepTolerance(0.00001);
        g.setR1MinimizeEpsilon(0.1, 0.01);
        g.setMaxIterations(50);
        g.setNormalize(true);
        g.showExitMessage(true);
        prob.gm = &g;

        g.calculate(x);

        IPrinter::printSeperatorLine(nullptr, '=');
    }
}

void example3()
{
    // Equation parameters
    EquationParameterH e_prm;
    e_prm.a = 1.0;
    e_prm.lambda = +0.00;

    // Pulse influences
    e_prm.Ns = 2;
    e_prm.q.resize(e_prm.Ns);
    e_prm.theta.resize(e_prm.Ns);

    e_prm.q[0] = +0.145; e_prm.theta[0].x = 0.2500; e_prm.theta[0].y = 0.7200;
    e_prm.q[1] = +0.157; e_prm.theta[1].x = 0.5400; e_prm.theta[1].y = 0.2700;

    e_prm.No = 2;
    e_prm.Nc = 2;

    // Optimization parameters
    OptimizeParameterH o_prm;
    o_prm.k.resize(e_prm.Nc, e_prm.No, 0.0);
    o_prm.z.resize(e_prm.Nc, e_prm.No, 0.0);
    o_prm.xi.resize(e_prm.No);
    o_prm.eta.resize(e_prm.Nc);

    o_prm.k[0][0]  = -0.9450; o_prm.k[0][1]  = -1.4520; o_prm.k[1][0]  = -0.1431; o_prm.k[1][1]  = -0.9184;
    o_prm.z[0][0]  = +1.4125; o_prm.z[0][1]  = -3.2285; o_prm.z[1][0]  = -0.4512; o_prm.z[1][1]  = +0.8468;
    o_prm.xi[0].x  = +0.5274; o_prm.xi[0].y  = +0.8248; o_prm.xi[1].x  = +0.2814; o_prm.xi[1].y  = +0.2384;
    o_prm.eta[0].x = +0.6458; o_prm.eta[0].y = +0.4635; o_prm.eta[1].x = +0.8926; o_prm.eta[1].y = +0.5368;


    // Regularization parameters
    OptimizeParameterH r_prm;
    r_prm.k.resize(e_prm.Nc, e_prm.No, 0.0);
    r_prm.z.resize(e_prm.Nc, e_prm.No, 0.0);
    r_prm.xi.resize(e_prm.No);
    r_prm.eta.resize(e_prm.Nc);

    r_prm.k[0][0]  = -0.0163; r_prm.k[0][1]  = -0.0274; r_prm.k[1][0]  = -0.0174; r_prm.k[1][1]  = -0.0474;
    r_prm.z[0][0]  = +0.0874; r_prm.z[0][1]  = +0.0985; r_prm.z[1][0]  = +0.1984; r_prm.z[1][1]  = +0.3856;
    r_prm.xi[0].x  = +0.7524; r_prm.xi[0].y  = +0.4982; r_prm.xi[1].x  = +0.1247; r_prm.xi[1].y  = +0.8243;
    r_prm.eta[0].x = +0.5638; r_prm.eta[0].y = +0.3865; r_prm.eta[1].x = +0.2946; r_prm.eta[1].y = +0.6465;

    //o_prm = r_prm;

    // Grid parameters
    double hx = 0.010; int Nx = 100;
    double hy = 0.010; int Ny = 100;
    double ht = 0.010; int Nt = 200;

    Dimension time(ht, 0, Nt);
    Dimension dimx(hx, 0, Nx);
    Dimension dimy(hy, 0, Ny);

    // Penalty paramteres
    DoubleVector r; r << 0.10 << 1.0 << 10.0 << 100.00;
    // Regularization coefficients
    DoubleVector e; e << 1.00 << 0.10 << 0.010 << 0.00100;

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
        prob.optimizeO = true;
        prob.optimizeC = true;
        prob.vmin.resize(e_prm.Nc, -0.005);
        prob.vmax.resize(e_prm.Nc, +0.005);
        prob.LD = 20;

        prob.regEpsilon = e[i];
        prob.r = r[i];
        if (i==0)
        {
            prob.PrmToVector(o_prm, x);
            //prob.checkGradient1(prob);
            IPrinter::printSeperatorLine();

            //            std::vector<DoubleMatrix> u;
            //            spif_vectorH u_info;
            //            prob.solveForwardIBVP(u, u_info, false);
            //            return;
        }

        //ConjugateGradient g;
        SteepestDescentGradient g;
        g.setFunction(&prob);
        g.setGradient(&prob);
        g.setPrinter(&prob);
        g.setProjection(&prob);
        //g.setGradientNormalizer(&prob);
        g.setOptimalityTolerance(0.00001);
        g.setFunctionTolerance(0.00001);
        g.setStepTolerance(0.00001);
        g.setR1MinimizeEpsilon(0.1, 0.01);
        g.setMaxIterations(20);
        g.setNormalize(true);
        g.showExitMessage(true);
        prob.gm = &g;

        g.calculate(x);

        IPrinter::printSeperatorLine(nullptr, '=');
    }
}

#include "problem2h_example.h"

void Problem2HDirichlet::Main(int argc UNUSED_PARAM, char *argv[] UNUSED_PARAM)
{
#ifdef USE_IMAGING
    QGuiApplication app(argc, argv);
#endif
    prod_example1();
    //example1();
    //example2();
    //example3();
}

void prod_example1()
{
    // Equation parameters
    EquationParameterH e_prm;
    e_prm.a = 1.0;
    e_prm.alpha = +0.00;

    e_prm.Q1 << 0.21;// << 0.22 << 0.24;
    e_prm.Q2 << 0.25;// << 0.27 << 0.29;

#if defined(DISCRETE_DELTA_TIME)
    e_prm.Nt = 10;
    e_prm.timeMoments.resize(e_prm.Nt);
    for (unsigned int s=0; s<e_prm.Nt; s++) e_prm.timeMoments[s] = (s)*0.3;
#endif

    // Pulse influences
    e_prm.Ns = 2;
    e_prm.pulses.resize(e_prm.Ns);

    // Optimization parameters
    OptimizeParameterH o_prm;
#if defined(DISCRETE_DELTA_TIME)
    e_prm.Nc = 2;
    e_prm.No = 2;
    o_prm.k = new DoubleMatrix[e_prm.Nt]; for (unsigned int s=0; s<e_prm.Nt; s++) o_prm.k[s].resize(e_prm.Nc, e_prm.No, 0.00);
    o_prm.z = new DoubleMatrix[e_prm.Nt]; for (unsigned int s=0; s<e_prm.Nt; s++) o_prm.z[s].resize(e_prm.Nc, e_prm.No, 0.00);
    o_prm.xi.resize(e_prm.No);
    o_prm.eta.resize(e_prm.Nc);
    for (unsigned int s=0; s<e_prm.Nt; s++)
    {
        for (unsigned int r=0; r<e_prm.Nc; r++)
        {
            for (unsigned int c=0; c<e_prm.No; c++)
            {
                //o_prm.k[s][r][c] = -static_cast<double>((rand() % 1000))/1000.0;
                //o_prm.k[s][r][c] = 1.0-static_cast<double>((rand() % 2000))/1000.0;
                //o_prm.z[s][r][c] = +static_cast<double>((rand() % 1000))/100000.0;
                //o_prm.k[s][r][c] = -fabs(sin((c+1)*10.0)*cos((r+1)*20.0)*sin((s+1)*0.1));
                //o_prm.z[s][r][c] = cos((c+1)*10.0)*sin((r+1)*20.0)*sin((s+1)*0.2);

                o_prm.k[s][r][c] = -0.05;
                o_prm.z[s][r][c] = +0.05;
            }
        }
    }
#else
    e_prm.Nc = 2;
    e_prm.No = 2;
    o_prm.k.resize(e_prm.Nc, e_prm.No, 0.0);
    o_prm.z.resize(e_prm.Nc, e_prm.No, 0.0);
    o_prm.xi.resize(e_prm.No);
#endif

    // 0.003691
    e_prm.pulses[0] = InitialPulse2D(SpacePoint(0.2800, 0.3200), 0.21);
    e_prm.pulses[1] = InitialPulse2D(SpacePoint(0.7300, 0.6500), 0.25);
    //o_prm.k[0][0]  = -2.0610; o_prm.k[0][1]  = -2.9376; o_prm.k[1][0]  = -2.1707; o_prm.k[1][1]  = -2.8527;
    //o_prm.k[0][0]  = -1.0000; o_prm.k[0][1]  = -1.0000; o_prm.k[1][0]  = -1.0000; o_prm.k[1][1]  = -1.0000;
    //o_prm.z[0][0]  = -0.0461; o_prm.z[0][1]  = -0.0246; o_prm.z[1][0]  = +0.0319; o_prm.z[1][1]  = +0.0161;
    o_prm.xi[0].x  = +0.5400; o_prm.xi[0].y  = +0.3500; o_prm.xi[1].x  = +0.8200; o_prm.xi[1].y  = +0.9400;
    o_prm.eta[0].x = +0.2000; o_prm.eta[0].y = +0.6000; o_prm.eta[1].x = +0.7000; o_prm.eta[1].y = +0.4000;

    // Regularization parameters
    OptimizeParameterH r_prm;
#if defined(DISCRETE_DELTA_TIME)
    r_prm.k = new DoubleMatrix[e_prm.Nt]; for (unsigned int s=0; s<e_prm.Nt; s++) r_prm.k[s].resize(e_prm.Nc, e_prm.No, 0.0);
    r_prm.z = new DoubleMatrix[e_prm.Nt]; for (unsigned int s=0; s<e_prm.Nt; s++) r_prm.z[s].resize(e_prm.Nc, e_prm.No, 0.0);
    r_prm.xi.resize(e_prm.No);
    r_prm.eta.resize(e_prm.Nc);
    for (unsigned int s=0; s<e_prm.Nt; s++)
    {
        for (unsigned int r=0; r<e_prm.Nc; r++)
        {
            for (unsigned int c=0; c<e_prm.No; c++)
            {
                //r_prm.k[s][r][c] = 1.0-static_cast<double>((rand() % 2000))/1000.0;
                //r_prm.z[s][r][c] = +static_cast<double>((rand() % 1000))/100000.0;
            }
        }
    }
#else
    r_prm.k.resize(e_prm.Nc, e_prm.No, 0.0);
    r_prm.z.resize(e_prm.Nc, e_prm.No, 0.0);
    r_prm.xi.resize(e_prm.No);
    r_prm.eta.resize(e_prm.Nc);
#endif

    // Penalty paramteres
    DoubleVector r; r << 0.0000;// << 0.0000 << 0.000;
    // Regularization coefficients
    DoubleVector e; e << 0.0000 << 0.0000 << 0.0000;
    DoubleVector e1; e1 << 0.1000 << 0.0100 << 0.0010;
    DoubleVector e2; e2 << 0.0100 << 0.0010 << 0.0001;

    DoubleVector x;
    for (unsigned int i=0; i<r.length(); i++)
    {
        Problem2HDirichlet1 prob;
        prob.setGridDimensions(Dimension(0.010, 0, 300),
                               Dimension(0.010, 0, 100),
                               Dimension(0.010, 0, 100));
        prob.mEquParameter = e_prm;
        prob.mOptParameter = o_prm;
        prob.mRegParameter = r_prm;
        prob.optimizeK = true;
        prob.optimizeZ = true;
        prob.optimizeO = true;
        prob.optimizeC = true;
        prob.vmin.resize(e_prm.Nc, -0.5);
        prob.vmax.resize(e_prm.Nc, +0.5);
        prob.LD = 30;
        prob.noise = 0.0;

        prob.regEpsilon = e[i];
        prob.r = r[i];

        if (i==0)
        {
            prob.PrmToVector(o_prm, x);
            //prob.checkGradient1(prob);
            //prob.checkGradient2(prob);
            //return;

            //prob.printLayers = false;
            //if (prob.printLayers)
            //{
            //    std::vector<DoubleVector> u;
            //    spif_vector1H u_info, p_info;
            //    prob.solveForwardIBVP(u, u_info, true);
            //    return;
            //}
            //IPrinter::printSeperatorLine();
            //prob.solveBackwardIBVP(u, p_info, true, u_info);
            //return;
        }

        //ConjugateGradient g;
        SteepestDescentGradient g;
        g.setFunction(&prob);
        g.setGradient(&prob);
        g.setPrinter(&prob);
        g.setProjection(&prob);
        //g.setProjection(new ProjectionEx1);
        //g.setGradientNormalizer(&prob);
        g.setOptimalityTolerance(0.00001);
        g.setFunctionTolerance(0.00001);
        g.setStepTolerance(0.00001);
        g.setR1MinimizeEpsilon(e1[i], e2[i]);
        g.setMaxIterations(50);
        g.setNormalize(false);
        g.showExitMessage(true);
        //prob.gm = &g;
        g.calculate(x);

        IPrinter::printSeperatorLine(nullptr, '=');

        IPrinter::print(x, x.length(), 6, 4);
        prob.VectorToPrm(x, o_prm);
        prob.mOptParameter = o_prm;
        prob.printLayers = true;
        if (prob.printLayers)
        {
            prob.setGridDimensions(Dimension(0.010, 0, 800), Dimension(0.010, 0, 100), Dimension(0.010, 0, 100));
            std::vector<DoubleMatrix> u;
            spif_vectorH u_info, p_info;
            prob.solveForwardIBVP(u, u_info, true, x, 0.25);
        }
        puts("End");
    }
}

void example1() {}

void example2() {}

void example3() {}

void prod_figure1(Problem2HDirichlet &prob, double ht, const DoubleVector &x)
{
    for (int i=0; i<=1000; i++)
    {
        //if (i%10==0)
        {
            prob.setTimeDimension(Dimension(ht, 0, i));
            printf("ln:%d fx:%8.6f\n", i, prob.fx(x));
        }
    }
}

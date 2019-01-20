#include "problem1h_example.h"

struct ProjectionEx1 : public IProjection
{
    virtual void project(DoubleVector &x, unsigned int index);
    virtual void project(DoubleVector &) const;
};

void ProjectionEx1::project(DoubleVector &x) const
{
    const double mn = 0.05;
    const double mx = 0.95;
    for (unsigned int i=8; i<=15; i++)
    {
        if (x[i]<mn) x[i] = mn;
        if (x[i]>mx) x[i] = mx;
    }

    return;


    //    if ( x[ 8] < 0.05 ) x[ 8] = 0.05; if ( x[ 8] > 0.45 ) x[ 8] = 0.45;
    //    if ( x[ 9] < 0.55 ) x[ 9] = 0.55; if ( x[ 9] > 0.95 ) x[ 9] = 0.95;
    //    if ( x[10] < 0.55 ) x[10] = 0.55; if ( x[10] > 0.95 ) x[10] = 0.95;
    //    if ( x[11] < 0.05 ) x[11] = 0.05; if ( x[11] > 0.45 ) x[11] = 0.45;

    //    if ( x[12] < 0.05 ) x[12] = 0.05; if ( x[12] > 0.45 ) x[12] = 0.45;
    //    if ( x[13] < 0.05 ) x[13] = 0.05; if ( x[13] > 0.45 ) x[13] = 0.45;
    //    if ( x[14] < 0.55 ) x[14] = 0.55; if ( x[14] > 0.95 ) x[14] = 0.95;
    //    if ( x[15] < 0.55 ) x[15] = 0.55; if ( x[15] > 0.95 ) x[15] = 0.95;

    //    DoubleVector x0 = x.mid(8,x.length()-1);
    //    IPrinter::print(x0);
}

void ProjectionEx1::project(DoubleVector &, unsigned int) {}

void Problem1HDirichletBase::Main(int argc UNUSED_PARAM, char *argv[] UNUSED_PARAM)
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
    EquationParameter1H e_prm;
    e_prm.a = 1.0;
    e_prm.alpha = +0.00;

    // Pulse influences
    e_prm.Ns = 2;
    e_prm.theta.resize(e_prm.Ns);

    // Optimization parameters
    OptimizeParameter1H o_prm;
    e_prm.Nc = 2;
    e_prm.No = 2;
    o_prm.k.resize(e_prm.Nc, e_prm.No, 0.0);
    o_prm.z.resize(e_prm.Nc, e_prm.No, 0.0);
    o_prm.ksi.resize(e_prm.No);
    o_prm.eta.resize(e_prm.Nc);

    e_prm.theta[0].q = +0.114; e_prm.theta[0].x = 0.2845;
    e_prm.theta[1].q = +0.128; e_prm.theta[1].x = 0.7382;

    // (70.749111 200)
    //o_prm.k[0][0]  = -2.0610; o_prm.k[0][1]  = -2.9376; o_prm.k[1][0]  = -2.1707; o_prm.k[1][1]  = -2.8527;
    //o_prm.z[0][0]  = -0.0461; o_prm.z[0][1]  = -0.0246; o_prm.z[1][0]  = +0.0319; o_prm.z[1][1]  = +0.0161;
    //o_prm.k[0][0]  = -2.5900; o_prm.k[0][1]  = -1.8100; o_prm.k[1][0]  = -1.5700; o_prm.k[1][1]  = -2.76000;
    o_prm.ksi[0].x = +0.5400; o_prm.ksi[1].x = +0.8200; o_prm.eta[0].x = +0.1200; o_prm.eta[1].x = +0.3400;
    // 0.011497
    //o_prm.k[0][0]  = -0.7099; o_prm.k[0][1]  = -1.7252; o_prm.k[1][0]  = -0.0766; o_prm.k[1][1]  = -0.7024;
    //o_prm.z[0][0]  = +0.9573; o_prm.z[0][1]  = +1.4058; o_prm.z[1][0]  = +1.7626; o_prm.z[1][1]  = +2.2906;
    //o_prm.ksi[0].x = +0.0500; o_prm.ksi[1].x = +0.9500; o_prm.eta[0].x = +0.0713; o_prm.eta[1].x = +0.9500;
    // 0.013009
    //o_prm.k[0][0]  = -1.0584; o_prm.k[0][1]  = -0.5405; o_prm.k[1][0]  = -1.4563; o_prm.k[1][1]  = -0.9511;
    //o_prm.z[0][0]  = +0.0148; o_prm.z[0][1]  = +0.0033; o_prm.z[1][0]  = +0.0116; o_prm.z[1][1]  = -0.0467;
    //o_prm.ksi[0].x = +0.8140; o_prm.ksi[1].x = +0.1057; o_prm.eta[0].x = +0.5469; o_prm.eta[1].x = +0.5550;
    // 0.013374
    //o_prm.k[0][0]  = -2.4222; o_prm.k[0][1]  = -3.2657; o_prm.k[1][0]  = -2.7579; o_prm.k[1][1]  = -3.4396;
    //o_prm.z[0][0]  = -0.0116; o_prm.z[0][1]  = +0.0072; o_prm.z[1][0]  = +0.0181; o_prm.z[1][1]  = -0.0154;
    //o_prm.ksi[0].x = +0.8291; o_prm.ksi[1].x = +0.1942; o_prm.eta[0].x = +0.8983; o_prm.eta[1].x = +0.1548;
    // 0.002373 50
    //o_prm.k[0][0]  = -1.6920; o_prm.k[0][1]  = -2.6234; o_prm.k[1][0]  = -1.5632; o_prm.k[1][1]  = -2.2822;
    //o_prm.z[0][0]  = +0.2738; o_prm.z[0][1]  = +0.4353; o_prm.z[1][0]  = +0.3810; o_prm.z[1][1]  = +0.4680;
    //o_prm.ksi[0].x = +0.0500; o_prm.ksi[1].x = +0.9500; o_prm.eta[0].x = +0.0537; o_prm.eta[1].x = +0.9500;

    //-----------------------------------------------------------------------------------------------------

    // Regularization parameters
    OptimizeParameter1H r_prm;
    r_prm.k.resize(e_prm.Nc, e_prm.No, 0.0);
    r_prm.z.resize(e_prm.Nc, e_prm.No, 0.0);
    r_prm.ksi.resize(e_prm.No);
    r_prm.eta.resize(e_prm.Nc);

    // Penalty paramteres
    DoubleVector r; r << 0.0000;// << 0.0000 << 0.000;
    // Regularization coefficients
    DoubleVector e; e << 0.0000 << 0.0000 << 0.0000;
    DoubleVector e1; e1 << 0.1000 << 0.0100 << 0.0010;
    DoubleVector e2; e2 << 0.0100 << 0.0010 << 0.0001;

    DoubleVector x;
    for (unsigned int i=0; i<r.length(); i++)
    {
        Problem1HDirichlet1 prob;
        prob.setGridDimensions(Dimension(0.010, 0, 500), Dimension(0.010, 0, 100));
        prob.mEquParameter = e_prm;
        prob.mOptParameter = o_prm;
        prob.mRegParameter = r_prm;
        prob.optimizeK = true;
        prob.optimizeZ = true;
        prob.optimizeO = true;
        prob.optimizeC = true;
        prob.vmin.resize(e_prm.Nc, -0.5);
        prob.vmax.resize(e_prm.Nc, +0.5);
        prob.LD = 50;
        prob.noise = 0.0;

        prob.regEpsilon = e[i];
        prob.r = r[i];

        if (i==0)
        {
            prob.PrmToVector(o_prm, x);
            //prob.checkGradient1(prob);
            //prob.checkGradient2(prob);
            //IPrinter::printSeperatorLine();
            //return;

            prob.printLayers = false;
            if (prob.printLayers)
            {
                std::vector<DoubleVector> u;
                spif_vector1H u_info, p_info;
                prob.solveForwardIBVP(u, u_info, true);
                return;
            }
            //IPrinter::printSeperatorLine();
            //prob.solveBackwardIBVP(u, p_info, true, u_info);
            //return;
        }

        IPrinter::print(x, x.length(), 6, 4);
        ConjugateGradient g;
        //SteepestDescentGradient g;
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
        //g.setMaxIterations(50);
        g.setNormalize(false);
        g.showExitMessage(true);
        //prob.gm = &g;
        g.calculate(x);
        //IPrinter::print(x, x.length(), 6, 4);

        IPrinter::printSeperatorLine(nullptr, '=');
        //return;

        //x.clear();
        //0.012558
        //x << -0.9747 << -0.2779 << 1.1427 << -0.1328 << 0.4921 << 0.2416 << -0.1491 << 0.2960 << 0.1881 << 0.9500 << 0.2076 << 0.6641;
        //0.003731
        //x << -1.5687 << -0.6123 << 1.5319 << -0.4288 << 1.2693 << 0.4855 << -1.0917 << 0.4616 <<  0.5381 <<  0.3260 <<  0.4519 << 0.7312;
        //1.720225
        //x << -2.5900 << -1.8100 << -1.5700 << -2.7600 << 0.0000  << 0.0000 <<  0.0000 <<  0.0000 <<  0.5400 << 0.8200 <<  0.1200  << 0.3400;

        IPrinter::print(x, x.length(), 6, 4);
        prob.VectorToPrm(x, o_prm);
        prob.mOptParameter = o_prm;
        prob.printLayers = true;
        if (prob.printLayers)
        {
            prob.setGridDimensions(Dimension(0.010, 0, 1000), Dimension(0.010, 0, 100));
            std::vector<DoubleVector> u;
            spif_vector1H u_info, p_info;
            prob.solveForwardIBVP(u, u_info, true);
            return;
        }
    }
}

void prod_figure1(Problem1HDirichletBase &prob, double ht, const DoubleVector &x)
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

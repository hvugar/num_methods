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

    e_prm.theta[0].q = +0.214; e_prm.theta[0].x = 0.2845; e_prm.theta[0].y = 0.2845;
    e_prm.theta[0].q = +0.228; e_prm.theta[1].x = 0.7382; e_prm.theta[1].y = 0.6518;
    //o_prm.k[0][0]  = -2.0610; o_prm.k[0][1]  = -2.9376; o_prm.k[1][0]  = -2.1707; o_prm.k[1][1]  = -2.8527;
    //o_prm.z[0][0]  = -0.0461; o_prm.z[0][1]  = -0.0246; o_prm.z[1][0]  = +0.0319; o_prm.z[1][1]  = +0.0161;
    o_prm.ksi[0].x = +0.5400; o_prm.ksi[1].x = +0.8200; o_prm.eta[0].x = +0.1200; o_prm.eta[1].x = +0.3400;

    // 0.010242
    //e_prm.q[0] = +0.314; e_prm.theta[0].x = 0.2845; e_prm.theta[0].y = 0.2845;
    //e_prm.q[1] = +0.328; e_prm.theta[1].x = 0.7382; e_prm.theta[1].y = 0.6518;
    //o_prm.k[0][0]  = -2.0610; o_prm.k[0][1]  = -1.9376; o_prm.k[1][0]  = -2.1707; o_prm.k[1][1]  = -1.8527;
    //o_prm.z[0][0]  = -0.0461; o_prm.z[0][1]  = -0.0246; o_prm.z[1][0]  = +0.0319; o_prm.z[1][1]  = +0.0161;
    //o_prm.xi[0].x  = +0.6402; o_prm.xi[0].y  = +0.4899; o_prm.xi[1].x  = +0.8148; o_prm.xi[1].y  = +0.9500;
    //o_prm.eta[0].x = +0.9500; o_prm.eta[0].y = +0.4609; o_prm.eta[1].x = +0.2512; o_prm.eta[1].y = +0.8932;

    // 0.006407
    //e_prm.q[0] = +0.314; e_prm.theta[0].x = 0.2845; e_prm.theta[0].y = 0.2845;
    //e_prm.q[1] = +0.328; e_prm.theta[1].x = 0.7382; e_prm.theta[1].y = 0.6518;
    //o_prm.k[0][0]  = +0.7700; o_prm.k[0][1]  = -0.5100; o_prm.k[1][0]  = +0.8444; o_prm.k[1][1]  = -0.6901;
    //o_prm.z[0][0]  = +0.0351; o_prm.z[0][1]  = +0.0238; o_prm.z[1][0]  = +0.1098; o_prm.z[1][1]  = -0.0302;
    //o_prm.xi[0].x  = +0.4605; o_prm.xi[0].y  = +0.4316; o_prm.xi[1].x  = +0.7556; o_prm.xi[1].y  = +0.6413;
    //o_prm.eta[0].x = +0.7516; o_prm.eta[0].y = +0.8028; o_prm.eta[1].x = +0.4595; o_prm.eta[1].y = +0.4302;

    // 0.008546
    //e_prm.q[0] = +0.314; e_prm.theta[0].x = 0.2845; e_prm.theta[0].y = 0.2845;
    //e_prm.q[1] = +0.328; e_prm.theta[1].x = 0.7382; e_prm.theta[1].y = 0.6518;
    //o_prm.k[0][0]  = +0.4062; o_prm.k[0][1]  = -0.1118; o_prm.k[1][0]  = +0.6411; o_prm.k[1][1]  = -0.5024;
    //o_prm.z[0][0]  = +0.0207; o_prm.z[0][1]  = -0.0022; o_prm.z[1][0]  = +0.0460; o_prm.z[1][1]  = -0.0277;
    //o_prm.xi[0].x  = +0.4985; o_prm.xi[0].y  = +0.5005; o_prm.xi[1].x  = +0.7629; o_prm.xi[1].y  = +0.6422;
    //o_prm.eta[0].x = +0.7180; o_prm.eta[0].y = +0.7701; o_prm.eta[1].x = +0.4996; o_prm.eta[1].y = +0.4999;

    // 0.016549 --------------------------------------------------------------------------------------------
    //e_prm.q[0] = +0.314; e_prm.theta[0].x = 0.2845; e_prm.theta[0].y = 0.2845;
    //e_prm.q[1] = +0.328; e_prm.theta[1].x = 0.7382; e_prm.theta[1].y = 0.6518;
    //o_prm.k[0][0]  = +0.0000; o_prm.k[0][1]  = +0.0000; o_prm.k[1][0]  = +0.0000; o_prm.k[1][1]  = +0.0000;
    //o_prm.z[0][0]  = +0.0000; o_prm.z[0][1]  = +0.0000; o_prm.z[1][0]  = +0.0000; o_prm.z[1][1]  = +0.0000;
    //o_prm.xi[0].x  = +0.3849; o_prm.xi[0].y  = +0.5442; o_prm.xi[1].x  = +0.7861; o_prm.xi[1].y  = +0.6785;
    //o_prm.eta[0].x = +0.6656; o_prm.eta[0].y = +0.7909; o_prm.eta[1].x = +0.4956; o_prm.eta[1].y = +0.3810;
    //-----------------------------------------------------------------------------------------------------

    // 0.012435 --------------------------------------------------------------------------------------------
    //e_prm.q[0] = +0.314; e_prm.theta[0].x = 0.2845; e_prm.theta[0].y = 0.2845;
    //e_prm.q[1] = +0.328; e_prm.theta[1].x = 0.7382; e_prm.theta[1].y = 0.6518;
    //o_prm.k[0][0]  = -1.1820; o_prm.k[0][1]  = -1.1250; o_prm.k[1][0]  = -1.1550; o_prm.k[1][1]  = -1.1310;
    //o_prm.z[0][0]  = +0.0125; o_prm.z[0][1]  = +0.0268; o_prm.z[1][0]  = +0.0359; o_prm.z[1][1]  = +0.0186;
    //o_prm.xi[0].x  = +0.3849; o_prm.xi[0].y  = +0.5442; o_prm.xi[1].x  = +0.7861; o_prm.xi[1].y  = +0.6785;
    //o_prm.eta[0].x = +0.6656; o_prm.eta[0].y = +0.7909; o_prm.eta[1].x = +0.4956; o_prm.eta[1].y = +0.3810;
    //-----------------------------------------------------------------------------------------------------

    // 0.009912 -------------------------------------------------------------------------------------------
    //    e_prm.q[0] = +0.314; e_prm.theta[0].x = 0.2845; e_prm.theta[0].y = 0.2845;
    //    e_prm.q[1] = +0.328; e_prm.theta[1].x = 0.7382; e_prm.theta[1].y = 0.6518;
    //    o_prm.k[0][0]  = -2.1866; o_prm.k[0][1]  = -2.1325; o_prm.k[1][0]  = -2.1321; o_prm.k[1][1]  = -2.1712;
    //    o_prm.z[0][0]  = +0.0359; o_prm.z[0][1]  = +0.0496; o_prm.z[1][0]  = +0.1320; o_prm.z[1][1]  = +0.1145;
    //    o_prm.xi[0].x  = +0.2732; o_prm.xi[0].y  = +0.6012; o_prm.xi[1].x  = +0.7578; o_prm.xi[1].y  = +0.7590;
    //    o_prm.eta[0].x = +0.6968; o_prm.eta[0].y = +0.6892; o_prm.eta[1].x = +0.3987; o_prm.eta[1].y = +0.3898;
    //o_prm.eta[0].x = +0.2500; o_prm.eta[0].y = +0.2500; o_prm.eta[1].x = +0.7500; o_prm.eta[1].y = +0.7500;
    //-----------------------------------------------------------------------------------------------------

    // Regularization parameters
    OptimizeParameter1H r_prm;
    r_prm.k.resize(e_prm.Nc, e_prm.No, 0.0);
    r_prm.z.resize(e_prm.Nc, e_prm.No, 0.0);
    r_prm.ksi.resize(e_prm.No);
    r_prm.eta.resize(e_prm.Nc);

    //r_prm.k[0][0]  = -0.0182; r_prm.k[0][1]  = -0.0125; r_prm.k[1][0]  = -0.0155; r_prm.k[1][1]  = -0.0131;
    //r_prm.z[0][0]  = -0.0262; r_prm.z[0][1]  = -0.0773; r_prm.z[1][0]  = -0.0570; r_prm.z[1][1]  = +0.0653;
    //r_prm.xi[0].x  = +0.3849; r_prm.xi[0].y  = +0.5442; r_prm.xi[1].x  = +0.7661; r_prm.xi[1].y  = +0.6785;
    //r_prm.eta[0].x = +0.6656; r_prm.eta[0].y = +0.7909; r_prm.eta[1].x = +0.4856; r_prm.eta[1].y = +0.3810;

    //r_prm.k[0][0]  = -1.1820; r_prm.k[0][1]  = -1.1250; r_prm.k[1][0]  = -1.1550; r_prm.k[1][1]  = -1.1310;
    //r_prm.k[0][0]  = -2.1820; r_prm.k[0][1]  = -2.1250; r_prm.k[1][0]  = -2.1550; r_prm.k[1][1]  = -2.1310;
    //r_prm.k[0][0]  = -3.1820; r_prm.k[0][1]  = -3.1250; r_prm.k[1][0]  = -3.1550; r_prm.k[1][1]  = -3.1310;
    //r_prm.k[0][0]  = +0.1820; r_prm.k[0][1]  = +0.1250; r_prm.k[1][0]  = +0.1550; r_prm.k[1][1]  = +0.1310;
    //r_prm.z[0][0]  = -0.0262; r_prm.z[0][1]  = -0.0773; r_prm.z[1][0]  = -0.0570; r_prm.z[1][1]  = +0.0653;
    //r_prm.xi[0].x  = +0.3849; r_prm.xi[0].y  = +0.5442; r_prm.xi[1].x  = +0.7661; r_prm.xi[1].y  = +0.6785;
    //r_prm.eta[0].x = +0.6656; r_prm.eta[0].y = +0.7909; r_prm.eta[1].x = +0.4856; r_prm.eta[1].y = +0.3810;

    //r_prm = o_prm;
    //o_prm = r_prm;

    // Penalty paramteres
    DoubleVector r; r << 0.0000 << 0.0000 << 0.000;
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
            prob.checkGradient1(prob);
            //prob.checkGradient2(prob);
            //return;

//            printf("fx: %.10f\n", prob.fx(x));
//            for (int n=-100; n<=100; n++)
//            {
//                x[3] = n*0.01;
//                printf("%f %f\n", n*0.01, prob.fx(x));
//            }
//
//            return;

            //            std::vector<DoubleMatrix> u;
            //            spif_vectorH u_info, p_info;
            //            prob.solveForwardIBVP(u, u_info, true);
            //            IPrinter::printSeperatorLine();
            //            prob.solveBackwardIBVP(u, p_info, true, u_info);
            //            return;

            //            spif_vectorH p_info;
            //            prob.solveBackwardIBVP(u, p_info, false, u_info);

            //            return;


            //FILE* cfile = fopen("e:/Google Drive/Problems/Problem2_H/Experiments/control2.txt", "w");
            //for (unsigned int ln=0; ln<=3000; ln++)
            //{
            //    fprintf(cfile, "%d,%.10f,%.10f\n", ln, prob.v(0, o_prm, e_prm, u_info, ln), prob.v(1, o_prm, e_prm, u_info, ln));
            //}
            //fclose(cfile);

            //            return;

            //prod_figure1(prob, ht, x);
            //return;
            //            printf("Fx: %10.6f\n", prob.fx(x)); return;
            //            IPrinter::printSeperatorLine();
            //            IPrinter::printMatrix(8,4,u[0]);
            //            prob.fx(x);
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

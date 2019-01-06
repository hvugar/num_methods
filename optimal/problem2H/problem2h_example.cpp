#include "problem2h_example.h"

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

void Problem2HNDirichlet::Main(int argc UNUSED_PARAM, char *argv[] UNUSED_PARAM)
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
    e_prm.lambda = +0.01;

    // Pulse influences
    e_prm.Ns = 2;
    e_prm.q.resize(e_prm.Ns);
    e_prm.theta.resize(e_prm.Ns);

    e_prm.q[0] = +0.514; e_prm.theta[0].x = 0.2500; e_prm.theta[0].y = 0.2500;
    e_prm.q[1] = +0.528; e_prm.theta[1].x = 0.7500; e_prm.theta[1].y = 0.7500;

    e_prm.No = 2;
    e_prm.Nc = 2;

    // Optimization parameters
    OptimizeParameterH o_prm;
    o_prm.k.resize(e_prm.Nc, e_prm.No, 0.0);
    o_prm.z.resize(e_prm.Nc, e_prm.No, 0.0);
    o_prm.xi.resize(e_prm.No);
    o_prm.eta.resize(e_prm.Nc);

    // Initial I:0.03533
    //o_prm.k[0][0]  = +0.0000; o_prm.k[0][1]  = +0.0000; o_prm.k[1][0]  = +0.0000; o_prm.k[1][1]  = +0.0000;
    //o_prm.z[0][0]  = +0.0000; o_prm.z[0][1]  = +0.0000; o_prm.z[1][0]  = +0.0000; o_prm.z[1][1]  = +0.0000;
    o_prm.k[0][0]  = -1.1820; o_prm.k[0][1]  = -1.1250; o_prm.k[1][0]  = -1.1550; o_prm.k[1][1]  = -1.1310;
    //o_prm.k[0][0]  = -2.1820; o_prm.k[0][1]  = -2.1250; o_prm.k[1][0]  = -2.1550; o_prm.k[1][1]  = -2.1310;
    o_prm.z[0][0]  = +0.0125; o_prm.z[0][1]  = +0.0268; o_prm.z[1][0]  = +0.0359; o_prm.z[1][1]  = +0.0186;
    o_prm.xi[0].x  = +0.3849; o_prm.xi[0].y  = +0.5442; o_prm.xi[1].x  = +0.7861; o_prm.xi[1].y  = +0.6785;
    o_prm.eta[0].x = +0.6656; o_prm.eta[0].y = +0.7909; o_prm.eta[1].x = +0.4956; o_prm.eta[1].y = +0.3810;
    //o_prm.xi[0].x  = +0.1849; o_prm.xi[0].y  = +0.2442; o_prm.xi[1].x  = +0.4861; o_prm.xi[1].y  = +0.7785;
    //o_prm.eta[0].x = +0.3656; o_prm.eta[0].y = +0.3909; o_prm.eta[1].x = +0.9156; o_prm.eta[1].y = +0.9410;
    // Optimal I:0.00054           --------------------------------------------------------------------------
    //o_prm.k[0][0]  = -2.1866; o_prm.k[0][1]  = -2.1325; o_prm.k[1][0]  = -2.1321; o_prm.k[1][1]  = -2.1712;
    //o_prm.z[0][0]  = +0.0359; o_prm.z[0][1]  = +0.0496; o_prm.z[1][0]  = +0.1320; o_prm.z[1][1]  = +0.1145;
    //o_prm.xi[0].x  = +0.2732; o_prm.xi[0].y  = +0.6012; o_prm.xi[1].x  = +0.7578; o_prm.xi[1].y  = +0.7590;
    //o_prm.eta[0].x = +0.6968; o_prm.eta[0].y = +0.6892; o_prm.eta[1].x = +0.3987; o_prm.eta[1].y = +0.3898;
    //o_prm.eta[0].x = +0.2500; o_prm.eta[0].y = +0.2500; o_prm.eta[1].x = +0.7500; o_prm.eta[1].y = +0.7500;
    // ------------------------------------------------------------------------------------------------------
    // Optimal I:0.00093           --------------------------------------------------------------------------
    //o_prm.k[0][0]  = -1.2826; o_prm.k[0][1]  = -1.0983; o_prm.k[1][0]  = -1.4301; o_prm.k[1][1]  = -1.1559;
    //o_prm.z[0][0]  = -0.1559; o_prm.z[0][1]  = -0.1369; o_prm.z[1][0]  = -0.0990; o_prm.z[1][1]  = -0.1317;
    //o_prm.xi[0].x  = +0.1391; o_prm.xi[0].y  = +0.8715; o_prm.xi[1].x  = +0.7716; o_prm.xi[1].y  = +0.7722;
    //o_prm.eta[0].x = +0.6358; o_prm.eta[0].y = +0.5878; o_prm.eta[1].x = +0.3662; o_prm.eta[1].y = +0.3292;
    // ------------------------------------------------------------------------------------------------------
    // Optimal I:0.00046           --------------------------------------------------------------------------
    //o_prm.k[0][0]  = -2.2155; o_prm.k[0][1]  = -2.1080; o_prm.k[1][0]  = -2.0897; o_prm.k[1][1]  = -2.2124;
    //o_prm.z[0][0]  = +0.0359; o_prm.z[0][1]  = +0.0496; o_prm.z[1][0]  = +0.1208; o_prm.z[1][1]  = -0.1028;
    //o_prm.xi[0].x  = +0.3943; o_prm.xi[0].y  = +0.7138; o_prm.xi[1].x  = +0.7461; o_prm.xi[1].y  = +0.6914;
    //o_prm.eta[0].x = +0.6914; o_prm.eta[0].y = +0.6823; o_prm.eta[1].x = +0.3760; o_prm.eta[1].y = +0.3700;
    // ------------------------------------------------------------------------------------------------------

    // Regularization parameters
    OptimizeParameterH r_prm;
    r_prm.k.resize(e_prm.Nc, e_prm.No, 0.0);
    r_prm.z.resize(e_prm.Nc, e_prm.No, 0.0);
    r_prm.xi.resize(e_prm.No);
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

    // Grid parameters
    double hx = 0.010; int Nx = 100;
    double hy = 0.010; int Ny = 100;
    //double ht = 0.010; int Nt = 500;
    double ht = 0.005; int Nt = 1000;

    Dimension time(ht, 0, Nt);
    Dimension dimx(hx, 0, Nx);
    Dimension dimy(hy, 0, Ny);

    // Penalty paramteres
    DoubleVector r; r << 0.0000;// << 1.0000 << 20.000 << 50.000;
    // Regularization coefficients
    DoubleVector e; e << 0.0000;// << 0.1000 << 0.0100 << 0.0000;

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
        prob.vmin.resize(e_prm.Nc, -0.2);
        prob.vmax.resize(e_prm.Nc, +0.2);
        prob.LD = 100;
        prob.noise = 0.00;

        prob.regEpsilon = e[i];
        prob.r = r[i];

        if (i==0)
        {
            prob.PrmToVector(o_prm, x);
            prob.checkGradient1(prob);

            //std::vector<DoubleMatrix> u;
            //spif_vectorH u_info;
            //prob.solveForwardIBVP(u, u_info, false);
            //return;


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

        ConjugateGradient g;
        //SteepestDescentGradient g;
        g.setFunction(&prob);
        g.setGradient(&prob);
        g.setPrinter(&prob);
        //g.setProjection(&prob);
        g.setProjection(new ProjectionEx1);
        //g.setGradientNormalizer(&prob);
        g.setOptimalityTolerance(0.0001);
        g.setFunctionTolerance(0.0001);
        g.setStepTolerance(0.0001);
        g.setR1MinimizeEpsilon(0.1, 0.01);
        g.setMaxIterations(50);
        g.setNormalize(false);
        g.showExitMessage(true);
        //prob.gm = &g;

        g.calculate(x);

        IPrinter::printSeperatorLine(nullptr, '=');
    }
}

void prod_figure1(Problem2HNDirichlet &prob, double ht, const DoubleVector &x)
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

void example1()
{
    // Equation parameters
    EquationParameterH e_prm;
    e_prm.a = 1.0;
    e_prm.lambda = +0.001;

    // Pulse influences
    e_prm.Ns = 2;
    e_prm.q.resize(e_prm.Ns);
    e_prm.theta.resize(e_prm.Ns);

    //    e_prm.q[0] = +0.145; e_prm.theta[0].x = 0.2500; e_prm.theta[0].y = 0.7200;
    //    e_prm.q[1] = +0.157; e_prm.theta[1].x = 0.6400; e_prm.theta[1].y = 0.2700;
    e_prm.q[0] = +0.145; e_prm.theta[0].x = 0.2500; e_prm.theta[0].y = 0.7200;
    e_prm.q[1] = +0.157; e_prm.theta[1].x = 0.6400; e_prm.theta[1].y = 0.2700;

    e_prm.No = 2;
    e_prm.Nc = 2;

    // Optimization parameters
    OptimizeParameterH o_prm;
    o_prm.k.resize(e_prm.Nc, e_prm.No, 0.0);
    o_prm.z.resize(e_prm.Nc, e_prm.No, 0.0);
    o_prm.xi.resize(e_prm.No);
    o_prm.eta.resize(e_prm.Nc);

    //o_prm.k[0][0]  = -2.1062; o_prm.k[0][1]  = -2.1038; o_prm.k[1][0]  = -2.1031; o_prm.k[1][1]  = -2.1052;
    //o_prm.z[0][0]  = -0.0245; o_prm.z[0][1]  = -0.0784; o_prm.z[1][0]  = -0.0587; o_prm.z[1][1]  = -0.0641;
    //o_prm.xi[0].x  = +0.1527; o_prm.xi[0].y  = +0.8412; o_prm.xi[1].x  = +0.7412; o_prm.xi[1].y  = +0.2483;
    //o_prm.eta[0].x = +0.3254; o_prm.eta[0].y = +0.3654; o_prm.eta[1].x = +0.9462; o_prm.eta[1].y = +0.4966;

    o_prm.k[0][0]  = -0.1262; o_prm.k[0][1]  = -0.4038; o_prm.k[1][0]  = -1.7431; o_prm.k[1][1]  = -2.8052;
    o_prm.z[0][0]  = -2.0245; o_prm.z[0][1]  = +3.0784; o_prm.z[1][0]  = -5.0587; o_prm.z[1][1]  = +8.0641;
    o_prm.xi[0].x  = +0.6527; o_prm.xi[0].y  = +0.8412; o_prm.xi[1].x  = +0.7412; o_prm.xi[1].y  = +0.2483;
    o_prm.eta[0].x = +0.3254; o_prm.eta[0].y = +0.3654; o_prm.eta[1].x = +0.9462; o_prm.eta[1].y = +0.4966;

    //o_prm.k[0][0]  = -0.1262; o_prm.k[0][1]  = -0.4038; o_prm.k[1][0]  = -0.7431; o_prm.k[1][1]  = -0.8052;
    //o_prm.z[0][0]  = -0.0245; o_prm.z[0][1]  = +0.0784; o_prm.z[1][0]  = -0.0587; o_prm.z[1][1]  = +0.0641;
    //o_prm.xi[0].x  = +0.1527; o_prm.xi[0].y  = +0.8412; o_prm.xi[1].x  = +0.7412; o_prm.xi[1].y  = +0.2483;
    //o_prm.eta[0].x = +0.3254; o_prm.eta[0].y = +0.3654; o_prm.eta[1].x = +0.9462; o_prm.eta[1].y = +0.5966;


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

    o_prm = r_prm;

    // Grid parameters
    double hx = 0.010; int Nx = 100;
    double hy = 0.010; int Ny = 100;
    double ht = 0.010; int Nt = 500;

    Dimension time(ht, 0, Nt);
    Dimension dimx(hx, 0, Nx);
    Dimension dimy(hy, 0, Ny);

    // Penalty paramteres
    DoubleVector r; r << 1.0000 << 10.000 << 20.000 << 50.000;
    // Regularization coefficients
    DoubleVector e; e << 1.0000 << 0.0000 << 0.0000 << 0.0000;
    //DoubleVector e; e << 1.00 << 0.10 << 0.010 << 0.0010;

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
        prob.vmin.resize(e_prm.Nc, -0.05);
        prob.vmax.resize(e_prm.Nc, +0.05);
        prob.LD = 50;
        prob.noise = 0.00;

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

            //            for (int i=0; i<10000; i++)
            //            {
            //                prob.setTimeDimension(Dimension(ht, 0, i));
            //                printf("ln:%d fx:%8.6f\n", i, prob.fx(x));
            //                //IPrinter::printMatrix(8,4,u[0]);
            //            }

            //            IPrinter::printMatrix(8,4,u[0]);
            //            return;

        }

        ConjugateGradient g;
        //SteepestDescentGradient g;
        g.setFunction(&prob);
        g.setGradient(&prob);
        g.setPrinter(&prob);
        g.setProjection(&prob);
        g.setProjection(new ProjectionEx1);
        //g.setGradientNormalizer(&prob);
        g.setOptimalityTolerance(0.0001);
        g.setFunctionTolerance(0.0001);
        g.setStepTolerance(0.0001);
        g.setR1MinimizeEpsilon(0.1, 0.01);
        g.setMaxIterations(20);
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

    //e_prm.q[0] = +0.645; e_prm.theta[0].x = 0.2500; e_prm.theta[0].y = 0.7200;
    //e_prm.q[1] = +0.657; e_prm.theta[1].x = 0.5400; e_prm.theta[1].y = 0.2700;
    e_prm.q[0] = +0.245; e_prm.theta[0].x = 0.3500; e_prm.theta[0].y = 0.3500;
    e_prm.q[1] = +0.257; e_prm.theta[1].x = 0.6500; e_prm.theta[1].y = 0.6500;

    e_prm.No = 2;
    e_prm.Nc = 2;

    // Optimization parameters
    OptimizeParameterH o_prm;
    o_prm.k.resize(e_prm.Nc, e_prm.No, 0.0);
    o_prm.z.resize(e_prm.Nc, e_prm.No, 0.0);
    o_prm.xi.resize(e_prm.No);
    o_prm.eta.resize(e_prm.Nc);

    //o_prm.k[0][0]  = -1.0345; o_prm.k[0][1]  = -2.5401; o_prm.k[1][0]  = -1.1431; o_prm.k[1][1]  = -1.0984;
    //o_prm.k[0][0]  = -0.3450; o_prm.k[0][1]  = -0.5410; o_prm.k[1][0]  = -0.4310; o_prm.k[1][1]  = -0.9840;
    o_prm.k[0][0]  = -0.0345*(0.9); o_prm.k[0][1]  = -0.0541*(1.1); o_prm.k[1][0]  = -0.0431*(1.1); o_prm.k[1][1]  = -0.0984*(0.9);

    //o_prm.z[0][0]  = +3.1245; o_prm.z[0][1]  = +1.8532; o_prm.z[1][0]  = -2.4512; o_prm.z[1][1]  = +2.5421;
    //o_prm.z[0][0]  = +0.1245; o_prm.z[0][1]  = +0.0325; o_prm.z[1][0]  = -0.1452; o_prm.z[1][1]  = +0.2154;
    o_prm.z[0][0]  = +0.0124*(1.1); o_prm.z[0][1]  = +0.0325*(0.9); o_prm.z[1][0]  = -0.0452*(1.1); o_prm.z[1][1]  = +0.0154*(0.9);

    o_prm.xi[0].x  = +0.7524+(+0.01); o_prm.xi[0].y  = +0.4828+(-0.01); o_prm.xi[1].x  = +0.1274+(+0.01); o_prm.xi[1].y  = +0.8234+(-0.01);
    //o_prm.eta[0].x = +0.8635; o_prm.eta[0].y = +0.3654; o_prm.eta[1].x = +0.2496; o_prm.eta[1].y = +0.6536;
    //o_prm.eta[0].x = +0.3500; o_prm.eta[0].y = +0.3500; o_prm.eta[1].x = +0.6500; o_prm.eta[1].y = +0.6500;
    o_prm.eta[0].x = +0.3500+(+0.01); o_prm.eta[0].y = +0.6500+(-0.01); o_prm.eta[1].x = +0.6500+(+0.01); o_prm.eta[1].y = +0.3500+(-0.01);
    //o_prm.xi[0].x  = +0.2754; o_prm.xi[0].y  = +0.2438; o_prm.xi[1].x  = +0.4271; o_prm.xi[1].y  = +0.3824;
    //o_prm.eta[0].x = +0.3865; o_prm.eta[0].y = +0.4365; o_prm.eta[1].x = +0.6294; o_prm.eta[1].y = +0.5466;

    o_prm.k[0][0]  = +0.1200; o_prm.k[0][1]  = +0.2400; o_prm.k[1][0]  = +0.4500; o_prm.k[1][1]  = +0.1800;
    o_prm.z[0][0]  = +0.5000; o_prm.z[0][1]  = -0.4000; o_prm.z[1][0]  = +0.7000; o_prm.z[1][1]  = +0.5000;
    o_prm.xi[0].x  = +0.4274; o_prm.xi[0].y  = +0.6735; o_prm.xi[1].x  = +0.6710; o_prm.xi[1].y  = +0.3851;
    o_prm.eta[0].x = +0.5174; o_prm.eta[0].y = +0.7635; o_prm.eta[1].x = +0.5570; o_prm.eta[1].y = +0.4751;

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

    r_prm.k[0][0]  = -0.3839; r_prm.k[0][1]  = -0.7371; r_prm.k[1][0]  = -0.5532; r_prm.k[1][1]  = -1.3051;
    r_prm.z[0][0]  = +0.0238; r_prm.z[0][1]  = +0.0467; r_prm.z[1][0]  = -0.0352; r_prm.z[1][1]  = +0.0366;
    r_prm.xi[0].x  = +0.5203; r_prm.xi[0].y  = +0.2636; r_prm.xi[1].x  = +0.2521; r_prm.xi[1].y  = +0.7225;
    r_prm.eta[0].x = +0.5911; r_prm.eta[0].y = +0.4344; r_prm.eta[1].x = +0.4637; r_prm.eta[1].y = +0.5160;

    r_prm.k[0][0]  = -0.6010; r_prm.k[0][1]  = -0.8657; r_prm.k[1][0]  = -0.8284; r_prm.k[1][1]  = -1.4035;
    r_prm.z[0][0]  = -0.0108; r_prm.z[0][1]  = -0.0128; r_prm.z[1][0]  = -0.0671; r_prm.z[1][1]  = -0.0328;
    r_prm.xi[0].x  = +0.5117; r_prm.xi[0].y  = +0.2698; r_prm.xi[1].x  = +0.2562; r_prm.xi[1].y  = +0.7156;
    r_prm.eta[0].x = +0.5353; r_prm.eta[0].y = +0.4041; r_prm.eta[1].x = +0.4411; r_prm.eta[1].y = +0.5339;

    //k: -0.3839  -0.7371  -0.5532  -1.3051
    //z:  0.0238   0.0467  -0.0352   0.0366
    //o:  0.5203   0.2636   0.2521   0.7225
    //c:  0.5911   0.4344   0.4637   0.5160

    //r_prm = o_prm;
    //o_prm = r_prm;

    // Grid parameters
    double hx = 0.010; int Nx = 100;
    double hy = 0.010; int Ny = 100;
    double ht = 0.010; int Nt = 500;

    Dimension time(ht, 0, Nt);
    Dimension dimx(hx, 0, Nx);
    Dimension dimy(hy, 0, Ny);

    // Penalty paramteres
    //DoubleVector r; r << 0.10 << 1.0 << 10.0 << 100.00;
    DoubleVector r; r << 1.00;// << 10.00 << 100.00 << 0.00;
    // Regularization coefficients
    //DoubleVector e; e << 1.00 << 0.10 << 0.010 << 0.00100;
    DoubleVector e; e << 0.00;// << 0.00 << 0.000 << 0.000;

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
        prob.vmin.resize(e_prm.Nc, -0.05);
        prob.vmax.resize(e_prm.Nc, +0.05);
        prob.LD = 50;

        prob.regEpsilon = e[i];
        prob.r = r[i];
        if (i==0)
        {
            prob.PrmToVector(o_prm, x);
            //prob.checkGradient1(prob);
            IPrinter::printSeperatorLine();

            std::vector<DoubleMatrix> u;
            spif_vectorH u_info;
            prob.solveForwardIBVP(u, u_info, false);

            //            FILE *file = fopen("e:/data.txt", "w");
            //            for (int i=0; i<10000; i++)
            //            {
            //                prob.setTimeDimension(Dimension(ht, 0, i));
            //                fprintf(file, "ln:%d fx:%8.6f\n", i, prob.fx(x));
            //                fprintf(stdout, "ln:%d fx:%8.6f\n", i, prob.fx(x));
            //                fflush(file);
            //                IPrinter::printMatrix(8,4,u[0]);
            //            }
            //            fclose(file);
            //            IPrinter::printMatrix(8,4,u[0]);
            //            return;
        }

        //ConjugateGradient g; g.setAlgorithm(ConjugateGradient::Algorithm::FLETCHER_REEVES);
        SteepestDescentGradient g;
        g.setFunction(&prob);
        g.setGradient(&prob);
        g.setPrinter(&prob);
        g.setProjection(&prob);
        //g.setGradientNormalizer(&prob);
        g.setProjection(new ProjectionEx1);
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

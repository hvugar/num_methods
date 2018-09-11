#include "problem2hN.h"

#define SAVE_TO_IMG1s
#define EXAMPLE4_SAMPLE_2

void Problem2HNDirichlet::Main(int argc UNUSED_PARAM, char *argv[] UNUSED_PARAM)
{
    //example1();
    //example2();
    //example3();
    example4();
}

void example4()
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

    //y10
#ifdef EXAMPLE4_SAMPLE_1
    o_prm.k[0][0]  = +1.1200; o_prm.k[0][1]  = +1.2400; o_prm.k[1][0]  = +1.4500; o_prm.k[1][1]  = +1.1800;
    o_prm.z[0][0]  = +0.5000; o_prm.z[0][1]  = -0.4000; o_prm.z[1][0]  = +0.7000; o_prm.z[1][1]  = +0.5000;
    o_prm.xi[0].x  = +0.4274; o_prm.xi[0].y  = +0.6735; o_prm.xi[1].x  = +0.6710; o_prm.xi[1].y  = +0.3851;
    o_prm.eta[0].x = +0.5174; o_prm.eta[0].y = +0.7635; o_prm.eta[1].x = +0.5570; o_prm.eta[1].y = +0.4751;

    r_prm.k[0][0]  = +0.4639; r_prm.k[0][1]  = -0.0136; r_prm.k[1][0]  = +0.1977; r_prm.k[1][1]  = -0.5896;
    r_prm.z[0][0]  = +0.3014; r_prm.z[0][1]  = -0.6160; r_prm.z[1][0]  = -0.1914; r_prm.z[1][1]  = -0.2933;
    r_prm.xi[0].x  = +0.4679; r_prm.xi[0].y  = +0.5770; r_prm.xi[1].x  = +0.7140; r_prm.xi[1].y  = +0.2614;
    r_prm.eta[0].x = +0.5579; r_prm.eta[0].y = +0.8282; r_prm.eta[1].x = +0.8040; r_prm.eta[1].y = +0.7535;
#endif

    //y20
#ifdef EXAMPLE4_SAMPLE_2
    o_prm.k[0][0]  = -2.6400; o_prm.k[0][1]  = +3.7400; o_prm.k[1][0]  = -2.1800; o_prm.k[1][1]  = -2.0700;
    o_prm.z[0][0]  = -0.9500; o_prm.z[0][1]  = +0.8500; o_prm.z[1][0]  = -0.1400; o_prm.z[1][1]  = -0.4500;
    o_prm.xi[0].x  = +0.1486; o_prm.xi[0].y  = +0.1284; o_prm.xi[1].x  = +0.7525; o_prm.xi[1].y  = +0.7920;
    o_prm.eta[0].x = +0.8512; o_prm.eta[0].y = +0.3245; o_prm.eta[1].x = +0.2854; o_prm.eta[1].y = +0.6515;

    r_prm.k[0][0]  = -0.5520; r_prm.k[0][1]  = +0.4493; r_prm.k[1][0]  = -0.9090; r_prm.k[1][1]  = -1.2765;
    r_prm.z[0][0]  = -1.8472; r_prm.z[0][1]  = +1.2049; r_prm.z[1][0]  = -0.1053; r_prm.z[1][1]  = -0.3725;
    r_prm.xi[0].x  = +0.0500; r_prm.xi[0].y  = +0.0500; r_prm.xi[1].x  = +0.0500; r_prm.xi[1].y  = +0.9500;
    r_prm.eta[0].x = +0.6062; r_prm.eta[0].y = +0.6485; r_prm.eta[1].x = +0.1481; r_prm.eta[1].y = +0.2337;
#endif


    //    o_prm.k[0][0]  = +1.0667; o_prm.k[0][1]  = +1.1309; o_prm.k[1][0]  = +1.0974; o_prm.k[1][1]  = +0.8603;
    //    o_prm.z[0][0]  = +0.4870; o_prm.z[0][1]  = -0.4140; o_prm.z[1][0]  = +0.3969; o_prm.z[1][1]  = +0.2696;
    //    o_prm.xi[0].x  = +0.0500; o_prm.xi[0].y  = +0.6306; o_prm.xi[1].x  = +0.9500; o_prm.xi[1].y  = +0.5327;
    //    o_prm.eta[0].x = +0.8600; o_prm.eta[0].y = +0.9500; o_prm.eta[1].x = +0.5128; o_prm.eta[1].y = +0.7914;

    //    o_prm.k[0][0]  = +1.0667; o_prm.k[0][1]  = +1.1309; o_prm.k[1][0]  = +1.0974; o_prm.k[1][1]  = +0.8603;
    //    o_prm.z[0][0]  = +0.4870; o_prm.z[0][1]  = -0.4140; o_prm.z[1][0]  = +0.3969; o_prm.z[1][1]  = +0.2696;
    //    o_prm.xi[0].x  = +0.0500; o_prm.xi[0].y  = +0.6306; o_prm.xi[1].x  = +0.9500; o_prm.xi[1].y  = +0.5327;
    //    o_prm.eta[0].x = +0.8600; o_prm.eta[0].y = +0.9500; o_prm.eta[1].x = +0.5128; o_prm.eta[1].y = +0.7914;

    //k: -0.4154   0.3307  -0.8037  -1.0511 z: -1.8620   1.2195  -0.0170  -0.2481 o:  0.0500   0.0500   0.0500   0.9500 c:  0.5822   0.6233   0.1400   0.1541
    //k: -0.0012   0.0042  -0.0007   0.0004 z:  0.0000  -0.0000   0.0006   0.0007 o:  0.0216   0.0217   0.0103  -0.0187 c:  0.0001  -0.0002   0.0014   0.0011
    //k: -0.0308   0.1121  -0.0186   0.0111 z:  0.0007  -0.0005   0.0151   0.0197 o:  0.5748   0.5771   0.2729  -0.4952 c:  0.0014  -0.0066   0.0382   0.0281

//    o_prm = r_prm;

    // Grid parameters
    double hx = 0.010; unsigned int Nx = 100;
    double hy = 0.010; unsigned int Ny = 100;
    double ht = 0.010; unsigned int Nt = 500;

    Dimension time(ht, 0, Nt);
    Dimension dimx(hx, 0, Nx);
    Dimension dimy(hy, 0, Ny);

    // Penalty paramteres
    DoubleVector r; r << 1.0000 << 10.000 << 100.000 << 1000.0;// << 20.000 << 50.000 << 100.00;
    // Regularization coefficients
    DoubleVector e; e << 1.0000 << 1.0000 << 1.00000 << 0.0000;// << 0.0000 << 0.0000 << 0.0000;

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
        prob.LD = 10;

        if (i==0)
        {
            prob.PrmToVector(o_prm, x);
            prob.checkGradient(prob);
            IPrinter::printSeperatorLine();
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
        g.setR1MinimizeEpsilon(1.0, 0.001);
        g.setNormalize(true);
        g.showEndMessage(true);

        g.calculate(x);

        IPrinter::printSeperatorLine(NULL, '=');
    }
}

void Problem2HNDirichlet::experimentInfo(const Problem2HNDirichlet &prob) const
{
    printf("Experiment #: %d\n", (unsigned int) time(NULL));

    printf("Grid info: time: %.6f x: %.6f y: %.6f\n", prob.timeDimension().step(),
           prob.spaceDimension(Dimension::DimensionX).step(),
           prob.spaceDimension(Dimension::DimensionY).step());

    //printf("Optimization epsilon: %f %f %f\n");
}

void Problem2HNDirichlet::f_layerInfo(const DoubleMatrix &u UNUSED_PARAM, unsigned int ln UNUSED_PARAM) const
{
#ifdef SAVE_TO_IMG
    double min = u.min();
    double max = u.max();
    //    if (MIN>min) MIN = min;
    //    if (MAX<max) MAX = max;

    //    double norm = 0.0;
    //    for (unsigned int m=0; m<u.rows(); m++)
    //    {
    //        for (unsigned int n=0; n<u.cols(); n++)
    //        {
    //            norm += u[m][n]*u[m][n];
    //        }
    //    }

    QPixmap pic;
    visualizeMatrixHeat(u, min, max, pic);
    pic.save("images/f/500/pic"+QString("%1").arg(ln)+".png", "PNG");
    printf("Layer: %d min: %f max: %f min: %f max: %f norm: %f\n", ln, min, max, min, max, fabs(max-min));
#endif
}

void Problem2HNDirichlet::initParameters(EquationParameter &e_prm, OptimizeParameter &o_prm, OptimizeParameter &o_prm0)
{
    e_prm.a = 1.0;
    e_prm.lambda = +0.01;

    {
        e_prm.Ns = 3;
        e_prm.q.resize(e_prm.Ns);
        e_prm.theta.resize(e_prm.Ns);

        e_prm.q[0] = 2.2; e_prm.theta[0].x = 0.5000; e_prm.theta[0].y = 0.5000;
        e_prm.q[1] = 3.3; e_prm.theta[1].x = 0.2000; e_prm.theta[1].y = 0.2000;
        e_prm.q[2] = 5.5; e_prm.theta[2].x = 0.8000; e_prm.theta[2].y = 0.8000;

        //e_prm.q[0] = -0.75; e_prm.theta[0].x = 0.5000; e_prm.theta[0].y = 0.5000;
        //e_prm.q[1] = -0.31; e_prm.theta[1].x = 0.2500; e_prm.theta[1].y = 0.2500;
        //e_prm.q[2] = -0.25; e_prm.theta[2].x = 0.2500; e_prm.theta[2].y = 0.7500;
        //e_prm.q[3] = -0.56; e_prm.theta[3].x = 0.7500; e_prm.theta[3].y = 0.7500;
        //e_prm.q[4] = -0.15; e_prm.theta[4].x = 0.7500; e_prm.theta[4].y = 0.2500;
    }

    {
        //e_prm.Nc = 2;
        //e_prm.No = 2;

        //o_prm.k.resize(e_prm.Nc, e_prm.No, 0.0);
        //o_prm.z.resize(e_prm.Nc, e_prm.No, 0.0);
        //o_prm.xi.resize(e_prm.No);
        //o_prm.eta.resize(e_prm.Nc);

        //o_prm.k[0][0] = -0.9186; o_prm.k[0][1] = -0.9051; o_prm.k[1][0] = -1.0169; o_prm.k[1][1] = -1.0221;
        //o_prm.z[0][0] = +0.3232; o_prm.z[0][1] = +0.2617; o_prm.z[1][0] = +0.2469; o_prm.z[1][1] = +0.0996;
        //o_prm.xi[0].x = 0.2078; o_prm.xi[0].y = 0.2080; o_prm.xi[1].x = 0.8455; o_prm.xi[1].y = 0.8491;
        //o_prm.eta[0].x = 0.7623; o_prm.eta[0].y = 0.2708; o_prm.eta[1].x = 0.2636; o_prm.eta[1].y = 0.7636;

        //o_prm.k[0][0] = -1.7000; o_prm.k[0][1] = -1.3000; o_prm.k[1][0] = -1.6000; o_prm.k[1][1] = -1.5000;
        //o_prm.z[0][0] = +0.0000; o_prm.z[0][1] = +0.0000; o_prm.z[1][0] = +0.0000; o_prm.z[1][1] = +0.0000;
        //o_prm.xi[0].x = +0.4000; o_prm.xi[0].y = +0.6000; o_prm.xi[1].x = 0.9000; o_prm.xi[1].y = 0.2000;
        //o_prm.eta[0].x = +0.2000; o_prm.eta[0].y = +0.4500; o_prm.eta[1].x = 0.6500; o_prm.eta[1].y = 0.7500;

        //o_prm0 = o_prm;
    }

    {
        //e_prm.No = 3;
        //e_prm.Nc = 2;

        //o_prm.k.resize(e_prm.Nc, e_prm.No, 0.0);
        //o_prm.z.resize(e_prm.Nc, e_prm.No, 0.0);
        //o_prm.xi.resize(e_prm.No);
        //o_prm.eta.resize(e_prm.Nc);

        //o_prm.k[0][0] = -0.1200; o_prm.k[0][1] = -0.2400; o_prm.k[0][2] = -2.2400; o_prm.k[1][0] = -0.4500; o_prm.k[1][1] = -0.1800; o_prm.k[1][2] = -2.1800;
        //o_prm.z[0][0] = +0.5000; o_prm.z[0][1] = +0.4000; o_prm.z[0][2] = +0.4000; o_prm.z[1][0] = +0.7000; o_prm.z[1][1] = +0.5000; o_prm.z[1][2] = +0.5000;
        //o_prm.xi[0].x = 0.5000; o_prm.xi[0].y = 0.6000; o_prm.xi[1].x = 0.7000; o_prm.xi[1].y = 0.2000; o_prm.xi[2].x = 0.5000; o_prm.xi[2].y = 0.5000;
        //o_prm.eta[0].x = 0.2000; o_prm.eta[0].y = 0.7000; o_prm.eta[1].x = 0.8000; o_prm.eta[1].y = 0.3000;
    }

    {
        /**/
        e_prm.No = 2;
        e_prm.Nc = 2;

        o_prm.k.resize(e_prm.Nc, e_prm.No, 0.0);
        o_prm.z.resize(e_prm.Nc, e_prm.No, 0.0);
        o_prm.xi.resize(e_prm.No);
        o_prm.eta.resize(e_prm.Nc);

        o_prm.xi[0].x = 0.3000; o_prm.xi[0].y = 0.8000;
        o_prm.xi[1].x = 0.6000; o_prm.xi[1].y = 0.4000;

        o_prm.eta[0].x = 0.5000; o_prm.eta[0].y = 0.7000;
        o_prm.eta[1].x = 0.7000; o_prm.eta[1].y = 0.3000;

        o_prm.k[0][0] = +1.1200; o_prm.k[0][1] = +1.2400;
        o_prm.k[1][0] = +2.4500; o_prm.k[1][1] = +2.1800;

        o_prm.z[0][0] = +0.5000; o_prm.z[0][1] = +0.4000;
        o_prm.z[1][0] = +0.7000; o_prm.z[1][1] = +0.5000;

        o_prm0 = o_prm;
        /**/

        /*
        o_prm0.k.resize(e_prm.Nc, e_prm.No, 0.0);
        o_prm0.z.resize(e_prm.Nc, e_prm.No, 0.0);

        o_prm0.xi.resize(e_prm.No);
        o_prm0.xi[0].x = 0.4000; o_prm0.xi[0].y = 0.4000;
        o_prm0.xi[1].x = 0.6000; o_prm0.xi[1].y = 0.6000;

        o_prm0.eta.resize(e_prm.Nc);
        o_prm0.eta[0].x = 0.3000; o_prm0.eta[0].y = 0.4000;
        o_prm0.eta[1].x = 0.7000; o_prm0.eta[1].y = 0.3000;

        o_prm0.k[0][0] = -0.1200; o_prm0.k[0][1] = -0.2400;
        o_prm0.k[1][0] = -0.4500; o_prm0.k[1][1] = -0.1800;

        o_prm0.z[0][0] = +5.5000; o_prm0.z[0][1] = +4.4000;
        o_prm0.z[1][0] = +4.7000; o_prm0.z[1][1] = +5.5000;
        */
    }

    {
        //I[ 40]: 0.024768 0.300001 0.324769 R:0.00 e:0.000
        //k: -0.9095  -0.8781  -0.9951  -1.0015 z:  0.3378   0.2736   0.2531   0.1034   o:0.1995 0.2099 0.8449 0.8503   c:0.7553 0.2647 0.2676 0.7584
        //k:  0.0332   0.0427   0.0295   0.0409 z:  0.0167   0.0161   0.0137   0.0138   o:-0.1201 0.1606 -0.0308 0.0131   c:-0.1623 -0.1512 0.1166 -0.0317
        //0.004568 1
    }

    {
        /*
        OptimizeParameter o_prm0;
        o_prm0.k.resize(e_prm.Nc, e_prm.No, 0.0);
        o_prm0.z.resize(e_prm.Nc, e_prm.No, 0.0);

        o_prm0.k[0][0] = -2.6388; o_prm0.k[0][1] = -2.4754; o_prm0.k[1][0] = -2.8992; o_prm0.k[1][1] = -2.7987;
        o_prm0.z[0][0] = +0.0991; o_prm0.z[0][1] = +0.1211; o_prm0.z[1][0] = -0.0860; o_prm0.z[1][1] = +0.3928;

        o_prm0.xi.resize(e_prm.No);
        o_prm0.xi[0].x = 0.2402; o_prm0.xi[0].y = 0.1762; o_prm0.xi[1].x = 0.5736; o_prm0.xi[1].y = 0.3214;

        o_prm0.eta.resize(e_prm.Nc);
        o_prm0.eta[0].x = 0.2743; o_prm0.eta[0].y = 0.7485; o_prm0.eta[1].x = 0.7657; o_prm0.eta[1].y = 0.2613;
        */
    }

    {
        /*
        e_prm.No = 1;
        e_prm.Nc = 1;

        o_prm.k.resize(e_prm.Nc, e_prm.No, 0.0);
        o_prm.z.resize(e_prm.Nc, e_prm.No, 0.0);
        o_prm.xi.resize(e_prm.No);
        o_prm.eta.resize(e_prm.Nc);

        o_prm.xi[0].x = 0.3000; o_prm.xi[0].y = 0.5000;

        o_prm.eta[0].x = 0.7000; o_prm.eta[0].y = 0.5000;

        o_prm.k[0][0] = -2.1200;

        o_prm.z[0][0] = +1.5000;

        o_prm0 = o_prm;
        */
    }

    {
        /*
        e_prm.Nc = 1;
        e_prm.No = 2;

        o_prm.k.resize(e_prm.Nc, e_prm.No, 0.0);
        o_prm.z.resize(e_prm.Nc, e_prm.No, 0.0);
        o_prm.xi.resize(e_prm.No);
        o_prm.eta.resize(e_prm.Nc);

        o_prm.xi[0].x = 0.3000; o_prm.xi[0].y = 0.5000;
        o_prm.xi[1].x = 0.4000; o_prm.xi[1].y = 0.8000;

        o_prm.eta[0].x = 0.7000; o_prm.eta[0].y = 0.5000;

        o_prm.k[0][0] = -2.1200; o_prm.k[0][1] = -1.1200;
        o_prm.z[0][0] = +1.5000; o_prm.z[0][1] = +0.5000;

        o_prm0 = o_prm;
        */
    }

    {
        /*
        e_prm.Nc = 2;
        e_prm.No = 1;

        o_prm.k.resize(e_prm.Nc, e_prm.No, 0.0);
        o_prm.z.resize(e_prm.Nc, e_prm.No, 0.0);
        o_prm.xi.resize(e_prm.No);
        o_prm.eta.resize(e_prm.Nc);

        o_prm.xi[0].x = 0.3000; o_prm.xi[0].y = 0.5000;

        o_prm.eta[0].x = 0.7000; o_prm.eta[0].y = 0.5000;
        o_prm.eta[1].x = 0.4000; o_prm.eta[1].y = 0.8000;

        o_prm.k[0][0] = -2.1200; o_prm.k[1][0] = -1.1200;
        o_prm.z[0][0] = +1.5000; o_prm.z[1][0] = +1.5000;

        o_prm0 = o_prm;
        */
    }
}

void Problem2HNDirichlet::checkGradient(const Problem2HNDirichlet &prob)
{
    EquationParameter e_prm = prob.mEquParameter;
    OptimizeParameter o_prm = prob.mOptParameter;
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

void Problem2HNDirichlet::optimization1()
{
    //    EquationParameter e_prm;
    //    OptimizeParameter o_prm;
    //    OptimizeParameter o_prm0;
    //    initParameters(e_prm, o_prm, o_prm0);
    //    DoubleVector x;

    //    DoubleVector r; r << 0.1 << 1.0 << 2.0 << 10.0 << 100.0;
    //    for (unsigned int i=0; i<r.length(); i++)
    //    {
    //        calculateWithPenalty();
    //    }


    EquationParameter e_prm;
    OptimizeParameter o_prm;
    OptimizeParameter o_prm0;
    initParameters(e_prm, o_prm, o_prm0);

    double hx, hy; hx = hy = 0.01;
    unsigned Nx, Ny; Nx = Ny = 100;


    Problem2HNDirichlet prob(Dimension(0.01, 0, 200), Dimension(hx, 0, Nx), Dimension(hy, 0, Ny), e_prm, o_prm, o_prm0);
    prob.optimizeK = true;
    prob.optimizeZ = true;
    prob.optimizeC = true;
    prob.optimizeO = true;

    prob.V0.resize(Ny+1, Nx+1, 0.0);

    prob.regEpsilon = 0.0;

    prob.r = 1.00;
    prob.vmin.resize(e_prm.Nc, -5.0);
    prob.vmax.resize(e_prm.Nc, +5.0);

    ConjugateGradient g;
    //SteepestDescentGradient g;
    g.setFunction(&prob);
    g.setGradient(&prob);
    g.setPrinter(&prob);
    g.setProjection(&prob);
    g.setEpsilon1(0.01);
    g.setEpsilon2(0.01);
    g.setEpsilon3(0.01);
    g.setR1MinimizeEpsilon(0.1, 0.001);
    g.setNormalize(true);
    g.showEndMessage(true);
    g.setResetIteration(true);

    DoubleVector x;
    prob.PrmToVector(o_prm, x);
    g.calculate(x);
}

void Problem2HNDirichlet::example1()
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
    DoubleVector e; e << 0.0000 << 0.0000 << 0.0000 << 0.0000 << 0.0000 << 0.0000 << 0.0000;

    double hx, hy; hx = hy = 0.01;
    unsigned Nx, Ny; Nx = Ny = 100;
    double ht = 0.01;
    unsigned int Nt = 500;

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
        prob.LD = 10;

        if (i==0)
        {
            prob.PrmToVector(o_prm, x);


            //DoubleVector cx = x;
            //for (unsigned int i=0; i<=100; i++)
            //{
            //    cx[10] = i*hx;
            //    printf("%f %f\n", cx[10], prob.fx(cx));
            //}

            //checkGradient(prob);
            //IPrinter::printSeperatorLine();
        }

        ConjugateGradient g;
        //SteepestDescentGradient g;
        //g.R1Minimizer().setCallback(new Problem2HNDirichletR1MinimizeCallback);
        g.setFunction(&prob);
        g.setGradient(&prob);
        g.setPrinter(&prob);
        g.setProjection(&prob);
        g.setEpsilon1(0.00001);
        g.setEpsilon2(0.00001);
        g.setEpsilon3(0.00001);
        g.setR1MinimizeEpsilon(1.0, 0.001);
        g.setNormalize(true);
        g.showEndMessage(true);
        g.setResetIteration(true);

        g.calculate(x);

        IPrinter::printSeperatorLine(NULL, '=');
    }
}

void Problem2HNDirichlet::example2()
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
        g.showEndMessage(true);
        g.setResetIteration(true);

        g.calculate(x);

        IPrinter::printSeperatorLine(NULL, '=');
    }
}

void Problem2HNDirichlet::example3()
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

Problem2HNDirichlet::Problem2HNDirichlet()
{
    r = 0.0;
    regEpsilon = 0.0;
}

Problem2HNDirichlet::Problem2HNDirichlet(const Dimension &time, const Dimension &dimx, const Dimension &dimy, const EquationParameter &mEquParameter, const OptimizeParameter &mOptParameter, const OptimizeParameter &mOptParameter0)
{
    setTimeDimension(time);
    addSpaceDimension(dimx);
    addSpaceDimension(dimy);

    this->mEquParameter = mEquParameter;
    this->mOptParameter = mOptParameter;
    this->mRegParameter = mOptParameter0;

    V0.resize(dimy.sizeN()+1, dimx.sizeN()+1, 0.0);

    r = 0.0;
    regEpsilon = 0.0;
}

Problem2HNDirichlet::~Problem2HNDirichlet()
{}

double Problem2HNDirichlet::fx(const DoubleVector &pv) const
{
    OptimizeParameter o_prm;

    VectorToPrm(pv, o_prm);

    Problem2HNDirichlet* prob = const_cast<Problem2HNDirichlet*>(this);
    prob->mOptParameter = o_prm;

    std::vector<DoubleMatrix> u;
    spif_vector u_info;
    prob->solveForwardIBVP(u, u_info, true);

    double intgrl = integral(u);

    for (unsigned int i=0; i<=LD; i++) u[i].clear();
    u.clear();

    double nrm = norm(mEquParameter, o_prm, mRegParameter);
    double pnt = penalty(u_info, o_prm);

    double sum = intgrl + regEpsilon*nrm + r*pnt;

    for (unsigned int j=0; j<u_info.size(); j++)
    {
        u_info[j].clearWeights();
    }
    u_info.clear();

    return sum;
}

double Problem2HNDirichlet::mu(double x UNUSED_PARAM, double y UNUSED_PARAM) const
{
    return 1.0;//sin(M_PI*x) * sin(M_PI*y);
}

double Problem2HNDirichlet::integral(const std::vector<DoubleMatrix> &vu) const
{
    const double ht = timeDimension().step();

    const unsigned int L = timeDimension().sizeN();
    const unsigned int LLD = L + LD;

    double sum = 0.0;

    for (unsigned int l=L; l<=LLD; l++)
    {
        const DoubleMatrix &u = vu.at(2*(l-L));
        if (l == L || l == LLD) sum += 0.5*integralU(u); else sum += integralU(u);
    }

    return sum*ht;
}

double Problem2HNDirichlet::integralU(const DoubleMatrix &u) const
{
    const double hx = spaceDimension(Dimension::DimensionX).step();
    const double hy = spaceDimension(Dimension::DimensionY).step();
    const unsigned int N = spaceDimension(Dimension::DimensionX).sizeN();
    const unsigned int M = spaceDimension(Dimension::DimensionY).sizeN();

    double udiff = 0.0;
    double usum = 0.0;

    udiff = u[0][0] - V0[0][0]; usum += 0.25 * udiff * udiff;// * mu(0.0, 0.0);
    udiff = u[0][N] - V0[0][N]; usum += 0.25 * udiff * udiff;// * mu(1.0, 0.0);
    udiff = u[M][0] - V0[M][0]; usum += 0.25 * udiff * udiff;// * mu(0.0, 1.0);
    udiff = u[M][N] - V0[M][N]; usum += 0.25 * udiff * udiff;// * mu(1.0, 1.0);

    for (unsigned int n=1; n<=N-1; n++)
    {
        udiff = u[0][n] - V0[0][n]; usum += 0.5 * udiff * udiff;// * mu(n*hx, 0.0);
        udiff = u[M][n] - V0[M][n]; usum += 0.5 * udiff * udiff;// * mu(n*hx, 1.0);
    }

    for (unsigned int m=1; m<=M-1; m++)
    {
        udiff = u[m][0] - V0[m][0]; usum += 0.5 * udiff * udiff;// * mu(0.0, m*hy);
        udiff = u[m][N] - V0[m][N]; usum += 0.5 * udiff * udiff;// * mu(1.0, m*hy);
    }

    for (unsigned int m=1; m<=M-1; m++)
    {
        for (unsigned int n=1; n<=N-1; n++)
        {
            udiff = u[m][n] - V0[m][n];
            usum += udiff * udiff;// * mu(n*hx, n*hy);
        }
    }

    return usum*(hx*hy);
}

double Problem2HNDirichlet::norm(const EquationParameter& e_prm, const OptimizeParameter &o_prm, const OptimizeParameter &o_prm0) const
{
    double norm = 0.0;

    for (unsigned int i=0; i<e_prm.Nc; i++)
    {
        norm += (o_prm.eta[i].x - o_prm0.eta[i].x)*(o_prm.eta[i].x - o_prm0.eta[i].x);
        norm += (o_prm.eta[i].y - o_prm0.eta[i].y)*(o_prm.eta[i].y - o_prm0.eta[i].y);

        for (unsigned int j=0; j<e_prm.No; j++)
        {
            norm += (o_prm.k[i][j] - o_prm0.k[i][j])*(o_prm.k[i][j] - o_prm0.k[i][j]);
            norm += (o_prm.z[i][j] - o_prm0.z[i][j])*(o_prm.z[i][j] - o_prm0.z[i][j]);

            norm += (o_prm.xi[j].x - o_prm0.xi[j].x)*(o_prm.xi[j].x - o_prm0.xi[j].x);
            norm += (o_prm.xi[j].y - o_prm0.xi[j].y)*(o_prm.xi[j].y - o_prm0.xi[j].y);
        }
    }

    return norm;
}

double Problem2HNDirichlet::penalty(const spif_vector &info, const OptimizeParameter &o_prm) const
{
    double ht = mtimeDimension.step();
    unsigned int L = mtimeDimension.sizeN();

    double p_sum = 0.0;

    for (unsigned int i=0; i<mEquParameter.Nc; i++)
    {
        double _gp0 = gpi(i, 0, info, o_prm); p_sum += 0.5*_gp0*_gp0;
        for (unsigned int l=1; l<=L-1; l++)
        {
            double _gpi = gpi(i, l, info, o_prm); p_sum += _gpi*_gpi;
        }
        double _gpL = gpi(i, L, info, o_prm); p_sum += 0.5*_gpL*_gpL;
    }

    return p_sum*ht;
}

double Problem2HNDirichlet::gpi(unsigned int i, unsigned int layer, const spif_vector &info, const OptimizeParameter &o_prm) const
{
    double p = fabs(g0i(i, layer, info, o_prm)) - (vmax.at(i) - vmin.at(i))/2.0;
    return p > 0.0 ? p : 0.0;
}

double Problem2HNDirichlet::g0i(unsigned int i, unsigned int layer, const spif_vector &u_info, const OptimizeParameter &o_prm) const
{
    double vi = 0.0;
    for (unsigned int j=0; j<mEquParameter.No; j++)
    {
        const SpacePointInfo &node = u_info[j];
        vi += o_prm.k[i][j] * (node.value(layer)-o_prm.z[i][j]);
    }
    return (vmax.at(i) + vmin.at(i))/2.0 - vi;
}

void Problem2HNDirichlet::gradient(const DoubleVector & pv, DoubleVector &g) const
{
    const unsigned int L = mtimeDimension.sizeN();
    const double ht = mtimeDimension.step();
    const unsigned int Nc = mEquParameter.Nc;
    const unsigned int No = mEquParameter.No;
    const unsigned int LLD = L + LD;

    OptimizeParameter o_prm;
    VectorToPrm(pv, o_prm);

    Problem2HNDirichlet* prob = const_cast<Problem2HNDirichlet*>(this);
    prob->mOptParameter = o_prm;

    std::vector<DoubleMatrix> u;

    spif_vector u_info;
    solveForwardIBVP(u, u_info, true);
    spif_vector p_info;
    solveBackwardIBVP(u, p_info, true, u_info);

    g.clear();
    g.resize(pv.length(), 0.0);
    unsigned int gi = 0;

    // k
    if (optimizeK)
    {
        for (unsigned int i=0; i<Nc; i++)
        {
            const SpacePointInfo &pi = p_info[i];

            for (unsigned int j=0; j<No; j++)
            {
                const SpacePointInfo &uj = u_info[j];

                double grad_Kij = 0.0;
                double zij = o_prm.z[i][j];

                grad_Kij += 0.5 * (pi.value(0) + 2.0*r*gpi(i,0,u_info,o_prm)*sgn(g0i(i,0,u_info,o_prm))) * (uj.value(0) - zij);
                for (unsigned int m=1; m<=LLD-1; m++)
                {
                    grad_Kij += (pi.value(m) + 2.0*r*gpi(i,m,u_info,o_prm)*sgn(g0i(i,m,u_info,o_prm))) * (uj.value(m) - zij);
                }
                grad_Kij += 0.5 * (pi.value(LLD) + 2.0*r*gpi(i,LLD,u_info,o_prm)*sgn(g0i(i,LLD,u_info,o_prm))) * (uj.value(LLD) - zij);
                grad_Kij *= -ht;

                g[gi++] = grad_Kij + 2.0*regEpsilon*(o_prm.k[i][j] - mRegParameter.k[i][j]);
            }
        }
    }
    else
    {
        for (unsigned int i=0; i<Nc; i++)
        {
            for (unsigned int j=0; j<No; j++)
            {
                g[gi++] = 0.0;
            }
        }
    }

    // z
    if (optimizeZ)
    {
        for (unsigned int i=0; i<Nc; i++)
        {
            const SpacePointInfo &pi = p_info[i];

            for (unsigned int j=0; j<No; j++)
            {
                double grad_Zij = 0.0;
                double kij = o_prm.k[i][j];

                grad_Zij += 0.5 * (pi.value(0) + 2.0*r*gpi(i,0,u_info,o_prm)*sgn(g0i(i,0,u_info,o_prm))) * kij;
                for (unsigned int m=1; m<=LLD-1; m++)
                {
                    grad_Zij += (pi.value(m) + 2.0*r*gpi(i,m,u_info,o_prm)*sgn(g0i(i,m,u_info,o_prm))) * kij;
                }
                grad_Zij += 0.5 * (pi.value(LLD) + 2.0*r*gpi(i,LLD,u_info,o_prm)*sgn(g0i(i,LLD,u_info,o_prm))) * kij;
                grad_Zij *= ht;

                g[gi++] = grad_Zij + 2.0*regEpsilon*(o_prm.z[i][j] - mRegParameter.z[i][j]);
            }
        }
    }
    else
    {
        for (unsigned int i=0; i<Nc; i++)
        {
            for (unsigned int j=0; j<No; j++)
            {
                g[gi++] = 0.0;
            }
        }
    }

    // xi
    if (optimizeO)
    {
        for (unsigned int j=0; j<No; j++)
        {
            const SpacePointInfo &uj = u_info[j];

            double gradXijX = 0.0;
            double gradXijY = 0.0;
            double vi = 0.0;

            vi = 0.0;
            for (unsigned int i=0; i<Nc; i++) vi += o_prm.k[i][j] * (p_info[i].value(0) + 2.0*r*gpi(i,0,u_info,o_prm)*sgn(g0i(i,0,u_info,o_prm)));
            gradXijX += 0.5 * uj.valueDx(0) * vi;
            gradXijY += 0.5 * uj.valueDy(0) * vi;

            for (unsigned int m=1; m<=LLD-1; m++)
            {
                vi = 0.0;
                for (unsigned int i=0; i<Nc; i++) vi += o_prm.k[i][j]*(p_info[i].value(m) + 2.0*r*gpi(i,m,u_info,o_prm)*sgn(g0i(i,m,u_info,o_prm)));
                gradXijX += uj.valueDx(m) * vi;
                gradXijY += uj.valueDy(m) * vi;
            }

            vi = 0.0;
            for (unsigned int i=0; i<Nc; i++) vi += o_prm.k[i][j]*(p_info[i].value(LLD) + 2.0*r*gpi(i,LLD,u_info,o_prm)*sgn(g0i(i,LLD,u_info,o_prm)));
            gradXijX += 0.5 * uj.valueDx(LLD) * vi;
            gradXijY += 0.5 * uj.valueDy(LLD) * vi;

            gradXijX *= -ht;
            gradXijY *= -ht;

            g[gi++] = gradXijX + 2.0*regEpsilon*(o_prm.xi[j].x - mRegParameter.xi[j].x);
            g[gi++] = gradXijY + 2.0*regEpsilon*(o_prm.xi[j].y - mRegParameter.xi[j].y);
        }
    }
    else
    {
        for (unsigned int j=0; j<No; j++)
        {
            g[gi++] = 0.0;
            g[gi++] = 0.0;
        }
    }

    // eta
    if (optimizeC)
    {
        for (unsigned int i=0; i<Nc; i++)
        {
            const SpacePointInfo &pi = p_info[i];

            double gradEtaiX = 0.0;
            double gradEtaiY = 0.0;
            double vi = 0.0;

            vi = 0.0;
            for (unsigned int j=0; j<No; j++) vi += o_prm.k[i][j] * (u_info[j].value(0) - o_prm.z[i][j]);
            gradEtaiX += 0.5 * pi.valueDx(0) * vi;
            gradEtaiY += 0.5 * pi.valueDy(0) * vi;

            for (unsigned int m=1; m<=LLD-1; m++)
            {
                vi = 0.0;
                for (unsigned int j=0; j<No; j++) vi += o_prm.k[i][j] * (u_info[j].value(m) - o_prm.z[i][j]);
                gradEtaiX += pi.valueDx(m) * vi;
                gradEtaiY += pi.valueDy(m) * vi;
            }

            vi = 0.0;
            for (unsigned int j=0; j<No; j++) vi += o_prm.k[i][j] * (u_info[j].value(LLD) - o_prm.z[i][j]);
            gradEtaiX += 0.5 * pi.valueDx(LLD) * vi;
            gradEtaiY += 0.5 * pi.valueDy(LLD) * vi;

            gradEtaiX *= -ht;
            gradEtaiY *= -ht;

            g[gi++] = gradEtaiX + 2.0*regEpsilon*(o_prm.eta[i].x - mRegParameter.eta[i].x);
            g[gi++] = gradEtaiY + 2.0*regEpsilon*(o_prm.eta[i].y - mRegParameter.eta[i].y);
        }
    }
    else
    {
        for (unsigned int i=0; i<Nc; i++)
        {
            g[gi++] = 0.0;
            g[gi++] = 0.0;
        }
    }

    for (unsigned int i=0; i<u_info.size(); i++)
    {
        u_info[i].clearWeights();
    }

    for (unsigned int i=0; i<p_info.size(); i++)
    {
        p_info[i].clearWeights();
    }

    u_info.clear();
    p_info.clear();
}

void Problem2HNDirichlet::print(unsigned int i UNUSED_PARAM, const DoubleVector &x, const DoubleVector &g, double f UNUSED_PARAM, GradientMethod::MethodResult result) const
{
    Problem2HNDirichlet* prob = const_cast<Problem2HNDirichlet*>(this);
    OptimizeParameter o_prm;
    VectorToPrm(x, o_prm);

    prob->mOptParameter = o_prm;

    std::vector<DoubleMatrix> u;

    spif_vector u_info;
    solveForwardIBVP(u, u_info, true);

    const char* msg = NULL;
    if (result == GradientMethod::BREAK_FIRST_ITERATION)    msg = "BREAK_FIRST_ITERATION";
    if (result == GradientMethod::FIRST_ITERATION)          msg = "FIRST_ITERATION";
    if (result == GradientMethod::BREAK_GRADIENT_NORM_LESS) msg = "BREAK_GRADIENT_NORM_LESS";
    if (result == GradientMethod::BREAK_DISTANCE_LESS)      msg = "BREAK_DISTANCE_LESS";
    if (result == GradientMethod::NEXT_ITERATION)           msg = "NEXT_ITERATION";

    double ing = integral(u);
    double pnt = penalty(u_info, o_prm);
    double nrm = norm(prob->mEquParameter, prob->mOptParameter, prob->mRegParameter);

    printf("I[%3d]: I1:%8.6f P:%8.6f N:%8.6f F:%8.6f R:%.2f e:%.3f %s\n", i, ing, pnt, nrm, f, r, regEpsilon, msg);
    printf("k:%8.4f %8.4f %8.4f %8.4f z:%8.4f %8.4f %8.4f %8.4f o:%8.4f %8.4f %8.4f %8.4f c:%8.4f %8.4f %8.4f %8.4f\n", x[0], x[1], x[2], x[3], x[4], x[5], x[6], x[7], x[8], x[9], x[10], x[11], x[12], x[13], x[14], x[15]);
    printf("k:%8.4f %8.4f %8.4f %8.4f z:%8.4f %8.4f %8.4f %8.4f o:%8.4f %8.4f %8.4f %8.4f c:%8.4f %8.4f %8.4f %8.4f\n", g[0], g[1], g[2], g[3], g[4], g[5], g[6], g[7], g[8], g[9], g[10], g[11], g[12], g[13], g[14], g[15]);
    DoubleVector n = g;
    n.L2Normalize();
    printf("k:%8.4f %8.4f %8.4f %8.4f z:%8.4f %8.4f %8.4f %8.4f o:%8.4f %8.4f %8.4f %8.4f c:%8.4f %8.4f %8.4f %8.4f\n", n[0], n[1], n[2], n[3], n[4], n[5], n[6], n[7], n[8], n[9], n[10], n[11], n[12], n[13], n[14], n[15]);

    u.clear();
    u_info.clear();

    //    prob->optimizeK = i%2==0;
    //    prob->optimizeZ = i%2==0;
    //    prob->optimizeO = i%2==1;
    //    prob->optimizeC = i%2==1;

    //    if (i==0)
    //    {
    //    prob->optimizeK = i%2==0;//i%4==3;
    //    prob->optimizeZ = i%2==0;//i%4==0;
    //    prob->optimizeO = i%2==1;//i%4==1;
    //    prob->optimizeC = i%2==1;//i%4==2;
    //    }

    C_UNUSED(prob);
    IPrinter::printSeperatorLine();
}

void Problem2HNDirichlet::project(DoubleVector &pv, unsigned int index)
{
    unsigned int Nc = mEquParameter.Nc;
    unsigned int No = mEquParameter.No;

    unsigned int offset = 2*Nc*No;

    // xi
    if ( offset <= index && index <= offset + 2*No - 1 )
    {
        if (pv[index] < 0.05) pv[index] = 0.05;
        if (pv[index] > 0.95) pv[index] = 0.95;
    }

    // eta
    if ( offset + 2*No <= index && index <= offset + 2*No + 2*Nc - 1 )
    {
        if (pv[index] < 0.05) pv[index] = 0.05;
        if (pv[index] > 0.95) pv[index] = 0.95;
    }
    //return;

    double dx = 0.09;

    if (index == 12 && fabs(pv[8] - pv[12])<dx)
    {
        pv[12] = pv[8] + dx;
        if (pv[12] > 0.95) pv[12] = pv[8] - dx;
    }

    if (index == 12 && fabs(pv[10] - pv[12])<dx)
    {
        pv[12] = pv[10] + dx;
        if (pv[12] > 0.95) pv[12] = pv[10] - dx;
    }

    if (index == 14 && fabs(pv[8] - pv[14])<dx)
    {
        pv[14] = pv[8] + dx;
        if (pv[14] > 0.95) pv[14] = pv[8] - dx;
    }
    if (index == 14 && fabs(pv[10] - pv[14])<dx)
    {
        pv[14] = pv[10] + dx;
        if (pv[14] > 0.95) pv[14] = pv[10] - dx;
    }

    if (index == 13 && fabs(pv[9] - pv[13])<dx)
    {
        pv[13] = pv[9] + dx;
        if (pv[13] > 0.95) pv[13] = pv[9] - dx;
    }
    if (index == 13 && fabs(pv[11] - pv[13])<dx)
    {
        pv[13] = pv[11] + dx;
        if (pv[13] > 0.95) pv[13] = pv[11] - dx;
    }

    if (index == 15 && fabs(pv[9] - pv[15])<dx)
    {
        pv[15] = pv[9] + dx;
        if (pv[15] > 0.95) pv[15] = pv[9] - dx;
    }
    if (index == 15 && fabs(pv[11] - pv[15])<dx)
    {
        pv[15] = pv[11] + dx;
        if (pv[15] > 0.95) pv[15] = pv[11] - dx;
    }
}

void Problem2HNDirichlet::solveForwardIBVP(std::vector<DoubleMatrix> &u, spif_vector &u_info, bool use) const
{
    solveForwardIBVP1(u, u_info, use);
}

void Problem2HNDirichlet::solveBackwardIBVP(const std::vector<DoubleMatrix> &u, spif_vector &p_info, bool use, const spif_vector &u_info) const
{
    solveBackwardIBVP1(u, p_info, use, u_info);
}
/**
 * @brief Problem2HNDirichlet::solveForwardIBVP
 * @param u
 * @param u_info
 * @param use
 */
void Problem2HNDirichlet::solveForwardIBVP0(std::vector<DoubleMatrix> &u, spif_vector &u_info, bool use) const
{
    const Dimension dimX = spaceDimension(Dimension::DimensionX);
    const Dimension dimY = spaceDimension(Dimension::DimensionY);
    const Dimension time = timeDimension();

    const unsigned int N = dimX.sizeN();
    const unsigned int M = dimY.sizeN();
    const unsigned int L = time.sizeN();
    const unsigned int LLD = L+LD;

    const double hx = dimX.step();
    const double hy = dimY.step();
    const double ht = time.step();

    const double a = mEquParameter.a;
    const double lambda = mEquParameter.lambda;
    const unsigned int No = mEquParameter.No;
    const unsigned int Nc = mEquParameter.Nc;
    const unsigned int Ns = mEquParameter.Ns;

    const double m_aa_htht__hxhx = -(a*a*ht*ht)/(hx*hx);
    const double p_aa_htht__hxhx___lambda_ht = +2.0 + 2.0*(a*a*ht*ht)/(hx*hx) + 3.0*(lambda*ht);
    const double p_aa_htht__hyhy = +(a*a*ht*ht)/(hy*hy);

    const double m_aa_htht__hyhy = -(a*a*ht*ht)/(hy*hy);
    const double p_aa_htht__hyhy___lambda_ht = +2.0 + 2.0*(a*a*ht*ht)/(hy*hy) + 3.0*(lambda*ht);
    const double p_aa_htht__hxhx = +(a*a*ht*ht)/(hx*hx);

    const double htht = ht*ht;
    const double lambda_ht = lambda*ht;

    const double aa__hxhx = (a*a)/(hx*hx);
    const double aa__hyhy = (a*a)/(hy*hy);

    DoubleMatrix u00(M+1, N+1);
    DoubleMatrix u05(M+1, N+1);
    DoubleMatrix u10(M+1, N+1);
    DoubleMatrix u15(M+1, N+1);
    DoubleMatrix u20(M+1, N+1);

    for (unsigned int l=0; l<u.size(); l++) u[l].clear(); u.clear();
    u.resize(2*LD+1); for (unsigned int l=0; l<u.size(); l++) u[l].resize(M+1, N+1);

    //----------------------------------------------------------------------------------------------//
    espn_vector obsPointNodes, cntDeltaNodes, qPointNodes;
    for (unsigned int j=0; j<No; j++) distributeDelta0(mOptParameter.xi[j], j, obsPointNodes, dimX, dimY, 4, 4);
    for (unsigned int i=0; i<Nc; i++) distributeDelta0(mOptParameter.eta[i], i, cntDeltaNodes, dimX, dimY, 4, 4);
    for (unsigned int s=0; s<Ns; s++) distributeDelta0(mEquParameter.theta[s], s, qPointNodes, dimX, dimY, 4, 4);

    //----------------------------------------------------------------------------------------------//
    uint_vector rows0, rows1, rows2, cols0, cols1, cols2;
    f_findRowsCols(rows0, rows1, rows2, cols0, cols1, cols2, obsPointNodes, cntDeltaNodes, N, M);

    //-------------------------------------------- info --------------------------------------------//
    if (use == true) f_prepareInfo(No, mOptParameter.xi, u_info, LLD, dimX, dimY);
    //----------------------------------------------------------------------------------------------//

    //------------------------------------- initial conditions -------------------------------------//
    f_initialLayers(u00, u05, u10, u_info, use, obsPointNodes, cntDeltaNodes, qPointNodes, N, M, hx, hy, ht, aa__hxhx, aa__hyhy, lambda);
    //------------------------------------- initial conditions -------------------------------------//

    double *ax = (double *) malloc(sizeof(double)*(N-1)); for (unsigned int n=1; n<=N-1; n++) ax[n-1] = m_aa_htht__hxhx;
    double *bx = (double *) malloc(sizeof(double)*(N-1)); for (unsigned int n=1; n<=N-1; n++) bx[n-1] = p_aa_htht__hxhx___lambda_ht;
    double *cx = (double *) malloc(sizeof(double)*(N-1)); for (unsigned int n=1; n<=N-1; n++) cx[n-1] = m_aa_htht__hxhx;
    double *dx = (double *) malloc(sizeof(double)*(N-1));
    double *rx = (double *) malloc(sizeof(double)*(N-1));
    ax[0] = cx[N-2] = 0.0;

    double *ay = (double *) malloc(sizeof(double)*(M-1)); for (unsigned int m=1; m<=M-1; m++) ay[m-1] = m_aa_htht__hyhy;
    double *by = (double *) malloc(sizeof(double)*(M-1)); for (unsigned int m=1; m<=M-1; m++) by[m-1] = p_aa_htht__hyhy___lambda_ht;
    double *cy = (double *) malloc(sizeof(double)*(M-1)); for (unsigned int m=1; m<=M-1; m++) cy[m-1] = m_aa_htht__hyhy;
    double *dy = (double *) malloc(sizeof(double)*(M-1));
    double *ry = (double *) malloc(sizeof(double)*(M-1));
    ay[0] = cy[M-2] = 0.0;

    SpaceNodePDE sn;

    for (unsigned int l=2; l<=LLD; l++)
    {
        TimeNodePDE tn20; tn20.i = l; tn20.t = l*ht;
        TimeNodePDE tn15; tn15.i = l; tn15.t = l*ht-0.5*ht;

        /**************************************************** border conditions ***************************************************/

        SpaceNodePDE sn0;
        SpaceNodePDE sn1;

        sn0.i = 0; sn0.x = 0.0;
        sn1.i = N; sn1.x = hx*N;
        for (unsigned int m=0; m<=M; m++)
        {
            sn0.j = m; sn0.y = m*hy; u15[m][0] = f_boundary(sn0, tn15); u20[m][0] = f_boundary(sn0, tn20);
            sn1.j = m; sn1.y = m*hy; u15[m][N] = f_boundary(sn1, tn15); u20[m][N] = f_boundary(sn1, tn20);
        }

        sn0.j = 0; sn0.y = 0.0;
        sn1.j = M; sn1.y = hy*M;
        for (unsigned int n=0; n<=N; n++)
        {
            sn0.i = n; sn0.x = n*hx; u15[0][n] = f_boundary(sn0, tn15); u20[0][n] = f_boundary(sn0, tn20);
            sn1.i = n; sn1.x = n*hx; u15[M][n] = f_boundary(sn1, tn15); u20[M][n] = f_boundary(sn1, tn20);
        }

        /**************************************************** border conditions ***************************************************/

        /**************************************************** x direction apprx ***************************************************/

        if (rows0.size() != 0)
        {
            for (unsigned int row=0; row<rows0.size(); row++)
            {
                unsigned int m = rows0.at(row);
                sn.j = m; sn.y = m*hy;
                for (unsigned int n=1; n<=N-1; n++)
                {
                    sn.i = n; sn.x = n*hx;

                    dx[n-1] = 0.0;
                    if (m>0 && m<M)  dx[n-1] = p_aa_htht__hyhy*(u10[m-1][n] - 2.0*u10[m][n]   + u10[m+1][n]);
                    else if (m == 0) dx[n-1] = p_aa_htht__hyhy*(u10[0][n]   - 2.0*u10[1][n]   + u10[2][n]);
                    else if (m == M) dx[n-1] = p_aa_htht__hyhy*(u10[M-2][n] - 2.0*u10[M-1][n] + u10[M][n]);

                    dx[n-1] += (u10[m][n]-u00[m][n]) + 2.0*u10[m][n];
                    dx[n-1] += lambda_ht*(4.0*u10[m][n]-u05[m][n]);

                    //if (l==2)
                    //{
                    //    for (unsigned int si=0; si<qPointNodes.size(); si++)
                    //    {
                    //        const ExtendedSpacePointNode &qNode = qPointNodes.at(si);
                    //        if (qNode.i == sn.i && qNode.j == sn.j)
                    //        {
                    //            dx[n-1] += (mEquParameter.q[qNode.id] * qNode.w * (2.0/ht))*((ht*ht)*0.5);
                    //        }
                    //    }
                    //}
                }

                dx[0]   -= m_aa_htht__hxhx * u15[m][0];
                dx[N-2] -= m_aa_htht__hxhx * u15[m][N];

                tomasAlgorithm(ax, bx, cx, dx, rx, N-1);
                for (unsigned int n=1; n<=N-1; n++) u15[m][n] = rx[n-1];
            }
        }

        if (rows1.size() != 0 && rows2.size() == 0)
        {
            //throw std::exception();

            double* U15 = (double *) malloc(sizeof(double)*No);
            for (unsigned int j=0; j<No; j++) U15[j] = 0.0;
            for (unsigned int odj=0; odj<obsPointNodes.size(); odj++)
            {
                const ExtendedSpacePointNode &opn = obsPointNodes[odj];
                U15[opn.id] += u15[opn.j][opn.i] * (opn.w * (hx*hy));
            }

            double *v = new double[Nc];
            for (unsigned int i=0; i<Nc; i++)
            {
                v[i] = 0.0;
                for (unsigned int j=0; j<No; j++)
                {
                    v[i] += mOptParameter.k[i][j] * (U15[j]-mOptParameter.z[i][j]);
                }
            }
            //printf("l %d v: %f %f\n", l, v[0], v[1]);

            for (unsigned int row=0; row<rows1.size(); row++)
            {
                unsigned int m = rows1.at(row);
                sn.j = m; sn.y = m*hy;
                for (unsigned int n=1; n<=N-1; n++)
                {
                    sn.i = n; sn.x = n*hx;

                    dx[n-1] = 0.0;
                    if (m>0 && m<M)  dx[n-1] = p_aa_htht__hyhy*(u10[m-1][n] - 2.0*u10[m][n]   + u10[m+1][n]);
                    else if (m == 0) dx[n-1] = p_aa_htht__hyhy*(u10[0][n]   - 2.0*u10[1][n]   + u10[2][n]);
                    else if (m == M) dx[n-1] = p_aa_htht__hyhy*(u10[M-2][n] - 2.0*u10[M-1][n] + u10[M][n]);

                    dx[n-1] += (u10[m][n]-u00[m][n]) + 2.0*u10[m][n];
                    dx[n-1] += lambda_ht*(4.0*u10[m][n]-u05[m][n]);

                    //------------------------------------- Adding delta part -------------------------------------//
                    for (unsigned int cni=0; cni<cntDeltaNodes.size(); cni++)
                    {
                        const ExtendedSpacePointNode &cdn = cntDeltaNodes.at(cni);
                        if (cdn.i == sn.i && cdn.j == sn.j)
                        {
                            dx[n-1] += htht * v[cdn.id] * cdn.w;
                            //for (unsigned int j=0; j<No; j++)
                            //{
                            //    dx[n-1] += htht * mOptParameter.k[cdn.id][j] * (U15[j]-mOptParameter.z[cdn.id][j]) * cdn.w;
                            //}
                        }
                    }
                    //------------------------------------- Adding delta part -------------------------------------//

                    //if (l==2)
                    //{
                    //    for (unsigned int si=0; si<qPointNodes.size(); si++)
                    //    {
                    //        const ExtendedSpacePointNode &qNode = qPointNodes.at(si);
                    //        if (qNode.i == sn.i && qNode.j == sn.j)
                    //        {
                    //            dx[n-1] += (mEquParameter.q[qNode.id] * qNode.w * (2.0/ht)) * htht_h;
                    //        }
                    //    }
                    //}
                }

                dx[0]   -= m_aa_htht__hxhx * u15[m][0];
                dx[N-2] -= m_aa_htht__hxhx * u15[m][N];

                tomasAlgorithm(ax, bx, cx, dx, rx, N-1);
                for (unsigned int n=1; n<=N-1; n++) u15[m][n] = rx[n-1];
            }

            delete [] v;
            free(U15);
        }

        if (rows1.size() != 0 && rows2.size() != 0)
        {
            //throw std::exception();

            unsigned int row1_size = rows1.size()*(N-1);
            double* a1 = (double*) malloc(sizeof(double)*row1_size);
            double* b1 = (double*) malloc(sizeof(double)*row1_size);
            double* c1 = (double*) malloc(sizeof(double)*row1_size);
            double* d1 = (double*) malloc(sizeof(double)*row1_size);
            double* x1 = (double*) malloc(sizeof(double)*row1_size);
            DoubleMatrix w1(row1_size, row1_size, 0.0);

            unsigned int offset = 0;
            for (unsigned int row=0; row<rows1.size(); row++)
            {
                unsigned int m = rows1.at(row);
                sn.j = m; sn.y = m*hy;

                for (unsigned int n=1; n<=N-1; n++)
                {
                    sn.i = n; sn.x = n*hx;

                    d1[offset+(n-1)] = 0.0;
                    if (m>0 && m<M) d1[offset+(n-1)] = p_aa_htht__hyhy*(u10[m-1][n] - 2.0*u10[m][n]   + u10[m+1][n]);
                    else if (m==0)  d1[offset+(n-1)] = p_aa_htht__hyhy*(u10[0][n]   - 2.0*u10[1][n]   + u10[2][n]);
                    else if (m==M)  d1[offset+(n-1)] = p_aa_htht__hyhy*(u10[M-2][n] - 2.0*u10[M-1][n] + u10[M][n]);

                    d1[offset+(n-1)] += (u10[m][n]-u00[m][n]) + 2.0*u10[m][n];
                    d1[offset+(n-1)] += lambda_ht*(4.0*u10[m][n]-u05[m][n]);

                    a1[offset+(n-1)] = m_aa_htht__hxhx;
                    b1[offset+(n-1)] = p_aa_htht__hxhx___lambda_ht;
                    c1[offset+(n-1)] = m_aa_htht__hxhx;

                    //------------------------------------- Adding delta part -------------------------------------//
                    for (unsigned int cni=0; cni<cntDeltaNodes.size(); cni++)
                    {
                        const ExtendedSpacePointNode &cdn = cntDeltaNodes.at(cni);
                        if (cdn.i == sn.i && cdn.j == sn.j)
                        {
                            for (unsigned int odj=0; odj<obsPointNodes.size(); odj++)
                            {
                                const ExtendedSpacePointNode &opn = obsPointNodes.at(odj);

                                bool found = false;
                                for (unsigned int rs=0; rs<rows1.size(); rs++)
                                {
                                    if (opn.j == rows1[rs])
                                    {
                                        found = true;
                                        w1[offset+(n-1)][rs*(N-1)+(opn.i-1)] -= htht * mOptParameter.k[cdn.id][opn.id] * (opn.w * (hx*hy)) * cdn.w;
                                        break;
                                    }
                                }

                                if (!found)
                                {
                                    d1[offset+(n-1)] += htht * mOptParameter.k[cdn.id][opn.id] * u15[opn.j][opn.i] * (opn.w * (hx*hy)) * cdn.w;
                                }
                            }

                            for (unsigned int j=0; j<No; j++)
                            {
                                d1[offset+(n-1)] -= htht * mOptParameter.k[cdn.id][j] * mOptParameter.z[cdn.id][j] * cdn.w;
                            }
                        }
                    }
                    //------------------------------------- Adding delta part -------------------------------------//

                    //if (l==2)
                    //{
                    //    for (unsigned int si=0; si<qPointNodes.size(); si++)
                    //    {
                    //        const ExtendedSpacePointNode &qNode = qPointNodes.at(si);
                    //        if (qNode.i == sn.i && qNode.j == sn.j)
                    //        {
                    //            d1[offset+(n-1)] += (mEquParameter.q[qNode.id] * qNode.w * (2.0/ht)) * htht_h;
                    //        }
                    //    }
                    //}
                }

                a1[offset+0]   = 0.0;
                c1[offset+N-2] = 0.0;

                d1[offset+0]   -= m_aa_htht__hxhx * u15[m][0];
                d1[offset+N-2] -= m_aa_htht__hxhx * u15[m][N];

                offset += N-1;
            }

            LinearEquation::func1(a1, b1, c1, d1, w1.data(), x1, row1_size);

            offset = 0;
            for (unsigned int row=0; row<rows1.size(); row++)
            {
                unsigned int m=rows1.at(row);
                for (unsigned int n=1; n<=N-1; n++)
                {
                    u15[m][n] = x1[offset+(n-1)];
                }
                offset += N-1;
            }

            w1.clear();
            free(x1);
            free(d1);
            free(c1);
            free(b1);
            free(a1);
        }

        /**************************************************** x direction apprx ***************************************************/

        /**************************************************** y direction apprx ***************************************************/

        if (cols0.size() != 0)
        {
            for (unsigned int col=0; col<cols0.size(); col++)
            {
                unsigned int n = cols0.at(col);
                sn.i = n; sn.x = n*hx;
                for (unsigned int m=1; m<=M-1; m++)
                {
                    sn.j = m; sn.y = m*hy;
                    \
                    dy[m-1] = 0.0;
                    if (n>0 && n<N) dy[m-1] = p_aa_htht__hxhx*(u15[m][n-1] - 2.0*u15[m][n]   + u15[m][n+1]);
                    else if (n==0)  dy[m-1] = p_aa_htht__hxhx*(u15[m][0]   - 2.0*u15[m][1]   + u15[m][2]);
                    else if (n==N)  dy[m-1] = p_aa_htht__hxhx*(u15[m][N-2] - 2.0*u15[m][N-1] + u15[m][N]);

                    dy[m-1] += (u10[m][n]-u00[m][n]) + 2.0*u15[m][n];
                    dy[m-1] += lambda_ht*(4.0*u15[m][n]-u10[m][n]);

                    if (l==2)
                    {
                        for (unsigned int si=0; si<qPointNodes.size(); si++)
                        {
                            const ExtendedSpacePointNode &qNode = qPointNodes.at(si);
                            if (qNode.i == sn.i && qNode.j == sn.j)
                            {
                                dy[m-1] += (mEquParameter.q[qNode.id] * qNode.w * (1.0/ht)) * htht;
                            }
                        }
                    }
                }

                dy[0]   -= m_aa_htht__hyhy * u20[0][n];
                dy[M-2] -= m_aa_htht__hyhy * u20[M][n];

                tomasAlgorithm(ay, by, cy, dy, ry, M-1);
                for (unsigned int m=1; m<=M-1; m++) u20[m][n] = ry[m-1];
            }
        }

        if (cols1.size() != 0 && cols2.size() == 0)
        {
            //throw std::exception();

            double* U20 = (double *) malloc(sizeof(double)*No);
            for (unsigned int j=0; j<No; j++) U20[j] = 0.0;
            for (unsigned int odj=0; odj<obsPointNodes.size(); odj++)
            {
                const ExtendedSpacePointNode &opn = obsPointNodes[odj];
                U20[opn.id] += u20[opn.j][opn.i] * (opn.w * (hx*hy));
            }
            //printf("%d %d U: %f %f %d\n", obsPointNodes.size(), cntDeltaNodes.size(), U20[0], U20[1], cols1.size());

            double *v = new double[Nc];
            for (unsigned int i=0; i<Nc; i++)
            {
                v[i] = 0.0;
                for (unsigned int j=0; j<No; j++)
                {
                    v[i] += mOptParameter.k[i][j] * (U20[j]-mOptParameter.z[i][j]);
                }
            }

            for (unsigned int col=0; col<cols1.size(); col++)
            {
                unsigned int n = cols1.at(col);
                sn.i = n; sn.x = n*hx;
                for (unsigned int m=1; m<=M-1; m++)
                {
                    sn.j = m; sn.y = m*hy;

                    dy[m-1] = 0.0;
                    if (n>0 && n<N) dy[m-1] = p_aa_htht__hxhx*(u15[m][n-1] - 2.0*u15[m][n]   + u15[m][n+1]);
                    else if (n==0)  dy[m-1] = p_aa_htht__hxhx*(u15[m][0]   - 2.0*u15[m][1]   + u15[m][2]);
                    else if (n==N)  dy[m-1] = p_aa_htht__hxhx*(u15[m][N-2] - 2.0*u15[m][N-1] + u15[m][N]);

                    dy[m-1] += (u10[m][n]-u00[m][n]) + 2.0*u15[m][n];
                    dy[m-1] += lambda_ht*(4.0*u15[m][n]-u10[m][n]);

                    //------------------------------------- Adding delta part -------------------------------------//
                    for (unsigned int cni=0; cni<cntDeltaNodes.size(); cni++)
                    {
                        const ExtendedSpacePointNode &cdn = cntDeltaNodes.at(cni);
                        if (cdn.i == sn.i && cdn.j == sn.j)
                        {
                            dy[m-1] += htht * v[cdn.id] * cdn.w;
                            //for (unsigned int j=0; j<No; j++)
                            //{
                            //    if (cols1[col]==48 && l==3) printf("%d %d %f %d %d\n", cdn.id, j, cdn.w, cdn.i, cdn.j);
                            //    dy[m-1] += htht * mOptParameter.k[cdn.id][j] * (U20[j]-mOptParameter.z[cdn.id][j]) * cdn.w;
                            //}
                        }
                    }
                    //------------------------------------- Adding delta part -------------------------------------//

                    if (l==2)
                    {
                        for (unsigned int si=0; si<qPointNodes.size(); si++)
                        {
                            const ExtendedSpacePointNode &qNode = qPointNodes.at(si);
                            if (qNode.i == sn.i && qNode.j == sn.j)
                            {
                                dy[m-1] += (mEquParameter.q[qNode.id] * qNode.w * (1.0/ht)) * htht;
                            }
                        }
                    }
                }

                dy[0]   -= m_aa_htht__hyhy * u20[0][n];
                dy[M-2] -= m_aa_htht__hyhy * u20[M][n];

                tomasAlgorithm(ay, by, cy, dy, ry, M-1);
                for (unsigned int m=1; m<=M-1; m++) u20[m][n] = ry[m-1];
            }

            delete [] v;
            free(U20);
        }

        if (cols1.size() != 0 && cols2.size() != 0)
        {
            //throw std::exception();

            unsigned int cols1_size = cols1.size()*(M-1);
            double* a2 = (double*) malloc(sizeof(double)*cols1_size);
            double* b2 = (double*) malloc(sizeof(double)*cols1_size);
            double* c2 = (double*) malloc(sizeof(double)*cols1_size);
            double* d2 = (double*) malloc(sizeof(double)*cols1_size);
            double* x2 = (double*) malloc(sizeof(double)*cols1_size);
            DoubleMatrix w2(cols1_size, cols1_size, 0.0);

            unsigned int offset = 0;
            for (unsigned int col=0; col<cols1.size(); col++)
            {
                unsigned int n = cols1.at(col);
                sn.i = n; sn.x = n*hx;

                for (unsigned int m=1; m<=M-1; m++)
                {
                    sn.j = m; sn.y = m*hy;

                    d2[offset+(m-1)] = 0.0;
                    if (n>0 && n<N) d2[offset+(m-1)] = p_aa_htht__hxhx*(u15[m][n-1] - 2.0*u15[m][n]   + u15[m][n+1]);
                    else if (n==0)  d2[offset+(m-1)] = p_aa_htht__hxhx*(u15[m][0]   - 2.0*u15[m][1]   + u15[m][2]);
                    else if (n==N)  d2[offset+(m-1)] = p_aa_htht__hxhx*(u15[m][N-2] - 2.0*u15[m][N-1] + u15[m][N]);

                    d2[offset+(m-1)] += (u10[m][n] - u00[m][n]) + 2.0*u15[m][n];
                    d2[offset+(m-1)] += lambda_ht*(4.0*u15[m][n] - u10[m][n]);

                    a2[offset+(m-1)] = m_aa_htht__hyhy;
                    b2[offset+(m-1)] = p_aa_htht__hyhy___lambda_ht;
                    c2[offset+(m-1)] = m_aa_htht__hyhy;

                    //------------------------------------- Adding delta part -------------------------------------//
                    for (unsigned int cni=0; cni<cntDeltaNodes.size(); cni++)
                    {
                        const ExtendedSpacePointNode &cdn = cntDeltaNodes.at(cni);
                        if (cdn.i == sn.i && cdn.j == sn.j)
                        {
                            for (unsigned int onj=0; onj<obsPointNodes.size(); onj++)
                            {
                                const ExtendedSpacePointNode &opn = obsPointNodes.at(onj);

                                bool found = false;
                                for (unsigned int cs=0; cs<cols1.size(); cs++)
                                {
                                    if (opn.i == cols1[cs])
                                    {
                                        found = true;
                                        w2[offset+(m-1)][cs*(M-1)+(opn.j-0)] -= htht * mOptParameter.k[cdn.id][opn.id] * (opn.w * (hx*hy)) * cdn.w;
                                        break;
                                    }
                                }

                                if (!found)
                                {
                                    d2[offset+(m-1)] += htht * mOptParameter.k[cdn.id][opn.id] * u20[opn.j][opn.i] * (opn.w * (hx*hy)) * cdn.w;
                                }
                            }
                            for (unsigned int j=0; j<No; j++)
                            {
                                d2[offset+(m-1)] -= htht * mOptParameter.k[cdn.id][j] * mOptParameter.z[cdn.id][j] * cdn.w;
                            }
                        }
                    }
                    //------------------------------------- Adding delta part -------------------------------------//

                    if (l==2)
                    {
                        for (unsigned int si=0; si<qPointNodes.size(); si++)
                        {
                            const ExtendedSpacePointNode &qNode = qPointNodes.at(si);
                            if (qNode.i == sn.i && qNode.j == sn.j)
                            {
                                d2[offset+(m-1)] += (mEquParameter.q[qNode.id] * qNode.w * (1.0/ht)) * htht;
                            }
                        }
                    }
                }

                a2[offset+0]   = 0.0;
                c2[offset+M-2] = 0.0;

                d2[offset+0]   -= m_aa_htht__hyhy * u20[0][n];
                d2[offset+M-2] -= m_aa_htht__hyhy * u20[M][n];

                offset += M-1;
            }

            LinearEquation::func1(a2, b2, c2, d2, w2.data(), x2, cols1.size()*(M-1));

            offset = 0;
            for (unsigned int col=0; col<cols1.size(); col++)
            {
                unsigned int n=cols1.at(col);
                for (unsigned int m=1; m<=M-1; m++)
                {
                    u20[m][n] = x2[offset+(m-1)];
                }
                offset += M-1;
            }

            w2.clear();
            free(x2);
            free(d2);
            free(c2);
            free(b2);
            free(a2);
        }

        /**************************************************** y direction apprx ***************************************************/

        for (unsigned int m=0; m<=M; m++)
        {
            for (unsigned int n=0; n<=N; n++)
            {
                u00[m][n] = u10[m][n];
                u10[m][n] = u20[m][n];
                u05[m][n] = u15[m][n];
            }
        }

        if (use == true) f_add2Info(u20, u_info, obsPointNodes, l, hx, hy);

        f_layerInfo(u20, l);

        if (L == l)
        {
            for (unsigned int m=0; m<=M; m++)
            {
                for (unsigned int n=0; n<=N; n++)
                {
                    u[l-L][m][n] = u20[m][n];
                }
            }
        }

        if ( L+1 <= l && l <= LLD )
        {
            for (unsigned int m=0; m<=M; m++)
            {
                for (unsigned int n=0; n<=N; n++)
                {
                    u[2*(l-L)-1][m][n] = u15[m][n];
                    u[2*(l-L)+0][m][n] = u20[m][n];
                }
            }
        }
    }

    free(rx);
    free(dx);
    free(cx);
    free(bx);
    free(ax);

    free(ry);
    free(dy);
    free(cy);
    free(by);
    free(ay);

    rows0.clear();
    rows1.clear();
    rows2.clear();

    cols0.clear();
    cols1.clear();
    cols2.clear();

    qPointNodes.clear();
    obsPointNodes.clear();
    cntDeltaNodes.clear();

    u00.clear();
    u05.clear();
    u10.clear();
    u15.clear();
    u20.clear();
}

void Problem2HNDirichlet::solveForwardIBVP1(std::vector<DoubleMatrix> &u, spif_vector &u_info, bool use) const
{
    const Dimension dimX = spaceDimension(Dimension::DimensionX);
    const Dimension dimY = spaceDimension(Dimension::DimensionY);
    const Dimension time = timeDimension();

    const unsigned int N = dimX.sizeN();
    const unsigned int M = dimY.sizeN();
    const unsigned int L = time.sizeN();
    const unsigned int LLD = L+LD;

    const double hx = dimX.step();
    const double hy = dimY.step();
    const double ht = time.step();

    const double a = mEquParameter.a;
    const double lambda = mEquParameter.lambda;
    const unsigned int No = mEquParameter.No;
    const unsigned int Nc = mEquParameter.Nc;
    const unsigned int Ns = mEquParameter.Ns;

    const double m_aa_htht__hxhx = -(a*a*ht*ht)/(hx*hx);
    const double p_aa_htht__hxhx___lambda_ht = +2.0 + 2.0*(a*a*ht*ht)/(hx*hx) + (lambda*ht);
    const double p_aa_htht__hyhy = +(a*a*ht*ht)/(hy*hy);

    const double m_aa_htht__hyhy = -(a*a*ht*ht)/(hy*hy);
    const double p_aa_htht__hyhy___lambda_ht = +2.0 + 2.0*(a*a*ht*ht)/(hy*hy) + (lambda*ht);
    const double p_aa_htht__hxhx = +(a*a*ht*ht)/(hx*hx);

    const double htht = ht*ht;
    const double lambda_ht = lambda*ht;

    const double aa__hxhx = (a*a)/(hx*hx);
    const double aa__hyhy = (a*a)/(hy*hy);

    DoubleMatrix u00(M+1, N+1);
    DoubleMatrix u10(M+1, N+1);
    DoubleMatrix u15(M+1, N+1);
    DoubleMatrix u20(M+1, N+1);

    for (unsigned int l=0; l<u.size(); l++) u[l].clear(); u.clear();
    u.resize(2*LD+1); for (unsigned int l=0; l<u.size(); l++) u[l].resize(M+1, N+1);

    //----------------------------------------------------------------------------------------------//
    espn_vector obsPointNodes, cntDeltaNodes, qPointNodes;
    for (unsigned int j=0; j<No; j++) distributeDelta0(mOptParameter.xi[j], j, obsPointNodes, dimX, dimY, 4, 4);
    for (unsigned int i=0; i<Nc; i++) distributeDelta0(mOptParameter.eta[i], i, cntDeltaNodes, dimX, dimY, 4, 4);
    for (unsigned int s=0; s<Ns; s++) distributeDelta0(mEquParameter.theta[s], s, qPointNodes, dimX, dimY, 4, 4);

    //----------------------------------------------------------------------------------------------//
    uint_vector rows0, rows1, rows2, cols0, cols1, cols2;
    f_findRowsCols(rows0, rows1, rows2, cols0, cols1, cols2, obsPointNodes, cntDeltaNodes, N, M);

    //-------------------------------------------- info --------------------------------------------//
    if (use == true) f_prepareInfo(No, mOptParameter.xi, u_info, LLD, dimX, dimY);
    //----------------------------------------------------------------------------------------------//

    //------------------------------------- initial conditions -------------------------------------//
    f_initialLayers1(u00, u10, u_info, use, obsPointNodes, cntDeltaNodes, qPointNodes, N, M, hx, hy, ht, aa__hxhx, aa__hyhy, lambda);
    //------------------------------------- initial conditions -------------------------------------//

    double *ax = (double *) malloc(sizeof(double)*(N-1)); for (unsigned int n=1; n<=N-1; n++) ax[n-1] = m_aa_htht__hxhx;
    double *bx = (double *) malloc(sizeof(double)*(N-1)); for (unsigned int n=1; n<=N-1; n++) bx[n-1] = p_aa_htht__hxhx___lambda_ht;
    double *cx = (double *) malloc(sizeof(double)*(N-1)); for (unsigned int n=1; n<=N-1; n++) cx[n-1] = m_aa_htht__hxhx;
    double *dx = (double *) malloc(sizeof(double)*(N-1));
    double *rx = (double *) malloc(sizeof(double)*(N-1));
    ax[0] = cx[N-2] = 0.0;

    double *ay = (double *) malloc(sizeof(double)*(M-1)); for (unsigned int m=1; m<=M-1; m++) ay[m-1] = m_aa_htht__hyhy;
    double *by = (double *) malloc(sizeof(double)*(M-1)); for (unsigned int m=1; m<=M-1; m++) by[m-1] = p_aa_htht__hyhy___lambda_ht;
    double *cy = (double *) malloc(sizeof(double)*(M-1)); for (unsigned int m=1; m<=M-1; m++) cy[m-1] = m_aa_htht__hyhy;
    double *dy = (double *) malloc(sizeof(double)*(M-1));
    double *ry = (double *) malloc(sizeof(double)*(M-1));
    ay[0] = cy[M-2] = 0.0;

    SpaceNodePDE sn;

    for (unsigned int l=2; l<=LLD; l++)
    {
        TimeNodePDE tn20; tn20.i = l; tn20.t = l*ht;
        TimeNodePDE tn15; tn15.i = l; tn15.t = l*ht-0.5*ht;

        /**************************************************** border conditions ***************************************************/

        SpaceNodePDE sn0;
        SpaceNodePDE sn1;

        sn0.i = 0; sn0.x = 0.0;
        sn1.i = N; sn1.x = hx*N;
        for (unsigned int m=0; m<=M; m++)
        {
            sn0.j = m; sn0.y = m*hy; u15[m][0] = f_boundary(sn0, tn15); u20[m][0] = f_boundary(sn0, tn20);
            sn1.j = m; sn1.y = m*hy; u15[m][N] = f_boundary(sn1, tn15); u20[m][N] = f_boundary(sn1, tn20);
        }

        sn0.j = 0; sn0.y = 0.0;
        sn1.j = M; sn1.y = hy*M;
        for (unsigned int n=0; n<=N; n++)
        {
            sn0.i = n; sn0.x = n*hx; u15[0][n] = f_boundary(sn0, tn15); u20[0][n] = f_boundary(sn0, tn20);
            sn1.i = n; sn1.x = n*hx; u15[M][n] = f_boundary(sn1, tn15); u20[M][n] = f_boundary(sn1, tn20);
        }

        /**************************************************** border conditions ***************************************************/

        /**************************************************** x direction apprx ***************************************************/

        if (rows0.size() != 0)
        {
            for (unsigned int row=0; row<rows0.size(); row++)
            {
                unsigned int m = rows0.at(row);
                sn.j = m; sn.y = m*hy;
                for (unsigned int n=1; n<=N-1; n++)
                {
                    sn.i = n; sn.x = n*hx;

                    dx[n-1] = 0.0;
                    if (m>0 && m<M)  dx[n-1] = p_aa_htht__hyhy*(u10[m-1][n] - 2.0*u10[m][n]   + u10[m+1][n]);
                    else if (m == 0) dx[n-1] = p_aa_htht__hyhy*(u10[0][n]   - 2.0*u10[1][n]   + u10[2][n]);
                    else if (m == M) dx[n-1] = p_aa_htht__hyhy*(u10[M-2][n] - 2.0*u10[M-1][n] + u10[M][n]);

                    dx[n-1] += lambda_ht*(u10[m][n] - 0.5*(u10[m][n]-u00[m][n]));
                    dx[n-1] += 2.0*u10[m][n] + (u10[m][n]-u00[m][n]);
                }

                dx[0]   -= m_aa_htht__hxhx * u15[m][0];
                dx[N-2] -= m_aa_htht__hxhx * u15[m][N];

                tomasAlgorithm(ax, bx, cx, dx, rx, N-1);
                for (unsigned int n=1; n<=N-1; n++) u15[m][n] = rx[n-1];
            }
        }

        if (rows1.size() != 0 && rows2.size() == 0)
        {
            //throw std::exception();

            double* U15 = (double *) malloc(sizeof(double)*No);
            for (unsigned int j=0; j<No; j++) U15[j] = 0.0;
            for (unsigned int odj=0; odj<obsPointNodes.size(); odj++)
            {
                const ExtendedSpacePointNode &opn = obsPointNodes[odj];
                U15[opn.id] += u15[opn.j][opn.i] * (opn.w * (hx*hy));
            }

            double *v = new double[Nc];
            for (unsigned int i=0; i<Nc; i++)
            {
                v[i] = 0.0;
                for (unsigned int j=0; j<No; j++)
                {
                    v[i] += mOptParameter.k[i][j] * (U15[j]-mOptParameter.z[i][j]);
                }
            }

            for (unsigned int row=0; row<rows1.size(); row++)
            {
                unsigned int m = rows1.at(row);
                sn.j = m; sn.y = m*hy;
                for (unsigned int n=1; n<=N-1; n++)
                {
                    sn.i = n; sn.x = n*hx;

                    dx[n-1] = 0.0;
                    if (m>0 && m<M)  dx[n-1] = p_aa_htht__hyhy*(u10[m-1][n] - 2.0*u10[m][n]   + u10[m+1][n]);
                    else if (m == 0) dx[n-1] = p_aa_htht__hyhy*(u10[0][n]   - 2.0*u10[1][n]   + u10[2][n]);
                    else if (m == M) dx[n-1] = p_aa_htht__hyhy*(u10[M-2][n] - 2.0*u10[M-1][n] + u10[M][n]);

                    dx[n-1] += lambda_ht*(u10[m][n] - 0.5*(u10[m][n]-u00[m][n]));
                    dx[n-1] += 2.0*u10[m][n] + (u10[m][n]-u00[m][n]);

                    //------------------------------------- Adding delta part -------------------------------------//
                    for (unsigned int cni=0; cni<cntDeltaNodes.size(); cni++)
                    {
                        const ExtendedSpacePointNode &cdn = cntDeltaNodes.at(cni);
                        if (cdn.i == sn.i && cdn.j == sn.j)
                        {
                            dx[n-1] += htht * v[cdn.id] * cdn.w;
                            //for (unsigned int j=0; j<No; j++)
                            //{
                            //    dx[n-1] += htht * mOptParameter.k[cdn.id][j] * (U15[j]-mOptParameter.z[cdn.id][j]) * cdn.w;
                            //}
                        }
                    }
                    //------------------------------------- Adding delta part -------------------------------------//

                    //if (l==2)
                    //{
                    //    for (unsigned int si=0; si<qPointNodes.size(); si++)
                    //    {
                    //        const ExtendedSpacePointNode &qNode = qPointNodes.at(si);
                    //        if (qNode.i == sn.i && qNode.j == sn.j)
                    //        {
                    //            dx[n-1] += (mEquParameter.q[qNode.id] * qNode.w * (2.0/ht)) * htht_h;
                    //        }
                    //    }
                    //}
                }

                dx[0]   -= m_aa_htht__hxhx * u15[m][0];
                dx[N-2] -= m_aa_htht__hxhx * u15[m][N];

                tomasAlgorithm(ax, bx, cx, dx, rx, N-1);
                for (unsigned int n=1; n<=N-1; n++) u15[m][n] = rx[n-1];
            }

            delete [] v;
            free(U15);
        }

        if (rows1.size() != 0 && rows2.size() != 0)
        {
            //throw std::exception();

            unsigned int row1_size = rows1.size()*(N-1);
            double* a1 = (double*) malloc(sizeof(double)*row1_size);
            double* b1 = (double*) malloc(sizeof(double)*row1_size);
            double* c1 = (double*) malloc(sizeof(double)*row1_size);
            double* d1 = (double*) malloc(sizeof(double)*row1_size);
            double* x1 = (double*) malloc(sizeof(double)*row1_size);
            DoubleMatrix w1(row1_size, row1_size, 0.0);

            unsigned int offset = 0;
            for (unsigned int row=0; row<rows1.size(); row++)
            {
                unsigned int m = rows1.at(row);
                sn.j = m; sn.y = m*hy;

                for (unsigned int n=1; n<=N-1; n++)
                {
                    sn.i = n; sn.x = n*hx;

                    d1[offset+(n-1)] = 0.0;
                    if (m>0 && m<M) d1[offset+(n-1)] = p_aa_htht__hyhy*(u10[m-1][n] - 2.0*u10[m][n]   + u10[m+1][n]);
                    else if (m==0)  d1[offset+(n-1)] = p_aa_htht__hyhy*(u10[0][n]   - 2.0*u10[1][n]   + u10[2][n]);
                    else if (m==M)  d1[offset+(n-1)] = p_aa_htht__hyhy*(u10[M-2][n] - 2.0*u10[M-1][n] + u10[M][n]);

                    d1[offset+(n-1)] += lambda_ht*(u10[m][n] - 0.5*(u10[m][n]-u00[m][n]));
                    d1[offset+(n-1)] += 2.0*u10[m][n] + (u10[m][n]-u00[m][n]);

                    a1[offset+(n-1)] = m_aa_htht__hxhx;
                    b1[offset+(n-1)] = p_aa_htht__hxhx___lambda_ht;
                    c1[offset+(n-1)] = m_aa_htht__hxhx;

                    //------------------------------------- Adding delta part -------------------------------------//
                    for (unsigned int cni=0; cni<cntDeltaNodes.size(); cni++)
                    {
                        const ExtendedSpacePointNode &cdn = cntDeltaNodes.at(cni);
                        if (cdn.i == sn.i && cdn.j == sn.j)
                        {
                            for (unsigned int odj=0; odj<obsPointNodes.size(); odj++)
                            {
                                const ExtendedSpacePointNode &opn = obsPointNodes.at(odj);

                                bool found = false;
                                for (unsigned int rs=0; rs<rows1.size(); rs++)
                                {
                                    if (opn.j == rows1[rs])
                                    {
                                        found = true;
                                        w1[offset+(n-1)][rs*(N-1)+(opn.i-1)] -= htht * mOptParameter.k[cdn.id][opn.id] * (opn.w * (hx*hy)) * cdn.w;
                                        break;
                                    }
                                }

                                if (!found)
                                {
                                    d1[offset+(n-1)] += htht * mOptParameter.k[cdn.id][opn.id] * u15[opn.j][opn.i] * (opn.w * (hx*hy)) * cdn.w;
                                }
                            }

                            for (unsigned int j=0; j<No; j++)
                            {
                                d1[offset+(n-1)] -= htht * mOptParameter.k[cdn.id][j] * mOptParameter.z[cdn.id][j] * cdn.w;
                            }
                        }
                    }
                    //------------------------------------- Adding delta part -------------------------------------//
                }

                a1[offset+0]   = 0.0;
                c1[offset+N-2] = 0.0;

                d1[offset+0]   -= m_aa_htht__hxhx * u15[m][0];
                d1[offset+N-2] -= m_aa_htht__hxhx * u15[m][N];

                offset += N-1;
            }

            LinearEquation::func1(a1, b1, c1, d1, w1.data(), x1, row1_size);

            offset = 0;
            for (unsigned int row=0; row<rows1.size(); row++)
            {
                unsigned int m=rows1.at(row);
                for (unsigned int n=1; n<=N-1; n++)
                {
                    u15[m][n] = x1[offset+(n-1)];
                }
                offset += N-1;
            }

            w1.clear();
            free(x1);
            free(d1);
            free(c1);
            free(b1);
            free(a1);
        }

        /**************************************************** x direction apprx ***************************************************/

        /**************************************************** y direction apprx ***************************************************/

        if (cols0.size() != 0)
        {
            for (unsigned int col=0; col<cols0.size(); col++)
            {
                unsigned int n = cols0.at(col);
                sn.i = n; sn.x = n*hx;
                for (unsigned int m=1; m<=M-1; m++)
                {
                    sn.j = m; sn.y = m*hy;

                    dy[m-1] = 0.0;
                    if (n>0 && n<N) dy[m-1] = p_aa_htht__hxhx*(u15[m][n-1] - 2.0*u15[m][n]   + u15[m][n+1]);
                    else if (n==0)  dy[m-1] = p_aa_htht__hxhx*(u15[m][0]   - 2.0*u15[m][1]   + u15[m][2]);
                    else if (n==N)  dy[m-1] = p_aa_htht__hxhx*(u15[m][N-2] - 2.0*u15[m][N-1] + u15[m][N]);

                    dy[m-1] += lambda_ht*(u15[m][n] - 0.5*(u10[m][n]-u00[m][n]));
                    dy[m-1] += 2.0*u15[m][n] + (u10[m][n]-u00[m][n]);
                }

                dy[0]   -= m_aa_htht__hyhy * u20[0][n];
                dy[M-2] -= m_aa_htht__hyhy * u20[M][n];

                tomasAlgorithm(ay, by, cy, dy, ry, M-1);
                for (unsigned int m=1; m<=M-1; m++) u20[m][n] = ry[m-1];
            }
        }

        if (cols1.size() != 0 && cols2.size() == 0)
        {
            //throw std::exception();

            double* U20 = (double *) malloc(sizeof(double)*No);
            for (unsigned int j=0; j<No; j++) U20[j] = 0.0;
            for (unsigned int odj=0; odj<obsPointNodes.size(); odj++)
            {
                const ExtendedSpacePointNode &opn = obsPointNodes[odj];
                U20[opn.id] += u20[opn.j][opn.i] * (opn.w * (hx*hy));
            }

            double *v = new double[Nc];
            for (unsigned int i=0; i<Nc; i++)
            {
                v[i] = 0.0;
                for (unsigned int j=0; j<No; j++)
                {
                    v[i] += mOptParameter.k[i][j] * (U20[j]-mOptParameter.z[i][j]);
                }
            }

            for (unsigned int col=0; col<cols1.size(); col++)
            {
                unsigned int n = cols1.at(col);
                sn.i = n; sn.x = n*hx;
                for (unsigned int m=1; m<=M-1; m++)
                {
                    sn.j = m; sn.y = m*hy;

                    dy[m-1] = 0.0;
                    if (n>0 && n<N) dy[m-1] = p_aa_htht__hxhx*(u15[m][n-1] - 2.0*u15[m][n]   + u15[m][n+1]);
                    else if (n==0)  dy[m-1] = p_aa_htht__hxhx*(u15[m][0]   - 2.0*u15[m][1]   + u15[m][2]);
                    else if (n==N)  dy[m-1] = p_aa_htht__hxhx*(u15[m][N-2] - 2.0*u15[m][N-1] + u15[m][N]);

                    dy[m-1] += lambda_ht*(u15[m][n] - 0.5*(u10[m][n]-u00[m][n]));
                    dy[m-1] += 2.0*u15[m][n] + (u10[m][n]-u00[m][n]);

                    //------------------------------------- Adding delta part -------------------------------------//
                    for (unsigned int cni=0; cni<cntDeltaNodes.size(); cni++)
                    {
                        const ExtendedSpacePointNode &cdn = cntDeltaNodes.at(cni);
                        if (cdn.i == sn.i && cdn.j == sn.j)
                        {
                            dy[m-1] += htht * v[cdn.id] * cdn.w;
                            //for (unsigned int j=0; j<No; j++)
                            //{
                            //    if (cols1[col]==48 && l==3) printf("%d %d %f %d %d\n", cdn.id, j, cdn.w, cdn.i, cdn.j);
                            //    dy[m-1] += htht * mOptParameter.k[cdn.id][j] * (U20[j]-mOptParameter.z[cdn.id][j]) * cdn.w;
                            //}
                        }
                    }
                    //------------------------------------- Adding delta part -------------------------------------//
                }

                dy[0]   -= m_aa_htht__hyhy * u20[0][n];
                dy[M-2] -= m_aa_htht__hyhy * u20[M][n];

                tomasAlgorithm(ay, by, cy, dy, ry, M-1);
                for (unsigned int m=1; m<=M-1; m++) u20[m][n] = ry[m-1];
            }

            delete [] v;
            free(U20);
        }

        if (cols1.size() != 0 && cols2.size() != 0)
        {
            //throw std::exception();

            unsigned int cols1_size = cols1.size()*(M-1);
            double* a2 = (double*) malloc(sizeof(double)*cols1_size);
            double* b2 = (double*) malloc(sizeof(double)*cols1_size);
            double* c2 = (double*) malloc(sizeof(double)*cols1_size);
            double* d2 = (double*) malloc(sizeof(double)*cols1_size);
            double* x2 = (double*) malloc(sizeof(double)*cols1_size);
            DoubleMatrix w2(cols1_size, cols1_size, 0.0);

            unsigned int offset = 0;
            for (unsigned int col=0; col<cols1.size(); col++)
            {
                unsigned int n = cols1.at(col);
                sn.i = n; sn.x = n*hx;

                for (unsigned int m=1; m<=M-1; m++)
                {
                    sn.j = m; sn.y = m*hy;

                    d2[offset+(m-1)] = 0.0;
                    if (n>0 && n<N) d2[offset+(m-1)] = p_aa_htht__hxhx*(u15[m][n-1] - 2.0*u15[m][n]   + u15[m][n+1]);
                    else if (n==0)  d2[offset+(m-1)] = p_aa_htht__hxhx*(u15[m][0]   - 2.0*u15[m][1]   + u15[m][2]);
                    else if (n==N)  d2[offset+(m-1)] = p_aa_htht__hxhx*(u15[m][N-2] - 2.0*u15[m][N-1] + u15[m][N]);

                    d2[offset+(m-1)] += lambda_ht*(u15[m][n] - 0.5*(u10[m][n]-u00[m][n]));
                    d2[offset+(m-1)] += 2.0*u15[m][n] + (u10[m][n]-u00[m][n]);

                    a2[offset+(m-1)] = m_aa_htht__hyhy;
                    b2[offset+(m-1)] = p_aa_htht__hyhy___lambda_ht;
                    c2[offset+(m-1)] = m_aa_htht__hyhy;

                    //------------------------------------- Adding delta part -------------------------------------//
                    for (unsigned int cni=0; cni<cntDeltaNodes.size(); cni++)
                    {
                        const ExtendedSpacePointNode &cdn = cntDeltaNodes.at(cni);
                        if (cdn.i == sn.i && cdn.j == sn.j)
                        {
                            for (unsigned int onj=0; onj<obsPointNodes.size(); onj++)
                            {
                                const ExtendedSpacePointNode &opn = obsPointNodes.at(onj);

                                bool found = false;
                                for (unsigned int cs=0; cs<cols1.size(); cs++)
                                {
                                    if (opn.i == cols1[cs])
                                    {
                                        found = true;
                                        w2[offset+(m-1)][cs*(M-1)+(opn.j-0)] -= htht * mOptParameter.k[cdn.id][opn.id] * (opn.w * (hx*hy)) * cdn.w;
                                        break;
                                    }
                                }

                                if (!found)
                                {
                                    d2[offset+(m-1)] += htht * mOptParameter.k[cdn.id][opn.id] * u20[opn.j][opn.i] * (opn.w * (hx*hy)) * cdn.w;
                                }
                            }
                            for (unsigned int j=0; j<No; j++)
                            {
                                d2[offset+(m-1)] -= htht * mOptParameter.k[cdn.id][j] * mOptParameter.z[cdn.id][j] * cdn.w;
                            }
                        }
                    }
                    //------------------------------------- Adding delta part -------------------------------------//
                }

                a2[offset+0]   = 0.0;
                c2[offset+M-2] = 0.0;

                d2[offset+0]   -= m_aa_htht__hyhy * u20[0][n];
                d2[offset+M-2] -= m_aa_htht__hyhy * u20[M][n];

                offset += M-1;
            }

            LinearEquation::func1(a2, b2, c2, d2, w2.data(), x2, cols1.size()*(M-1));

            offset = 0;
            for (unsigned int col=0; col<cols1.size(); col++)
            {
                unsigned int n=cols1.at(col);
                for (unsigned int m=1; m<=M-1; m++)
                {
                    u20[m][n] = x2[offset+(m-1)];
                }
                offset += M-1;
            }

            w2.clear();
            free(x2);
            free(d2);
            free(c2);
            free(b2);
            free(a2);
        }

        /**************************************************** y direction apprx ***************************************************/

        for (unsigned int m=0; m<=M; m++)
        {
            for (unsigned int n=0; n<=N; n++)
            {
                u00[m][n] = u10[m][n];
                u10[m][n] = u20[m][n];
            }
        }

        if (use == true) f_add2Info(u20, u_info, obsPointNodes, l, hx, hy);

        f_layerInfo(u20, l);

        if (L == l)
        {
            for (unsigned int m=0; m<=M; m++)
            {
                for (unsigned int n=0; n<=N; n++)
                {
                    u[l-L][m][n] = u20[m][n];
                }
            }
        }

        if ( L+1 <= l && l <= LLD )
        {
            for (unsigned int m=0; m<=M; m++)
            {
                for (unsigned int n=0; n<=N; n++)
                {
                    u[2*(l-L)-1][m][n] = u15[m][n];
                    u[2*(l-L)+0][m][n] = u20[m][n];
                }
            }
        }
    }

    free(rx);
    free(dx);
    free(cx);
    free(bx);
    free(ax);

    free(ry);
    free(dy);
    free(cy);
    free(by);
    free(ay);

    rows0.clear();
    rows1.clear();
    rows2.clear();

    cols0.clear();
    cols1.clear();
    cols2.clear();

    qPointNodes.clear();
    obsPointNodes.clear();
    cntDeltaNodes.clear();

    u00.clear();
    u10.clear();
    u15.clear();
    u20.clear();
}

void Problem2HNDirichlet::f_initialLayers1(DoubleMatrix &u00, DoubleMatrix &u10, spif_vector &u_info, bool use,
                                           espn_vector &obsPointNodes, espn_vector &cntDeltaNodes UNUSED_PARAM,
                                           espn_vector &qPointNodes UNUSED_PARAM, unsigned int N, unsigned int M, double hx, double hy, double ht, double aa__hxhx, double aa__hyhy, double lambda) const
{
    /************************************************************************/
    SpaceNodePDE sn;
    for (unsigned int m=0; m<=M; m++)
    {
        sn.j = m; sn.y = m*hy;
        for (unsigned int n=0; n<=N; n++)
        {
            sn.i = n; sn.x = n*hx;
            u00[m][n] = f_initial1(sn);
        }
    }

    /************************************************************************/
    TimeNodePDE tn10; tn10.i = 1; tn10.t = ht;

    SpaceNodePDE sn0, sn1;
    sn0.i = 0; sn0.x = 0.0;
    sn1.i = N; sn1.x = N*hx;
    for (unsigned int m=0; m<=M; m++)
    {
        sn0.j = m; sn0.y = m*hy; u10[m][0] = f_boundary(sn0, tn10);
        sn1.j = m; sn1.y = m*hy; u10[m][N] = f_boundary(sn1, tn10);
    }

    sn0.j = 0;  sn0.y = 0.0;
    sn1.j = M; sn1.y = M*hy;
    for (unsigned int n=0; n<=N; n++)
    {
        sn0.i = n; sn0.x = n*hx; u10[0][n] = f_boundary(sn0, tn10);
        sn1.i = n; sn1.x = n*hx; u10[M][n] = f_boundary(sn1, tn10);
    }

    /************************************************************************/

    for (unsigned int m=1; m<=M-1; m++)
    {
        sn.j = m; sn.y = m*hy;
        for (unsigned int n=1; n<=N-1; n++)
        {
            sn.i = n; sn.x = n*hx;

            double sum = 0.0;
            sum += aa__hxhx*(u00[m][n-1]-2.0*u00[m][n]+u00[m][n+1]);
            sum += aa__hyhy*(u00[m-1][n]-2.0*u00[m][n]+u00[m+1][n]);
            sum -= lambda*f_initial2(sn);

            double Q = 0.0;
            for (unsigned int si=0; si<qPointNodes.size(); si++)
            {
                const ExtendedSpacePointNode &qNode = qPointNodes.at(si);
                if (qNode.i == sn.i && qNode.j == sn.j)
                {
                    Q += mEquParameter.q[qNode.id] * qNode.w;
                }
            }

            u10[m][n] = u00[m][n] + f_initial2(sn)*ht + 0.5*ht*ht*sum;

            u10[m][n] += Q*ht*0.5;
        }
    }

    if (use == true) f_add2Info(u00, u_info, obsPointNodes, 0, hx, hy);
    f_layerInfo(u00, 0);

    if (use == true) f_add2Info(u10, u_info, obsPointNodes, 1, hx, hy);
    f_layerInfo(u10, 1);
}

double Problem2HNDirichlet::f_initial1(const SpaceNodePDE &sn UNUSED_PARAM) const
{
    return 0.0;
}

double Problem2HNDirichlet::f_initial2(const SpaceNodePDE &) const
{
    return 0.0;
}

double Problem2HNDirichlet::f_boundary(const SpaceNodePDE &sn UNUSED_PARAM, const TimeNodePDE &tn UNUSED_PARAM, BoundaryType) const
{
    return 0.0;
}

void Problem2HNDirichlet::f_findRowsCols(uint_vector &rows0, uint_vector &rows1, uint_vector &rows2, uint_vector &cols0, uint_vector &cols1, uint_vector &cols2,
                                         espn_vector &obsPointNodes, espn_vector &cntDeltaNodes, unsigned int N, unsigned int M) const
{
    for (unsigned int m=1; m<=M-1; m++)
    {
        bool found1 = false;
        bool found2 = false;
        for (unsigned int i=0; i<cntDeltaNodes.size(); i++)
        {
            const ExtendedSpacePointNode &cdn = cntDeltaNodes.at(i);

            if (cdn.j == m)
            {
                found1 = true;
                for (unsigned int j=0; j<obsPointNodes.size(); j++)
                {
                    const ExtendedSpacePointNode &opn = obsPointNodes.at(j);
                    if (opn.j == m)
                    {
                        found2 = true;
                        break;
                    }
                }
                break;
            }
        }
        if (found1 == false && found2 == false) if(std::find(rows0.begin(), rows0.end(), m) == rows0.end()) rows0.push_back(m);
        if (found1 == true  && found2 == true)  if(std::find(rows2.begin(), rows2.end(), m) == rows2.end()) rows2.push_back(m);
        if (found1 == true)                     if(std::find(rows1.begin(), rows1.end(), m) == rows1.end()) rows1.push_back(m);
    }

    for (unsigned int n=1; n<=N-1; n++)
    {
        bool found1 = false;
        bool found2 = false;
        for (unsigned int i=0; i<cntDeltaNodes.size(); i++)
        {
            const ExtendedSpacePointNode &cdn = cntDeltaNodes.at(i);
            if (cdn.i == n)
            {
                found1 = true;
                for (unsigned int j=0; j<obsPointNodes.size(); j++)
                {
                    const ExtendedSpacePointNode &opn = obsPointNodes.at(j);
                    if (opn.i == n)
                    {
                        found2 = true;
                        break;
                    }
                }
                break;
            }
        }
        if (found1 == false && found2 == false) if(std::find(cols0.begin(), cols0.end(), n) == cols0.end()) cols0.push_back(n);
        if (found1 == true  && found2 == true)  if(std::find(cols2.begin(), cols2.end(), n) == cols2.end()) cols2.push_back(n);
        if (found1 == true)                     if(std::find(cols1.begin(), cols1.end(), n) == cols1.end()) cols1.push_back(n);
    }
}

void Problem2HNDirichlet::f_initialLayers(DoubleMatrix &u00, DoubleMatrix &u05, DoubleMatrix &u10, spif_vector &u_info, bool use,
                                          espn_vector &obsPointNodes, espn_vector &cntDeltaNodes UNUSED_PARAM,
                                          espn_vector &qPointNodes UNUSED_PARAM, unsigned int N, unsigned int M, double hx, double hy, double ht, double aa__hxhx, double aa__hyhy, double lambda) const
{
    SpaceNodePDE sn;
    for (unsigned int m=0; m<=M; m++)
    {
        sn.j = m; sn.y = m*hy;
        for (unsigned int n=0; n<=N; n++)
        {
            sn.i = n; sn.x = n*hx;
            u00[m][n] = f_initial1(sn);
        }
    }

    TimeNodePDE tn05; tn05.i = 1; tn05.t = ht*0.5;
    TimeNodePDE tn10; tn10.i = 1; tn10.t = ht;

    SpaceNodePDE sn0, sn1;
    sn0.i = 0; sn0.x = 0.0;
    sn1.i = N; sn1.x = N*hx;
    for (unsigned int m=0; m<=M; m++)
    {
        sn0.j = m; sn0.y = m*hy; u05[m][0] = f_boundary(sn0, tn05); u10[m][0] = f_boundary(sn0, tn10);
        sn1.j = m; sn1.y = m*hy; u05[m][N] = f_boundary(sn1, tn05); u10[m][N] = f_boundary(sn1, tn10);
    }

    sn0.j = 0;  sn0.y = 0.0;
    sn1.j = M; sn1.y = M*hy;
    for (unsigned int n=0; n<=N; n++)
    {
        sn0.i = n; sn0.x = n*hx; u05[0][n] = f_boundary(sn0, tn05); u10[0][n] = f_boundary(sn0, tn10);
        sn1.i = n; sn1.x = n*hx; u05[M][n] = f_boundary(sn1, tn05); u10[M][n] = f_boundary(sn1, tn10);
    }

    for (unsigned int m=1; m<=M-1; m++)
    {
        sn.j = m; sn.y = m*hy;
        for (unsigned int n=1; n<=N-1; n++)
        {
            sn.i = n; sn.x = n*hx;

            double sum = 0.0;
            sum += aa__hxhx*(u00[m][n-1]-2.0*u00[m][n]+u00[m][n+1]);
            sum += aa__hyhy*(u00[m-1][n]-2.0*u00[m][n]+u00[m+1][n]);
            sum -= lambda*f_initial2(sn);

            u05[m][n] = u00[m][n] + f_initial2(sn)*ht*0.5 + 0.125*ht*ht*sum;
            u10[m][n] = u00[m][n] + f_initial2(sn)*ht     + 0.500*ht*ht*sum;
        }
    }

    /*
    for (unsigned int m=1; m<=M-1; m++)
    {
        sn.j = m; sn.y = m*hy;
        for (unsigned int n=1; n<=N-1; n++)
        {
            sn.i = n; sn.x = n*hx;

            for (unsigned int si=0; si<qPointNodes.size(); si++)
            {
                const ExtendedSpacePointNode &qNode = qPointNodes.at(si);
                if (qNode.i == sn.i && qNode.j == sn.j)
                {
                    u10[m][n] += (mEquParameter.q[qNode.id] * qNode.w)*(0.001*0.001);
                }
            }
        }
    }
    */

    if (use == true) f_add2Info(u00, u_info, obsPointNodes, 0, hx, hy);
    f_layerInfo(u00, 0);

    if (use == true) f_add2Info(u10, u_info, obsPointNodes, 1, hx, hy);
    f_layerInfo(u10, 1);
}

void Problem2HNDirichlet::f_borderLayer(DoubleMatrix &u, DoubleMatrix &um5, unsigned int ln) const
{
    const Dimension dimX = spaceDimension(Dimension::DimensionX);
    const Dimension dimY = spaceDimension(Dimension::DimensionY);
    const Dimension time = timeDimension();

    const unsigned int N = dimX.sizeN();
    const unsigned int M = dimY.sizeN();

    const double hx = dimX.step();
    const double hy = dimY.step();
    const double ht = time.step();

    TimeNodePDE tn00; tn00.i = ln; tn00.t = ln*ht;
    TimeNodePDE tnm5; tnm5.i = ln; tnm5.t = ln*ht-0.5*ht;

    SpaceNodePDE sn0;
    SpaceNodePDE sn1;

    sn0.i = 0; sn0.x = 0.0;
    sn1.i = N; sn1.x = hx*N;
    for (unsigned int m=0; m<=M; m++)
    {
        sn0.j = m; sn0.y = m*hy; um5[m][0] = f_boundary(sn0, tnm5); u[m][0] = f_boundary(sn0, tn00);
        sn1.j = m; sn1.y = m*hy; um5[m][N] = f_boundary(sn1, tnm5); u[m][N] = f_boundary(sn1, tn00);
    }

    sn0.j = 0; sn0.y = 0.0;
    sn1.j = M; sn1.y = hy*M;
    for (unsigned int n=0; n<=N; n++)
    {
        sn0.i = n; sn0.x = n*hx; um5[0][n] = f_boundary(sn0, tnm5); u[0][n] = f_boundary(sn0, tn00);
        sn1.i = n; sn1.x = n*hx; um5[M][n] = f_boundary(sn1, tnm5); u[M][n] = f_boundary(sn1, tn00);
    }
}

void Problem2HNDirichlet::f_prepareInfo(unsigned int No, const std::vector<SpacePoint> &points, spif_vector &u_info,
                                        unsigned int LLD, const Dimension &dimX, const Dimension &dimY) const
{
    u_info.resize(No);
    for (unsigned int j=0; j<No; j++)
    {
        SpacePointInfo &inf = u_info[j];
        const SpacePoint &sp = points[j];
        inf.id = j;
        inf.x = sp.x;
        inf.y = sp.y;

        inf.i = (unsigned int)( round( inf.x * dimX.sizeN() ) );
        inf.j = (unsigned int)( round( inf.y * dimY.sizeN() ) );

        inf.layerNumber = LLD+1;
        inf.u  = new double[LLD+1];
        inf.ux = new double[LLD+1];
        inf.uy = new double[LLD+1];
    }
}

void Problem2HNDirichlet::f_add2Info(const DoubleMatrix &u, spif_vector &u_info, const espn_vector &obsPointNodes, unsigned int ln, double hx, double hy, int method) const
{
    if (method == 1)
    {
        for (unsigned int j=0; j<u_info.size(); j++) u_info[j].u[ln] = u_info[j].ux[ln] = u_info[j].uy[ln] = 0.0;

        for (unsigned int opj=0; opj<obsPointNodes.size(); opj++)
        {
            const ExtendedSpacePointNode &opn = obsPointNodes[opj];
            const SpacePointInfo &ui = u_info[opn.id];

            ui.u[ln] += u[opn.j][opn.i] * (opn.w * (hx*hy));

            ui.ux[ln] = (u[ui.j][ui.i+1] - u[ui.j][ui.i-1])/(2.0*hx);
            ui.uy[ln] = (u[ui.j+1][ui.i] - u[ui.j-1][ui.i])/(2.0*hy);

            //ui.ux[ln] = (u[ui.j][ui.i-2]-8.0*u[ui.j][ui.i-1]+8.0*u[ui.j][ui.i+1]-u[ui.j][ui.i+2])/(12.0*hx);
            //ui.uy[ln] = (u[ui.j-2][ui.i]-8.0*u[ui.j-1][ui.i]+8.0*u[ui.j+1][ui.i]-u[ui.j+2][ui.i])/(12.0*hy);
        }
    }

    if (method == 2)
    {
        for (unsigned int j=0; j<u_info.size(); j++) u_info[j].u[ln] = u_info[j].ux[ln] = u_info[j].uy[ln] = 0.0;

        for (unsigned int opj=0; opj<obsPointNodes.size(); opj++)
        {
            const ExtendedSpacePointNode &opn = obsPointNodes[opj];
            const SpacePointInfo &ui = u_info[opn.id];

            ui.u[ln] += u[opn.j][opn.i] * (opn.w * (hx*hy));

            ui.ux[ln] = (u[ui.j][ui.i+1] - u[ui.j][ui.i-1])/(2.0*hx);
            ui.uy[ln] = (u[ui.j+1][ui.i] - u[ui.j-1][ui.i])/(2.0*hy);

            //ui.ux[ln] = (u[ui.j][ui.i-2]-8.0*u[ui.j][ui.i-1]+8.0*u[ui.j][ui.i+1]-u[ui.j][ui.i+2])/(12.0*hx);
            //ui.uy[ln] = (u[ui.j-2][ui.i]-8.0*u[ui.j-1][ui.i]+8.0*u[ui.j+1][ui.i]-u[ui.j+2][ui.i])/(12.0*hy);
        }
    }

    if (method == 4)
    {
        for (unsigned int j=0; j<u_info.size(); j++) u_info[j].u[ln] = u_info[j].ux[ln] = u_info[j].uy[ln] = 0.0;

        for (unsigned int opj=0; opj<obsPointNodes.size(); opj++)
        {
            const ExtendedSpacePointNode &opn = obsPointNodes[opj];
            const SpacePointInfo &ui = u_info[opn.id];

            ui.u[ln] += u[opn.j][opn.i] * (opn.w * (hx*hy));

            ui.ux[ln] = (u[ui.j][ui.i+1] - u[ui.j][ui.i-1])/(2.0*hx);
            ui.uy[ln] = (u[ui.j+1][ui.i] - u[ui.j-1][ui.i])/(2.0*hy);

            //ui.ux[ln] = (u[ui.j][ui.i-2]-8.0*u[ui.j][ui.i-1]+8.0*u[ui.j][ui.i+1]-u[ui.j][ui.i+2])/(12.0*hx);
            //ui.uy[ln] = (u[ui.j-2][ui.i]-8.0*u[ui.j-1][ui.i]+8.0*u[ui.j+1][ui.i]-u[ui.j+2][ui.i])/(12.0*hy);
        }
    }
}

//forward -------------------------------------

// backward -----------------------------------
void Problem2HNDirichlet::solveBackwardIBVP0(const std::vector<DoubleMatrix> &u, spif_vector &p_info, bool use, const spif_vector &u_info) const
{
    //puts("-void Problem2HNDirichlet::solveBackwardIBVP(const std::vector<DoubleMatrix> &u, spif_vector &p_info, bool use, const spif_vector &u_info) const");
    const Dimension dimX = spaceDimension(Dimension::DimensionX);
    const Dimension dimY = spaceDimension(Dimension::DimensionY);
    const Dimension time = timeDimension();

    const unsigned int N = dimX.sizeN();
    const unsigned int M = dimY.sizeN();
    const unsigned int L = time.sizeN();
    const unsigned int LLD = L+LD;

    const double hx = dimX.step();
    const double hy = dimY.step();
    const double ht = time.step();

    const double a = mEquParameter.a;
    const double lambda = mEquParameter.lambda;
    const unsigned int No = mEquParameter.No;
    const unsigned int Nc = mEquParameter.Nc;

    const double m_aa_htht__hxhx = -(a*a*ht*ht)/(hx*hx);
    const double p_aa_htht__hxhx___lambda_ht = +2.0 + 2.0*(a*a*ht*ht)/(hx*hx) + 3.0*lambda*ht;
    const double p_aa_htht__hxhx = +(a*a*ht*ht)/(hx*hx);

    const double m_aa_htht__hyhy = -(a*a*ht*ht)/(hy*hy);
    const double p_aa_htht__hyhy___lambda_ht = +2.0 + 2.0*(a*a*ht*ht)/(hy*hy) + 3.0*lambda*ht;
    const double p_aa_htht__hyhy = +(a*a*ht*ht)/(hy*hy);

    const double htht = ht*ht;
    const double lambda_ht = lambda*ht;

    const double aa__hxhx = (a*a)/(hx*hx);
    const double aa__hyhy = (a*a)/(hy*hy);

    DoubleMatrix p00(M+1, N+1);
    DoubleMatrix p05(M+1, N+1);
    DoubleMatrix p10(M+1, N+1);
    DoubleMatrix p15(M+1, N+1);
    DoubleMatrix p20(M+1, N+1);

    //--------------------------------------------------------------------------------------------//
    espn_vector obsDeltaNodes, cntPointNodes;
    for (unsigned int j=0; j<No; j++) distributeDelta0(mOptParameter.xi[j], j, obsDeltaNodes, dimX, dimY, 4, 4);
    for (unsigned int i=0; i<Nc; i++) distributeDelta0(mOptParameter.eta[i], i, cntPointNodes, dimX, dimY, 4, 4);

    //----------------------------------------------------------------------------------------------//
    uint_vector rows0, rows1, rows2, cols0, cols1, cols2;
    b_findRowsCols(rows0, rows1, rows2, cols0, cols1, cols2, cntPointNodes, obsDeltaNodes, N, M);
    //-------------------------------------------- info --------------------------------------------//
    if (use == true) b_prepareInfo(Nc, mOptParameter.eta, p_info, LLD, dimX, dimY);
    //-------------------------------------------- info --------------------------------------------//

    //------------------------------------- initial conditions -------------------------------------//
    b_initialLayers(p00, p05, p10, p_info, use, cntPointNodes, obsDeltaNodes, N, M, hx, hy, ht, aa__hxhx, aa__hyhy, lambda);
    //------------------------------------- initial conditions -------------------------------------//

    double *ax = (double *) malloc(sizeof(double)*(N-1)); for (unsigned int n=1; n<=N-1; n++) ax[n-1] = m_aa_htht__hxhx;
    double *bx = (double *) malloc(sizeof(double)*(N-1)); for (unsigned int n=1; n<=N-1; n++) bx[n-1] = p_aa_htht__hxhx___lambda_ht;
    double *cx = (double *) malloc(sizeof(double)*(N-1)); for (unsigned int n=1; n<=N-1; n++) cx[n-1] = m_aa_htht__hxhx;
    double *dx = (double *) malloc(sizeof(double)*(N-1));
    double *rx = (double *) malloc(sizeof(double)*(N-1));
    ax[0] = cx[N-2] = 0.0;

    double *ay = (double *) malloc(sizeof(double)*(M-1)); for (unsigned int m=1; m<=M-1; m++) ay[m-1] = m_aa_htht__hyhy;
    double *by = (double *) malloc(sizeof(double)*(M-1)); for (unsigned int m=1; m<=M-1; m++) by[m-1] = p_aa_htht__hyhy___lambda_ht;
    double *cy = (double *) malloc(sizeof(double)*(M-1)); for (unsigned int m=1; m<=M-1; m++) cy[m-1] = m_aa_htht__hyhy;
    double *dy = (double *) malloc(sizeof(double)*(M-1));
    double *ry = (double *) malloc(sizeof(double)*(M-1));
    ay[0] = cy[M-2] = 0.0;

    SpaceNodePDE sn;

    for (unsigned int l1=2,l=LLD-l1; l1<=LLD; l1++,l=LLD-l1)
    {
        TimeNodePDE tn20; tn20.i = l; tn20.t = l*ht;
        TimeNodePDE tn15; tn15.i = l; tn15.t = l*ht+0.5*ht;

        /**************************************************** border conditions ***************************************************/

        SpaceNodePDE sn0;
        SpaceNodePDE sn1;

        sn0.i = 0; sn0.x = 0.0;
        sn1.i = N; sn1.x = hx*N;
        for (unsigned int m=0; m<=M; m++)
        {
            sn0.j = m; sn0.y = m*hy; p15[m][0] = b_boundary(sn0, tn15); p20[m][0] = b_boundary(sn0, tn20);
            sn1.j = m; sn1.y = m*hy; p15[m][N] = b_boundary(sn1, tn15); p20[m][N] = b_boundary(sn1, tn20);
        }

        sn0.j = 0; sn0.y = 0.0;
        sn1.j = M; sn1.y = hy*M;
        for (unsigned int n=0; n<=N; n++)
        {
            sn0.i = n; sn0.x = n*hx; p15[0][n] = b_boundary(sn0, tn15); p20[0][n] = b_boundary(sn0, tn20);
            sn1.i = n; sn1.x = n*hx; p15[M][n] = b_boundary(sn1, tn15); p20[M][n] = b_boundary(sn1, tn20);
        }

        /**************************************************** border conditions ***************************************************/

        /**************************************************** x direction apprx ***************************************************/

        if (rows0.size() != 0)
        {
            for (unsigned int row=0; row<rows0.size(); row++)
            {
                unsigned int m = rows0.at(row);
                sn.j = m; sn.y = m*hy;
                for (unsigned int n=1; n<=N-1; n++)
                {
                    sn.i = n; sn.x = n*hx;

                    dx[n-1] = 0.0;
                    if (m>0 && m<M)  dx[n-1] = p_aa_htht__hyhy*(p10[m-1][n] - 2.0*p10[m][n]   + p10[m+1][n]);
                    else if (m == 0) dx[n-1] = p_aa_htht__hyhy*(p10[0][n]   - 2.0*p10[1][n]   + p10[2][n]);
                    else if (m == M) dx[n-1] = p_aa_htht__hyhy*(p10[M-2][n] - 2.0*p10[M-1][n] + p10[M][n]);

                    dx[n-1] += (p10[m][n]-p00[m][n]) + 2.0*p10[m][n];
                    dx[n-1] += lambda_ht*(4.0*p10[m][n]-p05[m][n]);

                    if (L <= l && l <= LLD) dx[n-1] -= 2.0*htht*(u[2*(l-L)+1][m][n] - V0[m][n]);
                }

                dx[0]   -= m_aa_htht__hxhx * p15[m][0];
                dx[N-2] -= m_aa_htht__hxhx * p15[m][N];

                tomasAlgorithm(ax, bx, cx, dx, rx, N-1);
                for (unsigned int n=1; n<=N-1; n++) p15[m][n] = rx[n-1];
            }
        }

        if (rows1.size() != 0 && rows2.size() == 0)
        {
            //throw std::exception();

            for (unsigned int row=0; row<rows1.size(); row++)
            {
                unsigned int m = rows1.at(row);
                sn.j = m; sn.y = m*hy;
                for (unsigned int n=1; n<=N-1; n++)
                {
                    sn.i = n; sn.x = n*hx;

                    dx[n-1] = 0.0;
                    if (m>0 && m<M)  dx[n-1] = p_aa_htht__hyhy*(p10[m-1][n] - 2.0*p10[m][n]   + p10[m+1][n]);
                    else if (m == 0) dx[n-1] = p_aa_htht__hyhy*(p10[0][n]   - 2.0*p10[1][n]   + p10[2][n]);
                    else if (m == M) dx[n-1] = p_aa_htht__hyhy*(p10[M-2][n] - 2.0*p10[M-1][n] + p10[M][n]);

                    dx[n-1] += (p10[m][n]-p00[m][n]) + 2.0*p10[m][n];
                    dx[n-1] += lambda_ht*(4.0*p10[m][n]-p05[m][n]);

                    if (L <= l && l <= LLD) dx[n-1] -= 2.0*htht*(u[2*(l-L)+1][m][n] - V0[m][n]);

                    //------------------------------------- Adding delta part -------------------------------------//
                    for (unsigned int onj=0; onj<obsDeltaNodes.size(); onj++)
                    {
                        const ExtendedSpacePointNode &odn = obsDeltaNodes.at(onj);
                        if (odn.i == sn.i && odn.j == sn.j)
                        {
                            for (unsigned int cni=0; cni<cntPointNodes.size(); cni++)
                            {
                                const ExtendedSpacePointNode &cpn = cntPointNodes[cni];
                                dx[n-1] += htht * mOptParameter.k[cpn.id][odn.id] * p15[cpn.j][cpn.i] * (cpn.w * (hx*hy)) * odn.w;
                            }

                            for (unsigned int i=0; i<Nc; i++)
                            {
                                dx[n-1] += 2.0*r * htht * mOptParameter.k[i][odn.id] * gpi(i, l, u_info, mOptParameter)*sgn(g0i(i, l, u_info, mOptParameter)) * odn.w;
                            }
                        }
                    }
                    //------------------------------------- Adding delta part -------------------------------------//
                }

                dx[0]   -= m_aa_htht__hxhx * p15[m][0];
                dx[N-2] -= m_aa_htht__hxhx * p15[m][N];

                tomasAlgorithm(ax, bx, cx, dx, rx, N-1);
                for (unsigned int n=1; n<=N-1; n++) p15[m][n] = rx[n-1];
            }
        }

        if (rows1.size() != 0 && rows2.size() != 0)
        {
            //throw std::exception();

            unsigned int row1_size = rows1.size()*(N-1);
            double* a1 = (double*) malloc(sizeof(double)*row1_size);
            double* b1 = (double*) malloc(sizeof(double)*row1_size);
            double* c1 = (double*) malloc(sizeof(double)*row1_size);
            double* d1 = (double*) malloc(sizeof(double)*row1_size);
            double* x1 = (double*) malloc(sizeof(double)*row1_size);
            DoubleMatrix w1(row1_size, row1_size, 0.0);

            unsigned int offset = 0;
            for (unsigned int row=0; row<rows1.size(); row++)
            {
                unsigned int m = rows1.at(row);
                sn.j = m; sn.y = m*hy;

                for (unsigned int n=1; n<=N-1; n++)
                {
                    sn.i = n; sn.x = n*hx;

                    d1[offset+(n-1)] = 0.0;
                    if (m>0 && m<M)  d1[offset+(n-1)] = p_aa_htht__hyhy*(p10[m-1][n] - 2.0*p10[m][n]   + p10[m+1][n]);
                    else if (m == 0) d1[offset+(n-1)] = p_aa_htht__hyhy*(p10[0][n]   - 2.0*p10[1][n]   + p10[2][n]);
                    else if (m == M) d1[offset+(n-1)] = p_aa_htht__hyhy*(p10[M-2][n] - 2.0*p10[M-1][n] + p10[M][n]);

                    d1[offset+(n-1)] += (p10[m][n]-p00[m][n]) + 2.0*p10[m][n];
                    d1[offset+(n-1)] += lambda_ht*(4.0*p10[m][n]-p05[m][n]);

                    if (L <= l && l <= LLD) d1[offset+(n-1)] -= 2.0*htht*(u[2*(l-L)+1][m][n] - V0[m][n]);

                    a1[offset+(n-1)] = m_aa_htht__hxhx;
                    b1[offset+(n-1)] = p_aa_htht__hxhx___lambda_ht;
                    c1[offset+(n-1)] = m_aa_htht__hxhx;

                    //------------------------------------- Adding delta part -------------------------------------//
                    for (unsigned int onj=0; onj<obsDeltaNodes.size(); onj++)
                    {
                        const ExtendedSpacePointNode &odn = obsDeltaNodes.at(onj);
                        if (odn.i == sn.i && odn.j == sn.j)
                        {
                            for (unsigned int cni=0; cni<cntPointNodes.size(); cni++)
                            {
                                const ExtendedSpacePointNode &cpn = cntPointNodes[cni];

                                bool found = false;
                                for (unsigned int rs=0; rs<rows1.size(); rs++)
                                {
                                    if (cpn.j == rows1[rs])
                                    {
                                        found = true;
                                        w1[offset+(n-1)][rs*(N-1)+(cpn.i-1)] -= htht * mOptParameter.k[cpn.id][odn.id] * (cpn.w * (hx*hy)) * odn.w;
                                        break;
                                    }
                                }

                                if (!found)
                                {
                                    d1[offset+(n-1)] += htht * mOptParameter.k[cpn.id][odn.id] * p15[cpn.j][cpn.i] * (cpn.w * (hx*hy)) * odn.w;
                                }
                            }

                            for (unsigned int i=0; i<Nc; i++)
                            {
                                d1[offset+(n-1)] += 2.0 * r * ht *  mOptParameter.k[i][odn.id] * gpi(i, l, u_info, mOptParameter)*sgn(g0i(i, l, u_info, mOptParameter)) * odn.w;
                            }
                        }
                    }
                    //------------------------------------- Adding delta part -------------------------------------//
                }

                a1[offset+0]   = 0.0;
                c1[offset+N-2] = 0.0;

                d1[offset+0]   -= m_aa_htht__hxhx * p15[m][0];
                d1[offset+N-2] -= m_aa_htht__hxhx * p15[m][N];

                offset += N-1;
            }

            LinearEquation::func1(a1, b1, c1, d1, w1.data(), x1, rows1.size()*(N-1));

            offset = 0;
            for (unsigned int row=0; row<rows1.size(); row++)
            {
                unsigned int m=rows1.at(row);
                for (unsigned int n=1; n<=N-1; n++)
                {
                    p15[m][n] = x1[offset+(n-1)];
                }
                offset += N-1;
            }

            w1.clear();
            free(x1);
            free(d1);
            free(c1);
            free(b1);
            free(a1);
        }

        /**************************************************** x direction apprx ***************************************************/

        /**************************************************** y direction apprx ***************************************************/

        if (cols0.size() != 0)
        {
            for (unsigned int col=0; col<cols0.size(); col++)
            {
                unsigned int n = cols0.at(col);
                sn.i = n; sn.x = n*hx;
                for (unsigned int m=1; m<=M-1; m++)
                {
                    sn.j = m; sn.y = m*hy;

                    dy[m-1] = 0.0;
                    if (n>0 && n<N) dy[m-1] = p_aa_htht__hxhx*(p15[m][n-1] - 2.0*p15[m][n]   + p15[m][n+1]);
                    else if (n==0)  dy[m-1] = p_aa_htht__hxhx*(p15[m][0]   - 2.0*p15[m][1]   + p15[m][2]);
                    else if (n==N)  dy[m-1] = p_aa_htht__hxhx*(p15[m][N-2] - 2.0*p15[m][N-1] + p15[m][N]);

                    dy[m-1] += (p10[m][n]-p00[m][n]) + 2.0*p15[m][n];
                    dy[m-1] += lambda_ht*(4.0*p15[m][n]-p10[m][n]);

                    if (L <= l && l <= LLD) dy[m-1] -= 2.0*htht*(u[2*(l-L)][m][n] - V0[m][n]);
                }

                dy[0]   -= m_aa_htht__hyhy * p20[0][n];
                dy[M-2] -= m_aa_htht__hyhy * p20[M][n];

                tomasAlgorithm(ay, by, cy, dy, ry, M-1);
                for (unsigned int m=1; m<=M-1; m++) p20[m][n] = ry[m-1];
            }
        }

        if (cols1.size() != 0 && cols2.size() == 0)
        {
            //throw std::exception();

            for (unsigned int col=0; col<cols1.size(); col++)
            {
                unsigned int n = cols1.at(col);
                sn.i = n; sn.x = n*hx;
                for (unsigned int m=1; m<=M-1; m++)
                {
                    sn.j = m; sn.y = m*hy;

                    dy[m-1] = 0.0;
                    if (n>0 && n<N) dy[m-1] = p_aa_htht__hxhx*(p15[m][n-1] - 2.0*p15[m][n]   + p15[m][n+1]);
                    else if (n==0)  dy[m-1] = p_aa_htht__hxhx*(p15[m][0]   - 2.0*p15[m][1]   + p15[m][2]);
                    else if (n==N)  dy[m-1] = p_aa_htht__hxhx*(p15[m][N-2] - 2.0*p15[m][N-1] + p15[m][N]);

                    dy[m-1] += (p10[m][n]-p00[m][n]) + 2.0*p15[m][n];
                    dy[m-1] += lambda_ht*(4.0*p15[m][n]-p10[m][n]);

                    if (L <= l && l <= LLD) dy[m-1] -= 2.0*htht*(u[2*(l-L)][m][n] - V0[m][n]);

                    //------------------------------------- Adding delta part -------------------------------------//
                    for (unsigned int onj=0; onj<obsDeltaNodes.size(); onj++)
                    {
                        const ExtendedSpacePointNode &odn = obsDeltaNodes.at(onj);
                        if (odn.i == sn.i && odn.j == sn.j)
                        {
                            for (unsigned int cni=0; cni<cntPointNodes.size(); cni++)
                            {
                                const ExtendedSpacePointNode &cpn = cntPointNodes.at(cni);
                                dy[m-1] += htht * mOptParameter.k[cpn.id][odn.id] * p20[cpn.j][cpn.i] * (cpn.w * (hx*hy)) * odn.w;
                            }

                            for (unsigned int i=0; i<Nc; i++)
                            {
                                dy[m-1] += 2.0 * r * htht * mOptParameter.k[i][odn.id] * gpi(i, l, u_info, mOptParameter)*sgn(g0i(i, l, u_info, mOptParameter)) * odn.w;
                            }
                        }
                    }

                    //------------------------------------- Adding delta part -------------------------------------//
                }

                dy[0]   -= m_aa_htht__hyhy * p20[0][n];
                dy[M-2] -= m_aa_htht__hyhy * p20[M][n];

                tomasAlgorithm(ay, by, cy, dy, ry, M-1);
                for (unsigned int m=1; m<=M-1; m++) p20[m][n] = ry[m-1];
            }
        }

        if (cols1.size() != 0 && cols2.size() != 0)
        {
            //throw std::exception();

            unsigned int cols1_size = cols1.size()*(M-1);
            double* a2 = (double*) malloc(sizeof(double)*cols1_size);
            double* b2 = (double*) malloc(sizeof(double)*cols1_size);
            double* c2 = (double*) malloc(sizeof(double)*cols1_size);
            double* d2 = (double*) malloc(sizeof(double)*cols1_size);
            double* x2 = (double*) malloc(sizeof(double)*cols1_size);
            DoubleMatrix w2(cols1_size, cols1_size, 0.0);

            unsigned int offset = 0;
            for (unsigned int col=0; col<cols1.size(); col++)
            {
                unsigned int n = cols1.at(col);
                sn.i = n; sn.x = n*hx;

                for (unsigned int m=1; m<=M-1; m++)
                {
                    sn.j = m; sn.y = m*hy;

                    d2[offset+(m-1)] = 0.0;
                    if (n>0 && n<N) d2[offset+(m-1)] = p_aa_htht__hxhx*(p15[m][n-1] - 2.0*p15[m][n]   + p15[m][n+1]);
                    else if (n==0)  d2[offset+(m-1)] = p_aa_htht__hxhx*(p15[m][0]   - 2.0*p15[m][1]   + p15[m][2]);
                    else if (n==N)  d2[offset+(m-1)] = p_aa_htht__hxhx*(p15[m][N-2] - 2.0*p15[m][N-1] + p15[m][N]);

                    d2[offset+(m-1)] += (p10[m][n]-p00[m][n]) + 2.0*p15[m][n];
                    d2[offset+(m-1)] += lambda*ht*(4.0*p15[m][n]-p10[m][n]);

                    if (L <= l && l <= LLD) d2[offset+(m-1)] -= 2.0*htht*(u[2*(l-L)][m][n] - V0[m][n]);

                    a2[offset+(m-1)] = m_aa_htht__hyhy;
                    b2[offset+(m-1)] = p_aa_htht__hyhy___lambda_ht;
                    c2[offset+(m-1)] = m_aa_htht__hyhy;

                    //------------------------------------- Adding delta part -------------------------------------//
                    for (unsigned int onj=0; onj<obsDeltaNodes.size(); onj++)
                    {
                        const ExtendedSpacePointNode &odn = obsDeltaNodes.at(onj);
                        if (odn.i == sn.i && odn.j == sn.j)
                        {
                            for (unsigned int cni=0; cni<cntPointNodes.size(); cni++)
                            {
                                const ExtendedSpacePointNode &cpn = cntPointNodes[cni];

                                bool found = false;
                                for (unsigned int cs=0; cs<cols1.size(); cs++)
                                {
                                    if (cpn.i == cols1[cs])
                                    {
                                        found = true;
                                        w2[offset+(m-1)][cs*(M-1)+(cpn.j-1)] -= htht * mOptParameter.k[cpn.id][odn.id] * (cpn.w * (hx*hy)) * odn.w;
                                        break;
                                    }
                                }

                                if (!found)
                                {
                                    d2[offset+(m-1)] += htht * mOptParameter.k[cpn.id][odn.id] * p20[cpn.j][cpn.i] * (cpn.w * (hx*hy)) * odn.w;
                                }
                            }

                            for (unsigned int i=0; i<Nc; i++)
                            {
                                d2[offset+(m-1)] += 2.0 * r * htht *  mOptParameter.k[i][odn.id] * gpi(i, l, u_info, mOptParameter)*sgn(g0i(i, l, u_info, mOptParameter)) * odn.w;
                            }
                        }
                    }
                    //------------------------------------- Adding delta part -------------------------------------//
                }

                a2[offset+0]   = 0.0;
                c2[offset+M-2] = 0.0;

                d2[offset+0]   -= m_aa_htht__hyhy * p20[0][n];
                d2[offset+M-2] -= m_aa_htht__hyhy * p20[M][n];

                offset += M-1;
            }

            LinearEquation::func1(a2, b2, c2, d2, w2.data(), x2, cols1.size()*(M-1));

            offset = 0;
            for (unsigned int col=0; col<cols1.size(); col++)
            {
                unsigned int n=cols1.at(col);
                for (unsigned int m=1; m<=M-1; m++)
                {
                    p20[m][n] = x2[offset+(m-1)];
                }
                offset += M-1;
            }

            w2.clear();
            free(x2);
            free(d2);
            free(c2);
            free(b2);
            free(a2);
        }

        /**************************************************** y direction apprx ***************************************************/

        for (unsigned int m=0; m<=M; m++)
        {
            for (unsigned int n=0; n<=N; n++)
            {
                p00[m][n] = p10[m][n];
                p10[m][n] = p20[m][n];
                p05[m][n] = p15[m][n];
            }
        }

        if (use == true) b_add2Info(p20, p_info, cntPointNodes, l, hx, hy);

        b_layerInfo(p20, l);
    }

    free(rx);
    free(dx);
    free(cx);
    free(bx);
    free(ax);

    free(ry);
    free(dy);
    free(cy);
    free(by);
    free(ay);

    rows0.clear();
    rows1.clear();
    rows2.clear();

    cols0.clear();
    cols1.clear();
    cols2.clear();

    obsDeltaNodes.clear();
    cntPointNodes.clear();

    p00.clear();
    p05.clear();
    p10.clear();
    p15.clear();
    p20.clear();
    //puts("+void Problem2HNDirichlet::solveBackwardIBVP(const std::vector<DoubleMatrix> &u, spif_vector &p_info, bool use, const spif_vector &u_info) const");
}

void Problem2HNDirichlet::solveBackwardIBVP1(const std::vector<DoubleMatrix> &u, spif_vector &p_info, bool use, const spif_vector &u_info) const
{
    //printf("eta: 0: %f %f 1: %f %f\n", mOptParameter.eta[0].x, mOptParameter.eta[0].y, mOptParameter.eta[1].x, mOptParameter.eta[1].y);
    //printf("xi : 0: %f %f 1: %f %f\n", mOptParameter.xi[0].x, mOptParameter.xi[0].y, mOptParameter.xi[1].x, mOptParameter.xi[1].y);
    const Dimension dimX = spaceDimension(Dimension::DimensionX);
    const Dimension dimY = spaceDimension(Dimension::DimensionY);
    const Dimension time = timeDimension();

    const unsigned int N = dimX.sizeN();
    const unsigned int M = dimY.sizeN();
    const unsigned int L = time.sizeN();
    const unsigned int LLD = L+LD;

    const double hx = dimX.step();
    const double hy = dimY.step();
    const double ht = time.step();

    const double a = mEquParameter.a;
    const double lambda = mEquParameter.lambda;
    const unsigned int No = mEquParameter.No;
    const unsigned int Nc = mEquParameter.Nc;

    const double m_aa_htht__hxhx = -(a*a*ht*ht)/(hx*hx);
    const double p_aa_htht__hxhx___lambda_ht = +2.0 + 2.0*(a*a*ht*ht)/(hx*hx) + (lambda*ht);
    const double p_aa_htht__hxhx = +(a*a*ht*ht)/(hx*hx);

    const double m_aa_htht__hyhy = -(a*a*ht*ht)/(hy*hy);
    const double p_aa_htht__hyhy___lambda_ht = +2.0 + 2.0*(a*a*ht*ht)/(hy*hy) + (lambda*ht);
    const double p_aa_htht__hyhy = +(a*a*ht*ht)/(hy*hy);

    const double htht = ht*ht;
    const double lambda_ht = lambda*ht;

    const double aa__hxhx = (a*a)/(hx*hx);
    const double aa__hyhy = (a*a)/(hy*hy);

    DoubleMatrix p00(M+1, N+1);
    DoubleMatrix p10(M+1, N+1);
    DoubleMatrix p15(M+1, N+1);
    DoubleMatrix p20(M+1, N+1);

    //--------------------------------------------------------------------------------------------//
    espn_vector obsDeltaNodes, cntPointNodes;
    for (unsigned int j=0; j<No; j++) distributeDelta0(mOptParameter.xi[j], j, obsDeltaNodes, dimX, dimY, 4, 4);
    for (unsigned int i=0; i<Nc; i++) distributeDelta0(mOptParameter.eta[i], i, cntPointNodes, dimX, dimY, 4, 4);

    //----------------------------------------------------------------------------------------------//
    uint_vector rows0, rows1, rows2, cols0, cols1, cols2;
    b_findRowsCols(rows0, rows1, rows2, cols0, cols1, cols2, cntPointNodes, obsDeltaNodes, N, M);

    //-------------------------------------------- info --------------------------------------------//
    if (use == true) b_prepareInfo(Nc, mOptParameter.eta, p_info, LLD, dimX, dimY);
    //-------------------------------------------- info --------------------------------------------//

    //------------------------------------- initial conditions -------------------------------------//
    b_initialLayers1(p00, p10, p_info, use, cntPointNodes, obsDeltaNodes, N, M, hx, hy, ht, aa__hxhx, aa__hyhy, lambda);
    //------------------------------------- initial conditions -------------------------------------//

    double *ax = (double *) malloc(sizeof(double)*(N-1)); for (unsigned int n=1; n<=N-1; n++) ax[n-1] = m_aa_htht__hxhx;
    double *bx = (double *) malloc(sizeof(double)*(N-1)); for (unsigned int n=1; n<=N-1; n++) bx[n-1] = p_aa_htht__hxhx___lambda_ht;
    double *cx = (double *) malloc(sizeof(double)*(N-1)); for (unsigned int n=1; n<=N-1; n++) cx[n-1] = m_aa_htht__hxhx;
    double *dx = (double *) malloc(sizeof(double)*(N-1));
    double *rx = (double *) malloc(sizeof(double)*(N-1));
    ax[0] = cx[N-2] = 0.0;

    double *ay = (double *) malloc(sizeof(double)*(M-1)); for (unsigned int m=1; m<=M-1; m++) ay[m-1] = m_aa_htht__hyhy;
    double *by = (double *) malloc(sizeof(double)*(M-1)); for (unsigned int m=1; m<=M-1; m++) by[m-1] = p_aa_htht__hyhy___lambda_ht;
    double *cy = (double *) malloc(sizeof(double)*(M-1)); for (unsigned int m=1; m<=M-1; m++) cy[m-1] = m_aa_htht__hyhy;
    double *dy = (double *) malloc(sizeof(double)*(M-1));
    double *ry = (double *) malloc(sizeof(double)*(M-1));
    ay[0] = cy[M-2] = 0.0;

    SpaceNodePDE sn;

    for (unsigned int l1=2,l=LLD-l1; l1<=LLD; l1++,l=LLD-l1)
    {
        TimeNodePDE tn20; tn20.i = l; tn20.t = l*ht;
        TimeNodePDE tn15; tn15.i = l; tn15.t = l*ht+0.5*ht;

        /**************************************************** border conditions ***************************************************/

        SpaceNodePDE sn0;
        SpaceNodePDE sn1;

        sn0.i = 0; sn0.x = 0.0;
        sn1.i = N; sn1.x = hx*N;
        for (unsigned int m=0; m<=M; m++)
        {
            sn0.j = m; sn0.y = m*hy; p15[m][0] = b_boundary(sn0, tn15); p20[m][0] = b_boundary(sn0, tn20);
            sn1.j = m; sn1.y = m*hy; p15[m][N] = b_boundary(sn1, tn15); p20[m][N] = b_boundary(sn1, tn20);
        }

        sn0.j = 0; sn0.y = 0.0;
        sn1.j = M; sn1.y = hy*M;
        for (unsigned int n=0; n<=N; n++)
        {
            sn0.i = n; sn0.x = n*hx; p15[0][n] = b_boundary(sn0, tn15); p20[0][n] = b_boundary(sn0, tn20);
            sn1.i = n; sn1.x = n*hx; p15[M][n] = b_boundary(sn1, tn15); p20[M][n] = b_boundary(sn1, tn20);
        }

        /**************************************************** border conditions ***************************************************/

        /**************************************************** x direction apprx ***************************************************/

        if (rows0.size() != 0)
        {
            for (unsigned int row=0; row<rows0.size(); row++)
            {
                unsigned int m = rows0.at(row);
                sn.j = m; sn.y = m*hy;
                for (unsigned int n=1; n<=N-1; n++)
                {
                    sn.i = n; sn.x = n*hx;

                    dx[n-1] = 0.0;
                    if (m>0 && m<M)  dx[n-1] = p_aa_htht__hyhy*(p10[m-1][n] - 2.0*p10[m][n]   + p10[m+1][n]);
                    else if (m == 0) dx[n-1] = p_aa_htht__hyhy*(p10[0][n]   - 2.0*p10[1][n]   + p10[2][n]);
                    else if (m == M) dx[n-1] = p_aa_htht__hyhy*(p10[M-2][n] - 2.0*p10[M-1][n] + p10[M][n]);

                    dx[n-1] += lambda_ht*(p10[m][n] - 0.5*(p10[m][n]-p00[m][n]));
                    dx[n-1] += 2.0*p10[m][n] + (p10[m][n]-p00[m][n]);

                    if (L <= l && l <= LLD) dx[n-1] -= 2.0*htht*(u[2*(l-L)+1][m][n] - V0[m][n]);
                }

                dx[0]   -= m_aa_htht__hxhx * p15[m][0];
                dx[N-2] -= m_aa_htht__hxhx * p15[m][N];

                tomasAlgorithm(ax, bx, cx, dx, rx, N-1);
                for (unsigned int n=1; n<=N-1; n++) p15[m][n] = rx[n-1];
            }
        }

        if (rows1.size() != 0 && rows2.size() == 0)
        {
            //throw std::exception();

            for (unsigned int row=0; row<rows1.size(); row++)
            {
                unsigned int m = rows1.at(row);
                sn.j = m; sn.y = m*hy;
                for (unsigned int n=1; n<=N-1; n++)
                {
                    sn.i = n; sn.x = n*hx;

                    dx[n-1] = 0.0;
                    if (m>0 && m<M)  dx[n-1] = p_aa_htht__hyhy*(p10[m-1][n] - 2.0*p10[m][n]   + p10[m+1][n]);
                    else if (m == 0) dx[n-1] = p_aa_htht__hyhy*(p10[0][n]   - 2.0*p10[1][n]   + p10[2][n]);
                    else if (m == M) dx[n-1] = p_aa_htht__hyhy*(p10[M-2][n] - 2.0*p10[M-1][n] + p10[M][n]);

                    dx[n-1] += lambda_ht*(p10[m][n] - 0.5*(p10[m][n]-p00[m][n]));
                    dx[n-1] += 2.0*p10[m][n] + (p10[m][n]-p00[m][n]);

                    if (L <= l && l <= LLD) dx[n-1] -= 2.0*htht*(u[2*(l-L)+1][m][n] - V0[m][n]);

                    //------------------------------------- Adding delta part -------------------------------------//
                    for (unsigned int onj=0; onj<obsDeltaNodes.size(); onj++)
                    {
                        const ExtendedSpacePointNode &odn = obsDeltaNodes.at(onj);
                        if (odn.i == sn.i && odn.j == sn.j)
                        {
                            for (unsigned int cni=0; cni<cntPointNodes.size(); cni++)
                            {
                                const ExtendedSpacePointNode &cpn = cntPointNodes[cni];
                                dx[n-1] += htht * mOptParameter.k[cpn.id][odn.id] * p15[cpn.j][cpn.i] * (cpn.w * (hx*hy)) * odn.w;
                            }

                            for (unsigned int i=0; i<Nc; i++)
                            {
                                dx[n-1] += 2.0*r * htht * mOptParameter.k[i][odn.id] * gpi(i, l, u_info, mOptParameter)*sgn(g0i(i, l, u_info, mOptParameter)) * odn.w;
                            }
                        }
                    }
                    //------------------------------------- Adding delta part -------------------------------------//
                }

                dx[0]   -= m_aa_htht__hxhx * p15[m][0];
                dx[N-2] -= m_aa_htht__hxhx * p15[m][N];

                tomasAlgorithm(ax, bx, cx, dx, rx, N-1);
                for (unsigned int n=1; n<=N-1; n++) p15[m][n] = rx[n-1];
            }
        }

        if (rows1.size() != 0 && rows2.size() != 0)
        {
            //throw std::exception();

            unsigned int row1_size = rows1.size()*(N-1);
            double* a1 = (double*) malloc(sizeof(double)*row1_size);
            double* b1 = (double*) malloc(sizeof(double)*row1_size);
            double* c1 = (double*) malloc(sizeof(double)*row1_size);
            double* d1 = (double*) malloc(sizeof(double)*row1_size);
            double* x1 = (double*) malloc(sizeof(double)*row1_size);
            DoubleMatrix w1(row1_size, row1_size, 0.0);

            unsigned int offset = 0;
            for (unsigned int row=0; row<rows1.size(); row++)
            {
                unsigned int m = rows1.at(row);
                sn.j = m; sn.y = m*hy;

                for (unsigned int n=1; n<=N-1; n++)
                {
                    sn.i = n; sn.x = n*hx;

                    d1[offset+(n-1)] = 0.0;
                    if (m>0 && m<M)  d1[offset+(n-1)] = p_aa_htht__hyhy*(p10[m-1][n] - 2.0*p10[m][n]   + p10[m+1][n]);
                    else if (m == 0) d1[offset+(n-1)] = p_aa_htht__hyhy*(p10[0][n]   - 2.0*p10[1][n]   + p10[2][n]);
                    else if (m == M) d1[offset+(n-1)] = p_aa_htht__hyhy*(p10[M-2][n] - 2.0*p10[M-1][n] + p10[M][n]);

                    dx[offset+(n-1)] += lambda_ht*(p10[m][n] - 0.5*(p10[m][n]-p00[m][n]));
                    dx[offset+(n-1)] += 2.0*p10[m][n] + (p10[m][n]-p00[m][n]);

                    if (L <= l && l <= LLD) d1[offset+(n-1)] -= 2.0*htht*(u[2*(l-L)+1][m][n] - V0[m][n]);

                    a1[offset+(n-1)] = m_aa_htht__hxhx;
                    b1[offset+(n-1)] = p_aa_htht__hxhx___lambda_ht;
                    c1[offset+(n-1)] = m_aa_htht__hxhx;

                    //------------------------------------- Adding delta part -------------------------------------//
                    for (unsigned int onj=0; onj<obsDeltaNodes.size(); onj++)
                    {
                        const ExtendedSpacePointNode &odn = obsDeltaNodes.at(onj);
                        if (odn.i == sn.i && odn.j == sn.j)
                        {
                            for (unsigned int cni=0; cni<cntPointNodes.size(); cni++)
                            {
                                const ExtendedSpacePointNode &cpn = cntPointNodes[cni];

                                bool found = false;
                                for (unsigned int rs=0; rs<rows1.size(); rs++)
                                {
                                    if (cpn.j == rows1[rs])
                                    {
                                        found = true;
                                        w1[offset+(n-1)][rs*(N-1)+(cpn.i-1)] -= htht * mOptParameter.k[cpn.id][odn.id] * (cpn.w * (hx*hy)) * odn.w;
                                        break;
                                    }
                                }

                                if (!found)
                                {
                                    d1[offset+(n-1)] += htht * mOptParameter.k[cpn.id][odn.id] * p15[cpn.j][cpn.i] * (cpn.w * (hx*hy)) * odn.w;
                                }
                            }

                            for (unsigned int i=0; i<Nc; i++)
                            {
                                d1[offset+(n-1)] += 2.0 * r * ht *  mOptParameter.k[i][odn.id] * gpi(i, l, u_info, mOptParameter)*sgn(g0i(i, l, u_info, mOptParameter)) * odn.w;
                            }
                        }
                    }
                    //------------------------------------- Adding delta part -------------------------------------//
                }

                a1[offset+0]   = 0.0;
                c1[offset+N-2] = 0.0;

                d1[offset+0]   -= m_aa_htht__hxhx * p15[m][0];
                d1[offset+N-2] -= m_aa_htht__hxhx * p15[m][N];

                offset += N-1;
            }
            printf("%d 1 1-1 ",l);
            LinearEquation::func1(a1, b1, c1, d1, w1.data(), x1, row1_size);
            printf("1-2 ");
            offset = 0;
            for (unsigned int row=0; row<rows1.size(); row++)
            {
                unsigned int m=rows1.at(row);
                for (unsigned int n=1; n<=N-1; n++)
                {
                    p15[m][n] = x1[offset+(n-1)];
                }
                offset += N-1;
            }

            w1.clear(); printf("1-3 ");
            free(x1);printf("1-4 ");
            free(d1);printf("1-5 ");
            free(c1);printf("1-6 ");
            free(b1);printf("1-7 ");
            free(a1);printf("1-8 ");
        }

        /**************************************************** x direction apprx ***************************************************/

        /**************************************************** y direction apprx ***************************************************/

        if (cols0.size() != 0)
        {
            for (unsigned int col=0; col<cols0.size(); col++)
            {
                unsigned int n = cols0.at(col);
                sn.i = n; sn.x = n*hx;
                for (unsigned int m=1; m<=M-1; m++)
                {
                    sn.j = m; sn.y = m*hy;

                    dy[m-1] = 0.0;
                    if (n>0 && n<N) dy[m-1] = p_aa_htht__hxhx*(p15[m][n-1] - 2.0*p15[m][n]   + p15[m][n+1]);
                    else if (n==0)  dy[m-1] = p_aa_htht__hxhx*(p15[m][0]   - 2.0*p15[m][1]   + p15[m][2]);
                    else if (n==N)  dy[m-1] = p_aa_htht__hxhx*(p15[m][N-2] - 2.0*p15[m][N-1] + p15[m][N]);

                    dy[m-1] += lambda_ht*(p15[m][n] - 0.5*(p10[m][n]-p00[m][n]));
                    dy[m-1] += 2.0*p15[m][n] + (p10[m][n]-p00[m][n]);

                    if (L <= l && l <= LLD) dy[m-1] -= 2.0*htht*(u[2*(l-L)][m][n] - V0[m][n]);
                }

                dy[0]   -= m_aa_htht__hyhy * p20[0][n];
                dy[M-2] -= m_aa_htht__hyhy * p20[M][n];

                tomasAlgorithm(ay, by, cy, dy, ry, M-1);
                for (unsigned int m=1; m<=M-1; m++) p20[m][n] = ry[m-1];
            }
        }

        if (cols1.size() != 0 && cols2.size() == 0)
        {
            //throw std::exception();

            for (unsigned int col=0; col<cols1.size(); col++)
            {
                unsigned int n = cols1.at(col);
                sn.i = n; sn.x = n*hx;
                for (unsigned int m=1; m<=M-1; m++)
                {
                    sn.j = m; sn.y = m*hy;

                    dy[m-1] = 0.0;
                    if (n>0 && n<N) dy[m-1] = p_aa_htht__hxhx*(p15[m][n-1] - 2.0*p15[m][n]   + p15[m][n+1]);
                    else if (n==0)  dy[m-1] = p_aa_htht__hxhx*(p15[m][0]   - 2.0*p15[m][1]   + p15[m][2]);
                    else if (n==N)  dy[m-1] = p_aa_htht__hxhx*(p15[m][N-2] - 2.0*p15[m][N-1] + p15[m][N]);

                    dy[m-1] += lambda_ht*(p15[m][n] - 0.5*(p10[m][n]-p00[m][n]));
                    dy[m-1] += 2.0*p15[m][n] + (p10[m][n]-p00[m][n]);

                    if (L <= l && l <= LLD) dy[m-1] -= 2.0*htht*(u[2*(l-L)][m][n] - V0[m][n]);

                    //------------------------------------- Adding delta part -------------------------------------//
                    for (unsigned int onj=0; onj<obsDeltaNodes.size(); onj++)
                    {
                        const ExtendedSpacePointNode &odn = obsDeltaNodes.at(onj);
                        if (odn.i == sn.i && odn.j == sn.j)
                        {
                            for (unsigned int cni=0; cni<cntPointNodes.size(); cni++)
                            {
                                const ExtendedSpacePointNode &cpn = cntPointNodes.at(cni);
                                dy[m-1] += htht * mOptParameter.k[cpn.id][odn.id] * p20[cpn.j][cpn.i] * (cpn.w * (hx*hy)) * odn.w;
                            }

                            for (unsigned int i=0; i<Nc; i++)
                            {
                                dy[m-1] += 2.0 * r * htht * mOptParameter.k[i][odn.id] * gpi(i, l, u_info, mOptParameter)*sgn(g0i(i, l, u_info, mOptParameter)) * odn.w;
                            }
                        }
                    }

                    //------------------------------------- Adding delta part -------------------------------------//
                }

                dy[0]   -= m_aa_htht__hyhy * p20[0][n];
                dy[M-2] -= m_aa_htht__hyhy * p20[M][n];

                tomasAlgorithm(ay, by, cy, dy, ry, M-1);
                for (unsigned int m=1; m<=M-1; m++) p20[m][n] = ry[m-1];
            }
        }

        if (cols1.size() != 0 && cols2.size() != 0)
        {
            //throw std::exception();

            unsigned int cols1_size = cols1.size()*(M-1);
            double* a2 = (double*) malloc(sizeof(double)*cols1_size);
            double* b2 = (double*) malloc(sizeof(double)*cols1_size);
            double* c2 = (double*) malloc(sizeof(double)*cols1_size);
            double* d2 = (double*) malloc(sizeof(double)*cols1_size);
            double* x2 = (double*) malloc(sizeof(double)*cols1_size);
            DoubleMatrix w2(cols1_size, cols1_size, 0.0);

            unsigned int offset = 0;
            for (unsigned int col=0; col<cols1.size(); col++)
            {
                unsigned int n = cols1.at(col);
                sn.i = n; sn.x = n*hx;

                for (unsigned int m=1; m<=M-1; m++)
                {
                    sn.j = m; sn.y = m*hy;

                    d2[offset+(m-1)] = 0.0;
                    if (n>0 && n<N) d2[offset+(m-1)] = p_aa_htht__hxhx*(p15[m][n-1] - 2.0*p15[m][n]   + p15[m][n+1]);
                    else if (n==0)  d2[offset+(m-1)] = p_aa_htht__hxhx*(p15[m][0]   - 2.0*p15[m][1]   + p15[m][2]);
                    else if (n==N)  d2[offset+(m-1)] = p_aa_htht__hxhx*(p15[m][N-2] - 2.0*p15[m][N-1] + p15[m][N]);

                    dy[offset+(m-1)] += lambda_ht*(p15[m][n] - 0.5*(p10[m][n]-p00[m][n]));
                    dy[offset+(m-1)] += 2.0*p15[m][n] + (p10[m][n]-p00[m][n]);

                    if (L <= l && l <= LLD) d2[offset+(m-1)] -= 2.0*htht*(u[2*(l-L)][m][n] - V0[m][n]);

                    a2[offset+(m-1)] = m_aa_htht__hyhy;
                    b2[offset+(m-1)] = p_aa_htht__hyhy___lambda_ht;
                    c2[offset+(m-1)] = m_aa_htht__hyhy;

                    //------------------------------------- Adding delta part -------------------------------------//
                    for (unsigned int onj=0; onj<obsDeltaNodes.size(); onj++)
                    {
                        const ExtendedSpacePointNode &odn = obsDeltaNodes.at(onj);
                        if (odn.i == sn.i && odn.j == sn.j)
                        {
                            for (unsigned int cni=0; cni<cntPointNodes.size(); cni++)
                            {
                                const ExtendedSpacePointNode &cpn = cntPointNodes[cni];

                                bool found = false;
                                for (unsigned int cs=0; cs<cols1.size(); cs++)
                                {
                                    if (cpn.i == cols1[cs])
                                    {
                                        found = true;
                                        w2[offset+(m-1)][cs*(M-1)+(cpn.j-1)] -= htht * mOptParameter.k[cpn.id][odn.id] * (cpn.w * (hx*hy)) * odn.w;
                                        break;
                                    }
                                }

                                if (!found)
                                {
                                    d2[offset+(m-1)] += htht * mOptParameter.k[cpn.id][odn.id] * p20[cpn.j][cpn.i] * (cpn.w * (hx*hy)) * odn.w;
                                }
                            }

                            for (unsigned int i=0; i<Nc; i++)
                            {
                                d2[offset+(m-1)] += 2.0 * r * htht *  mOptParameter.k[i][odn.id] * gpi(i, l, u_info, mOptParameter)*sgn(g0i(i, l, u_info, mOptParameter)) * odn.w;
                            }
                        }
                    }
                    //------------------------------------- Adding delta part -------------------------------------//
                }

                a2[offset+0]   = 0.0;
                c2[offset+M-2] = 0.0;

                d2[offset+0]   -= m_aa_htht__hyhy * p20[0][n];
                d2[offset+M-2] -= m_aa_htht__hyhy * p20[M][n];

                offset += M-1;
            }

            LinearEquation::func1(a2, b2, c2, d2, w2.data(), x2, cols1.size()*(M-1));

            offset = 0;
            for (unsigned int col=0; col<cols1.size(); col++)
            {
                unsigned int n=cols1.at(col);
                for (unsigned int m=1; m<=M-1; m++)
                {
                    p20[m][n] = x2[offset+(m-1)];
                }
                offset += M-1;
            }

            w2.clear();
            free(x2);
            free(d2);
            free(c2);
            free(b2);
            free(a2);
        }

        /**************************************************** y direction apprx ***************************************************/

        for (unsigned int m=0; m<=M; m++)
        {
            for (unsigned int n=0; n<=N; n++)
            {
                p00[m][n] = p10[m][n];
                p10[m][n] = p20[m][n];
            }
        }

        if (use == true) b_add2Info(p20, p_info, cntPointNodes, l, hx, hy);

        b_layerInfo(p20, l);

    }

    free(rx);
    free(dx);
    free(cx);
    free(bx);
    free(ax);

    free(ry);
    free(dy);
    free(cy);
    free(by);
    free(ay);

    rows0.clear();
    rows1.clear();
    rows2.clear();

    cols0.clear();
    cols1.clear();
    cols2.clear();

    obsDeltaNodes.clear();
    cntPointNodes.clear();

    p00.clear();
    p10.clear();
    p15.clear();
    p20.clear();
}

double Problem2HNDirichlet::b_initial1(const SpaceNodePDE &sn UNUSED_PARAM) const
{
    return 0.0;
}

double Problem2HNDirichlet::b_initial2(const SpaceNodePDE &sn UNUSED_PARAM) const
{
    return 0.0;
}

double Problem2HNDirichlet::b_boundary(const SpaceNodePDE &, const TimeNodePDE &, BoundaryType) const
{
    return 0.0;
}

void Problem2HNDirichlet::b_findRowsCols(uint_vector &rows0, uint_vector &rows1, uint_vector &rows2, uint_vector &cols0, uint_vector &cols1, uint_vector &cols2,
                                         espn_vector &cntPointNodes, espn_vector &obsDeltaNodes, unsigned int N, unsigned int M) const
{
    for (unsigned int m=1; m<=M-1; m++)
    {
        bool found1 = false;
        bool found2 = false;
        for (unsigned int j=0; j<obsDeltaNodes.size(); j++)
        {
            const ExtendedSpacePointNode &odn = obsDeltaNodes.at(j);
            if (odn.j == m)
            {
                found1 = true;
                for (unsigned int i=0; i<cntPointNodes.size(); i++)
                {
                    const ExtendedSpacePointNode &cpn = cntPointNodes.at(i);
                    if (cpn.j == m)
                    {
                        found2 = true;
                        break;
                    }
                }
                break;
            }
        }
        if (found1 == false && found2 == false) if(std::find(rows0.begin(), rows0.end(), m) == rows0.end()) rows0.push_back(m);
        if (found1 == true  && found2 == true)  if(std::find(rows2.begin(), rows2.end(), m) == rows2.end()) rows2.push_back(m);
        if (found1 == true)                     if(std::find(rows1.begin(), rows1.end(), m) == rows1.end()) rows1.push_back(m);
    }

    for (unsigned int n=1; n<=N-1; n++)
    {
        bool found1 = false;
        bool found2 = false;
        for (unsigned int j=0; j<obsDeltaNodes.size(); j++)
        {
            const ExtendedSpacePointNode &odn = obsDeltaNodes.at(j);
            if (odn.i == n)
            {
                found1 = true;
                for (unsigned int i=0; i<cntPointNodes.size(); i++)
                {
                    const ExtendedSpacePointNode &cpn = cntPointNodes.at(i);
                    if (cpn.i == n)
                    {
                        found2 = true;
                        break;
                    }
                }
                break;
            }
        }
        if (found1 == false && found2 == false) if(std::find(cols0.begin(), cols0.end(), n) == cols0.end()) cols0.push_back(n);
        if (found1 == true  && found2 == true)  if(std::find(cols2.begin(), cols2.end(), n) == cols2.end()) cols2.push_back(n);
        if (found1 == true)                     if(std::find(cols1.begin(), cols1.end(), n) == cols1.end()) cols1.push_back(n);
    }
}

void Problem2HNDirichlet::b_initialLayers(DoubleMatrix &p00, DoubleMatrix &p05, DoubleMatrix &p10, spif_vector &p_info, bool use, espn_vector &cntPointNodes, espn_vector &obsDeltaNodes,
                                          unsigned int N, unsigned int M, double hx, double hy, double ht, double aa__hxhx, double aa__hyhy, double lambda) const
{
    unsigned int L = mtimeDimension.sizeN();
    unsigned int LLD = L+LD;

    SpaceNodePDE sn;
    for (unsigned int m=0; m<=M; m++)
    {
        sn.j = m; sn.y = m*hy;
        for (unsigned int n=0; n<=N; n++)
        {
            sn.i = n; sn.x = n*hx;
            p00[m][n] = b_initial1(sn);
        }
    }

    TimeNodePDE tn05; tn05.i = LLD-1; tn05.t = LLD*ht + 0.5*ht;
    TimeNodePDE tn10; tn10.i = LLD-1; tn10.t = LLD*ht;

    SpaceNodePDE sn0, sn1;
    sn0.i = 0; sn0.x = 0.0;
    sn1.i = N; sn1.x = N*hx;
    for (unsigned int m=0; m<=M; m++)
    {
        sn0.j = m; sn0.y = m*hy; p05[m][0] = f_boundary(sn0, tn05); p10[m][0] = f_boundary(sn0, tn10);
        sn1.j = m; sn1.y = m*hy; p05[m][N] = f_boundary(sn1, tn05); p10[m][N] = f_boundary(sn1, tn10);
    }

    sn0.j = 0;  sn0.y = 0.0;
    sn1.j = M; sn1.y = M*hy;
    for (unsigned int n=0; n<=N; n++)
    {
        sn0.i = n; sn0.x = n*hx; p05[0][n] = f_boundary(sn0, tn05); p10[0][n] = f_boundary(sn0, tn10);
        sn1.i = n; sn1.x = n*hx; p05[M][n] = f_boundary(sn1, tn05); p10[M][n] = f_boundary(sn1, tn10);
    }

    for (unsigned int m=1; m<=M-1; m++)
    {
        sn.j = m; sn.y = m*hy;
        for (unsigned int n=1; n<=N-1; n++)
        {
            sn.i = n; sn.x = n*hx;

            double sum = 0.0;
            sum += aa__hxhx*(p00[m][n-1]-2.0*p00[m][n]+p00[m][n+1]);
            sum += aa__hyhy*(p00[m-1][n]-2.0*p00[m][n]+p00[m+1][n]);
            sum += lambda*b_initial2(sn);

            for (unsigned int odn=0; odn<obsDeltaNodes.size(); odn++)
            {
                const ExtendedSpacePointNode &obsNode = obsDeltaNodes.at(odn);
                if (obsNode.i == n && obsNode.j == m)
                {
                    for (unsigned int cpn=0; cpn<cntPointNodes.size(); cpn++)
                    {
                        const ExtendedSpacePointNode &cntNode = cntPointNodes.at(cpn);
                        if (cntNode.i == n && cntNode.j == m)
                        {
                            sum += mOptParameter.k[cntNode.id][obsNode.id] * p00[m][n]*(cntNode.w*(hx*hy)) * obsNode.w;
                        }
                    }
                }
            }

            p05[m][n] = p00[m][n] - b_initial2(sn)*ht*0.5 + 0.125*ht*ht*sum;
            p10[m][n] = p00[m][n] - b_initial2(sn)*ht     + 0.500*ht*ht*sum;
        }
    }

    if (use == true) b_add2Info(p00, p_info, cntPointNodes, LLD, hx, hy);
    b_layerInfo(p00, LLD);
    if (use == true) b_add2Info(p10, p_info, cntPointNodes, LLD-1, hx, hy);
    b_layerInfo(p10, LLD-1);
}

void Problem2HNDirichlet::b_initialLayers1(DoubleMatrix &p00, DoubleMatrix &p10, spif_vector &p_info, bool use, espn_vector &cntPointNodes, espn_vector &obsDeltaNodes,
                                           unsigned int N, unsigned int M, double hx, double hy, double ht, double aa__hxhx, double aa__hyhy, double lambda) const
{
    unsigned int L = mtimeDimension.sizeN();
    unsigned int LLD = L+LD;

    /************************************************************************/
    SpaceNodePDE sn;
    for (unsigned int m=0; m<=M; m++)
    {
        sn.j = m; sn.y = m*hy;
        for (unsigned int n=0; n<=N; n++)
        {
            sn.i = n; sn.x = n*hx;
            p00[m][n] = b_initial1(sn);
        }
    }

    /************************************************************************/
    TimeNodePDE tn10; tn10.i = LLD-1; tn10.t = LLD*ht;

    SpaceNodePDE sn0, sn1;
    sn0.i = 0; sn0.x = 0.0;
    sn1.i = N; sn1.x = N*hx;
    for (unsigned int m=0; m<=M; m++)
    {
        sn0.j = m; sn0.y = m*hy; p10[m][0] = f_boundary(sn0, tn10);
        sn1.j = m; sn1.y = m*hy; p10[m][N] = f_boundary(sn1, tn10);
    }

    sn0.j = 0;  sn0.y = 0.0;
    sn1.j = M; sn1.y = M*hy;
    for (unsigned int n=0; n<=N; n++)
    {
        sn0.i = n; sn0.x = n*hx; p10[0][n] = f_boundary(sn0, tn10);
        sn1.i = n; sn1.x = n*hx; p10[M][n] = f_boundary(sn1, tn10);
    }

    /************************************************************************/

    for (unsigned int m=1; m<=M-1; m++)
    {
        sn.j = m; sn.y = m*hy;
        for (unsigned int n=1; n<=N-1; n++)
        {
            sn.i = n; sn.x = n*hx;

            double sum = 0.0;
            sum += aa__hxhx*(p00[m][n-1]-2.0*p00[m][n]+p00[m][n+1]);
            sum += aa__hyhy*(p00[m-1][n]-2.0*p00[m][n]+p00[m+1][n]);
            sum += lambda*b_initial2(sn);

            for (unsigned int odn=0; odn<obsDeltaNodes.size(); odn++)
            {
                const ExtendedSpacePointNode &obsNode = obsDeltaNodes.at(odn);
                if (obsNode.i == n && obsNode.j == m)
                {
                    for (unsigned int cpn=0; cpn<cntPointNodes.size(); cpn++)
                    {
                        const ExtendedSpacePointNode &cntNode = cntPointNodes.at(cpn);
                        if (cntNode.i == n && cntNode.j == m)
                        {
                            sum += mOptParameter.k[cntNode.id][obsNode.id] * p00[m][n]*(cntNode.w*(hx*hy)) * obsNode.w;
                        }
                    }
                }
            }

            p10[m][n] = p00[m][n] - b_initial2(sn)*ht + 0.5*ht*ht*sum;
        }
    }

    if (use == true) b_add2Info(p00, p_info, cntPointNodes, LLD, hx, hy);
    b_layerInfo(p00, LLD);

    if (use == true) b_add2Info(p10, p_info, cntPointNodes, LLD-1, hx, hy);
    b_layerInfo(p10, LLD-1);
}

void Problem2HNDirichlet::b_prepareInfo(unsigned int Nc, const std::vector<SpacePoint> &points, spif_vector &p_info,
                                        unsigned int LLD, const Dimension &dimX, const Dimension &dimY) const
{
    p_info.resize(Nc);
    for (unsigned int i=0; i<Nc; i++)
    {
        SpacePointInfo &inf = p_info[i];
        const SpacePoint &sp = points[i];
        inf.id = i;
        inf.x = sp.x;
        inf.y = sp.y;

        inf.i = (unsigned int)( round( inf.x * dimX.sizeN() ) );
        inf.j = (unsigned int)( round( inf.y * dimY.sizeN() ) );

        inf.layerNumber = LLD+1;
        inf.u = new double[LLD+1];
        inf.ux = new double[LLD+1];
        inf.uy = new double[LLD+1];
    }
}

void Problem2HNDirichlet::b_add2Info(const DoubleMatrix &p, spif_vector &p_info, const espn_vector &cntPointNodes, unsigned int ln, double hx, double hy, int method) const
{
    if (method == 1)
    {
        for (unsigned int j=0; j<p_info.size(); j++) p_info[j].u[ln] = p_info[j].ux[ln] = p_info[j].uy[ln] = 0.0;

        for (unsigned int cpi=0; cpi<cntPointNodes.size(); cpi++)
        {
            const ExtendedSpacePointNode &cpn = cntPointNodes[cpi];
            const SpacePointInfo &pi = p_info[cpn.id];

            pi.u[ln] += p[cpn.j][cpn.i] * (cpn.w * (hx*hy));

            pi.ux[ln] = (p[pi.j][pi.i+1] - p[pi.j][pi.i-1])/(2.0*hx);
            pi.uy[ln] = (p[pi.j+1][pi.i] - p[pi.j-1][pi.i])/(2.0*hy);

            //ui.ux[ln] = (u[ui.j][ui.i-2]-8.0*u[ui.j][ui.i-1]+8.0*u[ui.j][ui.i+1]-u[ui.j][ui.i+2])/(12.0*hx);
            //ui.uy[ln] = (u[ui.j-2][ui.i]-8.0*u[ui.j-1][ui.i]+8.0*u[ui.j+1][ui.i]-u[ui.j+2][ui.i])/(12.0*hy);
        }
    }

    if (method == 2)
    {
        for (unsigned int j=0; j<p_info.size(); j++) p_info[j].u[ln] = p_info[j].ux[ln] = p_info[j].uy[ln] = 0.0;

        for (unsigned int cpi=0; cpi<cntPointNodes.size(); cpi++)
        {
            const ExtendedSpacePointNode &cpn = cntPointNodes[cpi];
            const SpacePointInfo &pi = p_info[cpn.id];

            pi.u[ln] += p[cpn.j][cpn.i] * (cpn.w * (hx*hy));

            pi.ux[ln] = (p[pi.j][pi.i+1] - p[pi.j][pi.i-1])/(2.0*hx);
            pi.uy[ln] = (p[pi.j+1][pi.i] - p[pi.j-1][pi.i])/(2.0*hy);

            //ui.ux[ln] = (u[ui.j][ui.i-2]-8.0*u[ui.j][ui.i-1]+8.0*u[ui.j][ui.i+1]-u[ui.j][ui.i+2])/(12.0*hx);
            //ui.uy[ln] = (u[ui.j-2][ui.i]-8.0*u[ui.j-1][ui.i]+8.0*u[ui.j+1][ui.i]-u[ui.j+2][ui.i])/(12.0*hy);
        }
    }

    if (method == 4)
    {
        for (unsigned int j=0; j<p_info.size(); j++) p_info[j].u[ln] = p_info[j].ux[ln] = p_info[j].uy[ln] = 0.0;

        for (unsigned int cpi=0; cpi<cntPointNodes.size(); cpi++)
        {
            const ExtendedSpacePointNode &cpn = cntPointNodes[cpi];
            const SpacePointInfo &pi = p_info[cpn.id];

            pi.u[ln] += p[cpn.j][cpn.i] * (cpn.w * (hx*hy));

            pi.ux[ln] = (p[pi.j][pi.i+1] - p[pi.j][pi.i-1])/(2.0*hx);
            pi.uy[ln] = (p[pi.j+1][pi.i] - p[pi.j-1][pi.i])/(2.0*hy);

            //pi.ux[ln] = (p[ui.j][ui.i-2]-8.0*u[ui.j][ui.i-1]+8.0*u[ui.j][ui.i+1]-u[ui.j][ui.i+2])/(12.0*hx);
            //pi.uy[ln] = (u[ui.j-2][ui.i]-8.0*u[ui.j-1][ui.i]+8.0*u[ui.j+1][ui.i]-u[ui.j+2][ui.i])/(12.0*hy);
        }
    }
}

void Problem2HNDirichlet::b_layerInfo(const DoubleMatrix &p UNUSED_PARAM, unsigned int ln UNUSED_PARAM) const
{
    //double min = p.min();
    //double max = p.max();

    //QPixmap pic;
    //visualizeMatrixHeat(p, min, max, pic);
    //pic.save("images/b/pic"+QString("%1").arg(ln)+".png", "PNG");
    //printf("Layer: %d min: %f max: %f min: %f max: %f norm: %f\n", ln, min, max, min, max, fabs(max-min));
}

// backward -----------------------------------

void Problem2HNDirichlet::distributeDelta0(const SpacePoint &pt, unsigned int id, espn_vector &nodes, const Dimension &dimX, const Dimension &dimY, unsigned int k, int method) const
{
    if (method == 1) distributeDeltaP(pt, id, nodes, dimX, dimY);

    if (method == 2) distributeDeltaR(pt, id, nodes, dimX, dimY);

    if (method == 4) distributeDeltaG(pt, id, nodes, dimX, dimY, k);
}

void Problem2HNDirichlet::distributeDeltaP(const SpacePoint &pt, unsigned int id, espn_vector &nodes, const Dimension &dimX, const Dimension &dimY, unsigned int) const
{
    double hx = dimX.step();
    double hy = dimY.step();

    unsigned int Nx = dimX.sizeN();
    unsigned int Ny = dimY.sizeN();

    unsigned int rx = (unsigned int) (round( pt.x * Nx ));
    unsigned int ry = (unsigned int) (round( pt.y * Ny ));

    ExtendedSpacePointNode node; node.id = id; node.pt = pt; node.i = rx; node.j = ry; node.w = 1.0/(hx*hy); nodes.push_back(node);
}

void Problem2HNDirichlet::distributeDeltaR(const SpacePoint &pt, unsigned int id, espn_vector &nodes, const Dimension &dimX, const Dimension &dimY, unsigned int) const
{
    double hx = dimX.step();
    double hy = dimY.step();

    unsigned int Nx = dimX.sizeN();
    unsigned int Ny = dimY.sizeN();

    unsigned int rx = (unsigned int) (ceil( pt.x * Nx ));
    unsigned int ry = (unsigned int) (ceil( pt.y * Ny ));

    double h1x = fabs(pt.x - rx*hx);
    double h1y = fabs(pt.y - ry*hy);
    double h2x = hx-fabs(pt.x - rx*hx);
    double h2y = hx-fabs(pt.y - ry*hy);

    ExtendedSpacePointNode node00; node00.id = id; node00.pt = pt; node00.i = rx+0; node00.j = ry+0; node00.w = ((h2x/hx)*(h2y/hy))/(hx*hy); nodes.push_back(node00);
    ExtendedSpacePointNode node01; node01.id = id; node01.pt = pt; node01.i = rx+0; node01.j = ry+1; node01.w = ((h2x/hx)*(h1y/hy))/(hx*hy); nodes.push_back(node01);
    ExtendedSpacePointNode node11; node11.id = id; node11.pt = pt; node11.i = rx+1; node11.j = ry+1; node11.w = ((h1x/hx)*(h1y/hy))/(hx*hy); nodes.push_back(node11);
    ExtendedSpacePointNode node10; node10.id = id; node10.pt = pt; node10.i = rx+1; node10.j = ry+0; node10.w = ((h1x/hx)*(h2y/hy))/(hx*hy); nodes.push_back(node10);
}

void Problem2HNDirichlet::distributeDeltaG(const SpacePoint &pt, unsigned int id, espn_vector &nodes, const Dimension &dimX, const Dimension &dimY, unsigned int k) const
{
    double hx = dimX.step();
    double hy = dimY.step();

    unsigned int Nx = dimX.sizeN();
    unsigned int Ny = dimY.sizeN();

    unsigned int rx = (unsigned int) (round(pt.x*Nx));
    unsigned int ry = (unsigned int) (round(pt.y*Ny));

    double sigmaX = hx;
    double sigmaY = hy;

    double sumX = 0.0;
    for (unsigned int n=rx-k; n<=rx+k; n++) sumX += exp(-((n*hx-pt.x)*(n*hx-pt.x))/(2.0*sigmaX*sigmaX));
    sumX *= hx;

    double sumY = 0.0;
    for (unsigned int m=ry-k; m<=ry+k; m++) sumY += exp(-((m*hy-pt.y)*(m*hy-pt.y))/(2.0*sigmaY*sigmaY));
    sumY *= hy;

    double sigma = (sumX*sumY) / (2.0*M_PI);
    double factor = 1.0/((2.0*M_PI)*sigma);

    for (unsigned int m=ry-k; m<=ry+k; m++)
    {
        for (unsigned int n=rx-k; n<=rx+k; n++)
        {
            ExtendedSpacePointNode node;
            node.i = n; node.x = n*hx;
            node.j = m; node.y = m*hy;
            node.pt = pt; node.id = id;
            node.w = factor*exp(-0.5*(((node.x-pt.x)*(node.x-pt.x))/(sigmaX*sigmaX)+((node.y-pt.y)*(node.y-pt.y))/(sigmaY*sigmaY)));
            nodes.push_back(node);
        }
    }
}

void Problem2HNDirichlet::PrmToVector(const OptimizeParameter &prm, DoubleVector &pv) const
{
    unsigned int Nc = mEquParameter.Nc;
    unsigned int No = mEquParameter.No;

    pv.clear();
    pv.resize(2*Nc*No+2*No+2*Nc);

    for (unsigned int i=0; i<Nc; i++)
    {
        for (unsigned int j=0; j<No; j++)
        {
            pv[i*No + j] = prm.k[i][j];
        }
    }

    for (unsigned int i=0; i<Nc; i++)
    {
        for (unsigned int j=0; j<No; j++)
        {
            pv[i*No + j + Nc*No] = prm.z[i][j];
        }
    }


    for (unsigned int j=0; j<No; j++)
    {
        pv[2*j + 0 + 2*Nc*No] = prm.xi[j].x;
        pv[2*j + 1 + 2*Nc*No] = prm.xi[j].y;
    }



    for (unsigned int i=0; i<Nc; i++)
    {
        pv[2*i + 0 + 2*No + 2*Nc*No] = prm.eta[i].x;
        pv[2*i + 1 + 2*No + 2*Nc*No] = prm.eta[i].y;
    }
}

void Problem2HNDirichlet::VectorToPrm(const DoubleVector &pv, OptimizeParameter &prm) const
{
    unsigned int Nc = mEquParameter.Nc;
    unsigned int No = mEquParameter.No;

    unsigned int index = 0;

    prm.k.clear();
    prm.k.resize(Nc, No);

    for (unsigned int i=0; i<Nc; i++)
    {
        for (unsigned int j=0; j<No; j++)
        {
            prm.k[i][j] = pv[index]; index++;
        }
    }

    prm.z.clear();
    prm.z.resize(Nc, No);

    for (unsigned int i=0; i<Nc; i++)
    {
        for (unsigned int j=0; j<No; j++)
        {
            prm.z[i][j] = pv[index]; index++;
        }
    }

    prm.xi.clear();
    prm.xi.resize(No);

    for (unsigned int j=0; j<No; j++)
    {
        prm.xi[j].x = pv[index]; index++;
        prm.xi[j].y = pv[index]; index++;
    }

    prm.eta.clear();
    prm.eta.resize(Nc);

    for (unsigned int i=0; i<Nc; i++)
    {
        prm.eta[i].x = pv[index]; index++;
        prm.eta[i].y = pv[index]; index++;
    }
}

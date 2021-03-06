#include "art_problem1.h"

#define RMAX 0.0

void ArtProblem1::Main(int argc UNUSED_PARAM, char *argv[] UNUSED_PARAM)
{
    ArtProblem1 p;
    p.initialize();
    p.startOptimize();

    //p.table2Generate();
    //p.image1L();
    //p.imager2L();
    //p.imager3L();

    //p.image1Generate();
    //p.image2Generate();
    //p.image3Generate();
}

ArtProblem1::ArtProblem1() {}

void ArtProblem1::initialize()
{
    optimizeK = true;
    optimizeZ = true;
    optimizeE = true;
    withError = false;
    DD = 1.0;

    L = 2;

    N = 1000;
    hx = 0.001;

    M = 1000;
    ht = 0.001;

    // initial temperatures
    vfi << 0.0;// << 1.0 << 2.0;
    // environment temperatures
    vtt << +5.0;// << +6.0 << +7.0;

    // initial temperature
    fi = 0.0;
    // environment temperature
    tt = 0.0;

    // золото
    //    {
    //        const double r0 = 19320.0; // kg/m^3   // плотность
    //        const double c0 = 130.0;   // C/(kg*S) // удельная теплоемкость
    //        const double k0 = 312.0;   //          // коэффициент теплопроводности

    //        const double h1 = 1000.0;      // коэффициент теплообмена ->
    //        const double h2 = 10.0;        // коэффициент теплообмена ->

    //        a = sqrt((k0/(c0*r0)));        // коэффициент температуропроворности
    //        lambda0 = h2/(c0*r0);          // коэффициент теплообмена ->
    //        lambda1 = h1/k0;               // коэффициент теплообмена ->
    //        lambda2 = h2/k0;               // коэффициент теплообмена ->
    //    }

    // мед
    //    {
    //        const double r0 = 8900.0;  // kg/m^3     // плотность
    //        const double c0 = 400.0;   // C/(kg*S)   // удельная теплоемкость
    //        const double k0 = 380.0;   // Vt/(m*S)   // коэффициент теплопроводности

    //        const double h1 = 1000.0;      // коэффициент теплообмена ->
    //        const double h2 = 10.0;        // коэффициент теплообмена ->

    //        a = sqrt((k0/(c0*r0)));        // коэффициент температуропроворности
    //        lambda0 = (h2/(c0*r0));          // коэффициент теплообмена ->
    //        lambda1 = (h1/k0);               // коэффициент теплообмена ->
    //        lambda2 = (h2/k0);               // коэффициент теплообмена ->
    //    }

    a = 1.0;
    lambda0 = 0.001;
    lambda1 = 100.0;
    lambda2 = 0.010;


    /* коэффициенты регуляризации */
    alpha0 = 1.0;
    alpha1 = 0.0000;
    alpha2 = 0.0000;
    alpha3 = 0.0000;

    //    if (!optimizeK) alpha1 = 0.0;
    //    if (!optimizeZ) alpha2 = 0.0;
    //    if (!optimizeE) alpha3 = 0.0;

    if (L==1)
    {
        k0 << 0.0;
        z0 << 0.0;
        e0 << 0.0;
    }

    if (L==2)
    {
        k0 << 0.0 << 0.0;
        z0 << 0.0 << 0.0;
        e0 << 0.0 << 0.0;
    }

    if (L==3)
    {
        k0 << 0.0 << 0.0 << 0.0;
        z0 << 0.0 << 0.0 << 0.0;
        e0 << 0.0 << 0.0 << 0.0;
    }

    if (L==5)
    {
        k0 << 0.0 << 0.0 << 0.0 << 0.0 << 0.0;
        z0 << 0.0 << 0.0 << 0.0 << 0.0 << 0.0;
        e0 << 0.0 << 0.0 << 0.0 << 0.0 << 0.0;
    }

    if (L==10)
    {
        k0 << 0.0 << 0.0 << 0.0 << 0.0 << 0.0 << 0.0 << 0.0 << 0.0 << 0.0 << 0.0;
        z0 << 0.0 << 0.0 << 0.0 << 0.0 << 0.0 << 0.0 << 0.0 << 0.0 << 0.0 << 0.0;
        e0 << 0.0 << 0.0 << 0.0 << 0.0 << 0.0 << 0.0 << 0.0 << 0.0 << 0.0 << 0.0;
    }

    /* шаги числовых производных */
    hk = 0.001;
    hz = 0.001;
    he = 0.001;

    R = 1.0;

    /* пределы z параметров */
    vmin = +10.0;
    vmax = +60.0;
    d0 = (vmax+vmin)/2.0;
    d1 = (vmax-vmin)/2.0;

    /* температура стержня */
    V.resize(N+1);
    for (unsigned int n=0; n<=N; n++) V[n] = 10.0;// - sin(2.0*M_PI*(n*hx))*sin(2.0*M_PI*(n*hx));
    IPrinter::printVector(V);
}

void ArtProblem1::startOptimize()
{
    DoubleVector x0;

    if (L==1)
    {
        if (optimizeK)
        {
            x0 << -4.5433;
            //x0 << -5.2547;
            //x0 << -10.2547;

            //x0 << -12.8746;
            //x0 << -8.8745;
        }
        if (optimizeZ)
        {
            x0 << 9.4343;
            //x0 << 8.2545;
            //x0 << 19.2545;

            //x0 << 15.8272;
            //x0 << -5.2481;
        }
        if (optimizeE)
        {
            x0 << 0.5435;
            //x0 << 0.1245;
            //x0 << 0.3245;

            //x0 << 0.7358;
            //x0 << 0.8745;
        }
    }

    if (L==2)
    {
        if (optimizeK)
        {
            //x0 << -5.5400 << -8.4500; //k
            //x0 << -8.1284 << -1.4928; //k
            //x0 << -2.8525 << -6.8774; //k
            x0 << -12.1248 << -4.4918; //k
            //x0 << -1.8142 << -2.9184; //k
            //x0 << -15.5744 << -8.4324;
        }
        else
        {
            K.clear();
            K << -8.5000 << -2.7000; //k
        }

        if (optimizeZ)
        {
            //x0 << +10.1500 << +12.9600; //z
            //x0 << +15.7465 << +7.0645; //k
            //x0 << +2.4548 << +2.4518; //k
            x0 << +20.8248 << +17.7545; //k
            //x0 << +10.4755 << +11.8428; //k
            //x0 << 15.4783 << 8.6603;
        }
        else
        {
            Z.clear();
            Z << +2.1000 << +4.9000; //z
        }

        if (optimizeE)
        {
            //x0 << 0.2500 << 0.7500; //e
            //x0 << 0.1055 << 0.8155; //e
            //x0 << 0.5587 << 0.7545; //e
            x0 << 0.2095 << 0.9511; //e
            //x0 << 0.5551 << 0.4915; //e
            //x0 << 0.2645 << 0.7832;
        }
        else
        {
            E.clear();
            E << +0.02000 << +0.08000; //e
        }
    }

    if (L==3)
    {
        if (optimizeK)
        {
            x0 << -1.1284 << -1.4928 << -2.5478; //k
            //x0 << -4.8525 << -6.8774 << -4.1428; //k
            //x0 << -4.4157 << -5.7698 << -3.1589; //k
            //x0 << -10.7845 << -8.8674 << -1.2457; //k
            //x0 << -5.8417 << -20.4157 << -14.1541; //k
        }

        if (optimizeZ)
        {
            x0 << +9.7465 << +8.0645 << +6.2471; //k
            //x0 << +4.5448 << +5.1852 << +8.2451;  //k
            //x0 << +20.874 << 18.2578 << +11.8428; //k
            //x0 << +11.4751 << +12.8541 << +11.2541; //z
            //x0 << +10.2568 << +8.4157 << 10.2144; //k
        }

        if (optimizeE)
        {
            x0 << 0.5517 << 0.8155 << 0.974; //e
            //x0 << 0.5587 << 0.7545 << 0.2145; //e
            //x0 << 0.9541 << 0.1225 << 0.5841; //e
            //x0 << 0.3541 << 0.6648 << 0.7574; //e
            //x0 << 0.2095 << 0.5142 << 0.9474; //e
        }
    }

    if (L==5)
    {
        if (optimizeK)
        {
            x0 << -1.1284 << -1.4928 << -2.5478 << -5.3234 << +0.0322; //k
        }

        if (optimizeZ)
        {
            x0 << +9.7465 << +8.0645 << +6.2471 << +9.2231 << +3.4322; //z
        }

        if (optimizeE)
        {
            x0 << +0.0050 << +0.2500 << +0.5000 << +0.7500 << +0.9950; //e
        }
    }

    if (L==10)
    {
        if (optimizeK)
        {
            x0 << -1.1284 << -1.4928 << -2.5478 << -5.3234 << +0.0322 << -8.1842 << -8.4892 << -3.5784 << -1.3342 << +1.0242; //k
        }

        if (optimizeZ)
        {
            x0 << +9.7465 << +8.0645 << +6.2471 << +9.2231 << +3.4322 << +5.7946 << +6.0485 << +7.2641 << +1.2923 << +2.4332; //z
        }

        if (optimizeE)
        {
            x0 << +0.1000 << +0.2000 << +0.3000 << +0.4000 << +0.5000 << +0.6000 << +0.7000 << +0.8000 << +0.9000 << +0.9950; //e
        }
    }

    DD = 1.0;
    R = 1.0;
    optimize(x0);
    while (R < 1000000.0)
    {
        R *= 2.0;
        IPrinter::printSeperatorLine();
        optimize(x0);
    }

    double y = fx(x0);

    IPrinter::printSeperatorLine();
    printf("Optimal k: "); for (unsigned int i=0*L; i<1*L; i++) { printf("%20.14f ", x0[i]); } printf("\n");
    printf("Optimal z: "); for (unsigned int i=1*L; i<2*L; i++) { printf("%20.14f ", x0[i]); } printf("\n");
    printf("Optimal e: "); for (unsigned int i=2*L; i<3*L; i++) { printf("%20.14f ", x0[i]); } printf("\n");
    IPrinter::printSeperatorLine();
}

void ArtProblem1::optimize(DoubleVector &x0) const
{
    ArtProblem1* p = const_cast<ArtProblem1*>(this);
    ConjugateGradient g;
    g.setFunction(p);
    g.setGradient(p);
    g.setPrinter(p);
    g.setProjection(p);
    g.setOptimalityTolerance(0.0001);//0.00000001
    g.setStepTolerance(0.0001);//0.00000001
    g.setFunctionTolerance(0.0001);//0.00000001
    g.setR1MinimizeEpsilon(10.0, 0.0001); //0.00000001
    g.setNormalize(true);
    g.showExitMessage(true);
    g.setResetIteration(false);
    g.calculate(x0);
}

void ArtProblem1::project(DoubleVector &x, size_t i)
{
    /* e lower/upper limits */
    if (optimizeE)
    {
        if (L==1)
        {
            if (i==2) { if (x.at(2) < 0.005) x.at(2) = 0.005; if (x.at(2) > 0.995) x.at(2) = 0.995; }
        }

        if (L==2)
        {
            if (i==4) { if (x.at(4) < 0.005) x.at(4) = 0.005; if (x.at(4) > 0.995) x.at(4) = 0.995; }
            if (i==5) { if (x.at(5) < 0.005) x.at(5) = 0.005; if (x.at(5) > 0.995) x.at(5) = 0.995; }
        }

        if (L==3)
        {
            if (i == 6) { if (x.at(6) < 0.005) x.at(6) = 0.005; if (x.at(6) > 0.995) x.at(6) = 0.995; }
            if (i == 7) { if (x.at(7) < 0.005) x.at(7) = 0.005; if (x.at(7) > 0.995) x.at(7) = 0.995; }
            if (i == 8) { if (x.at(8) < 0.005) x.at(8) = 0.005; if (x.at(8) > 0.995) x.at(8) = 0.995; }
        }

        if (L==5)
        {
            if (i == 10) { if (x.at(10) < 0.005) x.at(10) = 0.005; if (x.at(10) > 0.995) x.at(10) = 0.995; }
            if (i == 11) { if (x.at(11) < 0.005) x.at(11) = 0.005; if (x.at(11) > 0.995) x.at(11) = 0.995; }
            if (i == 12) { if (x.at(12) < 0.005) x.at(12) = 0.005; if (x.at(12) > 0.995) x.at(12) = 0.995; }
            if (i == 13) { if (x.at(13) < 0.005) x.at(13) = 0.005; if (x.at(13) > 0.995) x.at(13) = 0.995; }
            if (i == 14) { if (x.at(14) < 0.005) x.at(14) = 0.005; if (x.at(14) > 0.995) x.at(14) = 0.995; }
        }

        if (L==10)
        {
            if (i == 20) { if (x.at(20) < 0.005) x.at(20) = 0.005; if (x.at(20) > 0.995) x.at(20) = 0.995; }
            if (i == 21) { if (x.at(21) < 0.005) x.at(21) = 0.005; if (x.at(21) > 0.995) x.at(21) = 0.995; }
            if (i == 22) { if (x.at(22) < 0.005) x.at(22) = 0.005; if (x.at(22) > 0.995) x.at(22) = 0.995; }
            if (i == 23) { if (x.at(23) < 0.005) x.at(23) = 0.005; if (x.at(23) > 0.995) x.at(23) = 0.995; }
            if (i == 24) { if (x.at(24) < 0.005) x.at(24) = 0.005; if (x.at(24) > 0.995) x.at(24) = 0.995; }
            if (i == 25) { if (x.at(25) < 0.005) x.at(25) = 0.005; if (x.at(25) > 0.995) x.at(25) = 0.995; }
            if (i == 26) { if (x.at(26) < 0.005) x.at(26) = 0.005; if (x.at(26) > 0.995) x.at(26) = 0.995; }
            if (i == 27) { if (x.at(27) < 0.005) x.at(27) = 0.005; if (x.at(27) > 0.995) x.at(27) = 0.995; }
            if (i == 28) { if (x.at(28) < 0.005) x.at(28) = 0.005; if (x.at(28) > 0.995) x.at(28) = 0.995; }
            if (i == 29) { if (x.at(29) < 0.005) x.at(29) = 0.005; if (x.at(29) > 0.995) x.at(29) = 0.995; }
        }
    }
}

void ArtProblem1::image1L()
{
    DD = 1;
    DoubleVector y0;
    optimizeK = true;
    optimizeZ = true;
    optimizeE = true;
    withError = false;

    if (optimizeK)
    {
        y0 << 4.9688708458;
        //        y0 << -4.8090586524;
        //        y0 << -6.8182056823;
    }

    if (optimizeZ)
    {
        y0 << 8.2385291524;
        //        y0 << 12.0785195749;
        //        y0 << 11.9038403093;
    }
    if (optimizeE)
    {
        y0 << 0.0186136703;
        //        y0 << 0.6253314029;
        //        y0 << 0.0782326627;
    }

    L = 1;
    R = 1.0;
    //y0[2] = 19*hx;
    //    double y = fx(y0);
    //    printf("%.10f %.10f\n", y, y0[2]);

    FILE *file = fopen("L1Data.txt", "w");
    for (unsigned int n=5; n<=995; n++)
    {
        y0[2] = n*hx;
        double f = fx(y0);
        printf("%d %.10f %f\n", n, f, y0[2]);
        fprintf(file, "%.10f ", f);
    }
    fclose(file);
}

void ArtProblem1::imager2L()
{
    DoubleVector y0;
    optimizeK = true;
    optimizeZ = true;
    optimizeE = true;
    withError = false;
    DD = 1;
    L=2;

    if (optimizeK)
    {
        y0 << -6.0341111852 << 1.1251024127;
    }

    if (optimizeZ)
    {
        y0 << 12.5876519120 << 14.9819901001;
    }
    if (optimizeE)
    {
        y0 << 0.5889625969 << 0.7277517979;
    }

    DoubleMatrix u(101,101);
    R = 1.0;
    double max = -100000000.0;
    double min = +100000000.0;
    for (unsigned int k0=0; k0<=100; k0++)
    {
        y0[4] = 0.01*k0;
        //        y0[4] = 0.01*k0;
        for (unsigned int k1=0; k1<=100; k1++)
        {
            y0[5] = 0.01*k1;
            //            y0[5] = 0.01*k1;

            double f = 0.0;
            f = fx(y0);
            u[k0][k1] = f;

            if (f > max) max = f;
            if (f < min) min = f;

            printf("%d %d %f\n", k0, k1, f);
        }
    }

    for (unsigned int i=0; i<=100; i++) u[0][i] = u[2][i];
    for (unsigned int i=0; i<=100; i++) u[1][i] = u[2][i];
    for (unsigned int i=0; i<=100; i++) u[99][i] = u[98][i];
    for (unsigned int i=0; i<=100; i++) u[100][i] = u[98][i];
    for (unsigned int i=0; i<=100; i++) u[i][0] = u[i][2];
    for (unsigned int i=0; i<=100; i++) u[i][1] = u[i][2];
    for (unsigned int i=0; i<=100; i++) u[i][99] = u[i][98];
    for (unsigned int i=0; i<=100; i++) u[i][100] = u[i][98];

    printf("%f %f\n", min, max);

    FILE *file = fopen("data1.txt", "w");
    IPrinter::printMatrix(u,u.rows(),u.cols(),NULL,file);
    fclose(file);
}

void ArtProblem1::imager3L()
{
    DoubleVector y0;
    optimizeK = true;
    optimizeZ = true;
    optimizeE = true;

    if (optimizeK)
    {
        y0 << -1.51758851605847 << -3.38265533950958 << 0.15272285434824;
    }
    if (optimizeZ)
    {
        y0 << +11.23656775924430 << +13.09747632971720 << +8.95446642160735;
    }
    if (optimizeE)
    {
        y0 << +0.03200000000000 << +0.06500000000000 << +0.06700000000000;
    }

    DoubleMatrix u(101,101);
    R = 0.0;
    double max = -100000000.0;
    double min = +100000000.0;
    for (unsigned int k0=2; k0<=98; k0++)
    {
        y0[7] = 0.001*k0;
        for (unsigned int k1=2; k1<=98; k1++)
        {
            y0[8] = 0.001*k1;

            double f = 0.0;
            if (k0 == k1 || k0 == k1-1 || k0 == k1+1 || k0 == 32 || k1 == 32)
            {
                f = 0.0;
            }
            else
            {
                f = fx(y0);
            }
            u[k0][k1] = f;


            if (f > max) max = f;
            if (f < min) min = f;

            printf("%d %d %f\n", k0, k1, f);
        }
    }
    printf("%f %f\n", min, max);

    FILE *file = fopen("data1.txt", "w");
    IPrinter::printMatrix(u,u.rows(),u.cols(),NULL,file);
    fclose(file);
}

void ArtProblem1::table1Generate()
{
    DoubleVector y0;
    optimizeK = true;
    optimizeZ = true;
    optimizeE = true;

    if (optimizeK)
    {
        y0 << -8.50 << -2.70;
    }
    if (optimizeZ)
    {
        y0 << +2.10 << +4.90;
    }
    if (optimizeE)
    {
        y0 << +0.020 << +0.080;
    }
    R = 1.0;
    optimize(y0);
    while (R < 10000.0)
    {
        IPrinter::printSeperatorLine();
        R *= 10.0;
        optimize(y0);
    }

    DoubleVector k,z,e;
    getParameters(k,z,e,y0);
    IPrinter::printSeperatorLine();
    printf("Optimal k: %20.14f %20.14f %20.14f\n", k[0], k[1], k[2]);
    printf("Optimal z: %20.14f %20.14f %20.14f\n", z[0], z[1], z[2]);
    printf("Optimal e: %20.14f %20.14f %20.14f\n", e[0], e[1], e[2]);
    IPrinter::printSeperatorLine();
}

void ArtProblem1::table2Generate()
{
    DoubleVector y0;
    optimizeK = true;
    optimizeZ = true;
    optimizeE = true;

    if (optimizeK)
    {
        y0 << -4.54 << -7.45;// << -3.500;;
    }
    if (optimizeZ)
    {
        y0 << +12.10 << +14.90;// << +10.13;
    }
    if (optimizeE)
    {
        y0 << 0.02500 << 0.04500;// << 0.07500;
    }

    DD = 10;

    R = 1.0;
    optimize(y0);
    while (R < 10000.0)
    {
        IPrinter::printSeperatorLine();
        R *= 10.0;
        optimize(y0);
    }

    DoubleVector k,z,e;
    getParameters(k,z,e,y0);
    IPrinter::printSeperatorLine();

    if (L==2)
    {
        printf("Optimal k: %20.14f %20.14f\n", k[0], k[1]);
        printf("Optimal z: %20.14f %20.14f\n", z[0], z[1]);
        printf("Optimal e: %20.14f %20.14f\n", e[0], e[1]);
    }

    if (L==3)
    {
        printf("Optimal k: %20.14f %20.14f %20.14f\n", k[0], k[1], k[2]);
        printf("Optimal z: %20.14f %20.14f %20.14f\n", z[0], z[1], z[2]);
        printf("Optimal e: %20.14f %20.14f %20.14f\n", e[0], e[1], e[2]);
    }

    IPrinter::printSeperatorLine();

    DoubleMatrix u;
    calculateU(u, k, z, e);
    IPrinter::printVector(14,10,u.row(u.rows()-1));
}

void ArtProblem1::image1Generate()
{
    DoubleVector y0;
    optimizeK = true;
    optimizeZ = true;
    optimizeE = true;

    if (optimizeK)
    {
        y0 << -6.0341111852 << 1.1251024127;
    }

    if (optimizeZ)
    {
        y0 << 12.5876519120 << 14.9819901001;
    }
    if (optimizeE)
    {
        y0 << 0.5889625969 << 0.7277517979;
    }

    DD = 1;

    DoubleVector k,z,e;
    getParameters(k,z,e,y0);

    //    M *= 2;

    FILE *file1 = fopen("image_1_v.txt", "w");

    vfi.clear();
    vtt.clear();

    fi = +0.0;
    tt = +5.0;
    DoubleMatrix u1;
    calculateU(u1,k,z,e);
    for (int unsigned m=0; m<=M; m++)
    {
        double v = vf(m,k,z,e,u1);
        fprintf(file1, "%.10f ",v);
        fflush(file1);
    }
    fprintf(file1, "\n");

    fi = -2.0;
    tt = +9.0;
    DoubleMatrix u2;
    calculateU(u2,k,z,e);
    for (int unsigned m=0; m<=M; m++)
    {
        double v = vf(m,k,z,e,u2);
        fprintf(file1, "%.10f ",v);
        fflush(file1);
    }
    fprintf(file1, "\n");

    fclose(file1);

    //    M=1000;
}

void ArtProblem1::image2Generate()
{
    DoubleVector y0;
    optimizeK = true;
    optimizeZ = true;
    optimizeE = true;

    if (optimizeK)
    {
        y0 << -6.0341111852 << 1.1251024127;
    }

    if (optimizeZ)
    {
        y0 << 12.5876519120 << 14.9819901001;
    }
    if (optimizeE)
    {
        y0 << 0.5889625969 << 0.7277517979;
    }

    DoubleVector k,z,e;
    getParameters(k,z,e,y0);

    FILE *file1 = fopen("image_2_du.txt", "w");
    //    M *= 2;

    vfi.clear();
    vtt.clear();

    fi = +0.0;
    tt = +5.0;
    DoubleMatrix u1;
    calculateU(u1,k,z,e);
    for (int unsigned m=0; m<=M; m++)
    {
        double max = -1000000.0;
        DoubleVector ut = u1.row(m);
        for (unsigned int n=0; n<=N; n++)
        {
            if (max < fabs(ut[n]-V[n])) max = fabs(ut[n]-V[n]);
        }
        fprintf(file1, "%.10f ",max);
        fflush(file1);
    }
    fprintf(file1, "\n");

    fi = -2.0;
    tt = +9.0;
    DoubleMatrix u2;
    calculateU(u2,k,z,e);
    for (int unsigned m=0; m<=M; m++)
    {
        double max = -1000000.0;
        DoubleVector ut = u2.row(m);
        for (unsigned int n=0; n<=N; n++)
        {
            if (max < fabs(ut[n]-V[n])) max = fabs(ut[n]-V[n]);
        }
        fprintf(file1, "%.10f ",max);
        fflush(file1);
    }
    fprintf(file1, "\n");

    fclose(file1);
    //    M=1000;
}

void ArtProblem1::image3Generate()
{
    DoubleVector y0;
    optimizeK = true;
    optimizeZ = true;
    optimizeE = true;

    if (optimizeK)
    {
        y0 << -6.0341111852 << 1.1251024127;
    }

    if (optimizeZ)
    {
        y0 << 12.5876519120 << 14.9819901001;
    }
    if (optimizeE)
    {
        y0 << 0.5889625969 << 0.7277517979;
    }

    DoubleVector k,z,e;
    getParameters(k,z,e,y0);

    vfi.clear();
    vtt.clear();

    withError = true;
    persent = 0.05;
    FILE *file1 = fopen("image_3_du_5.txt", "w");

    fi = +0.0;
    tt = +5.0;
    M *= 2;
    DoubleMatrix u1;
    calculateU(u1,k,z,e);
    for (int unsigned m=0; m<=M; m++)
    {
        DoubleVector ut = u1.row(m);
        double max = -1000000.0;
        for (unsigned int n=0; n<=N; n++)
        {
            if (max < fabs(ut[n]-V[n])) max = fabs(ut[n]-V[n]);
        }

        fprintf(file1, "%.10f ",max);
        fflush(file1);
    }
    fprintf(file1, "\n");
    fclose(file1);
}

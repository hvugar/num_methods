#include "art_problem1L2.h"

#define RMAX 0.0

void ArtProblem1L2::Main(int argc UNUSED_PARAM, char *argv[] UNUSED_PARAM)
{
    ArtProblem1L2 p;
    p.initialize();
    p.startOptimize();

    //p.table2Generate();
    //p.imager3L();
    //p.imager2L();

//    p.image1Generate();
//    p.image2Generate();
//    p.image3Generate();
}

ArtProblem1L2::ArtProblem1L2() {}

void ArtProblem1L2::initialize()
{
    optimizeK = true;
    optimizeZ = true;
    optimizeE = true;
    withError = false;
    DD = 1;

    L = 2;

    N = 1000;
    hx = 0.001;

    M = 1000;
    ht = 0.001;

    // initial temperatures
    vfi << 0.0 << 1.0 << 2.0;
    // environment temperatures
    vtt << +5.0 << +6.0 << +7.0;

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
    {
        const double r0 = 8900.0;  // kg/m^3     // плотность
        const double c0 = 400.0;   // C/(kg*S)   // удельная теплоемкость
        const double k0 = 380.0;   // Vt/(m*S)   // коэффициент теплопроводности

        const double h1 = 1000.0;      // коэффициент теплообмена ->
        const double h2 = 10.0;        // коэффициент теплообмена ->

        a = sqrt((k0/(c0*r0)));        // коэффициент температуропроворности
        lambda0 = (h2/(c0*r0));          // коэффициент теплообмена ->
        lambda1 = (h1/k0);               // коэффициент теплообмена ->
        lambda2 = (h2/k0);               // коэффициент теплообмена ->
    }

    a = 1.0;
    lambda0 = 0.001;
    lambda1 = 1000.0;
    lambda2 = 0.0010;


    /* коэффициенты регуляризации */
    alpha0 = 1.0;
    alpha1 = 0.0000;
    alpha2 = 0.0000;
    alpha3 = 0.0000;

    if (!optimizeK) alpha1 = 0.0;
    if (!optimizeZ) alpha2 = 0.0;
    if (!optimizeE) alpha3 = 0.0;

    if (L==2)
    {
        k0 << 0.0 << 0.0;//-20.57571017653454 << -30.63314593795166;
        z0 << 0.0 << 0.0;//+10.33818417154749 << +10.47968970008047;
        e0 << 0.0 << 0.0;// +0.04500000000000 <<  +0.09500000000000;
    }

    if (L==3)
    {
        k0 << 0.0 << 0.0 << 0.0;
        z0 << 0.0 << 0.0 << 0.0;
        e0 << 0.0 << 0.0 << 0.0;
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
    for (unsigned int n=0; n<=N; n++) V[n] = 10.0;
    IPrinter::printVector(V);
}

void ArtProblem1L2::table1Generate()
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

void ArtProblem1L2::table2Generate()
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

void ArtProblem1L2::image1Generate()
{
    DoubleVector y0;
    optimizeK = true;
    optimizeZ = true;
    optimizeE = true;

    if (optimizeK)
    {
        y0 << -0.103894  <<  -1.837793;
    }

    if (optimizeZ)
    {
        y0 << 13.063753 << 15.239693;
    }
    if (optimizeE)
    {
        y0 << 0.067496 << 0.793735;
    }

    DD = 1;

    DoubleVector k,z,e;
    getParameters(k,z,e,y0);

    M *= 2;

    FILE *file1 = fopen("image_1_v.txt", "w");

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

    fi = -3.0;
    tt = +5.0;
    DoubleMatrix u2;
    calculateU(u2,k,z,e);
    for (int unsigned m=0; m<=M; m++)
    {
        double v = vf(m,k,z,e,u2);
        fprintf(file1, "%.10f ",v);
        fflush(file1);
    }
    fprintf(file1, "\n");

    fi = +2.0;
    tt = +5.0;
    DoubleMatrix u3;
    calculateU(u3,k,z,e);
    for (int unsigned m=0; m<=M; m++)
    {
        double v = vf(m,k,z,e,u3);
        fprintf(file1, "%.10f ",v);
        fflush(file1);
    }

    fclose(file1);

    M=1000;
}

void ArtProblem1L2::image2Generate()
{
    DoubleVector y0;
    optimizeK = true;
    optimizeZ = true;
    optimizeE = true;

    if (optimizeK)
    {
        y0 << -0.103894  <<  -1.837793;
    }

    if (optimizeZ)
    {
        y0 << 13.063753 << 15.239693;
    }
    if (optimizeE)
    {
        y0 << 0.067496 << 0.793735;
    }

    DoubleVector k,z,e;
    getParameters(k,z,e,y0);

    FILE *file1 = fopen("image_2_du.txt", "w");
    M *= 2;

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

    fi = -3.0;
    tt = +5.0;
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

    fi = +2.0;
    tt = +5.0;
    DoubleMatrix u3;
    calculateU(u3,k,z,e);
    for (int unsigned m=0; m<=M; m++)
    {
        double max = -1000000.0;
        DoubleVector ut = u3.row(m);
        for (unsigned int n=0; n<=N; n++)
        {
            if (max < fabs(ut[n]-V[n])) max = fabs(ut[n]-V[n]);
        }
        fprintf(file1, "%.10f ",max);
        fflush(file1);
    }

    fclose(file1);
    M=1000;
}

void ArtProblem1L2::image3Generate()
{
    DoubleVector y0;
    optimizeK = true;
    optimizeZ = true;
    optimizeE = true;

    if (optimizeK)
    {
        y0 << -0.103894  <<  -1.837793;
    }

    if (optimizeZ)
    {
        y0 << 13.063753 << 15.239693;
    }
    if (optimizeE)
    {
        y0 << 0.067496 << 0.793735;
    }

    DoubleVector k,z,e;
    getParameters(k,z,e,y0);

    withError = true;
    persent = 0.03;
    FILE *file1 = fopen("image_3_du.txt", "w");
    M *= 2;

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

//        double max = 0.5*(ut[0]-V[0])*(ut[0]-V[0]);
//        for (unsigned int n=1; n<=N-1; n++)
//        {
//            max += (ut[n]-V[n])*(ut[n]-V[n]);
//        }
//        max += 0.5*(ut[N]-V[N])*(ut[N]-V[N]);

        fprintf(file1, "%.10f ",max*hx);
        fflush(file1);
    }
    fprintf(file1, "\n");
    fclose(file1);
    M=1000;
}

void ArtProblem1L2::startOptimize()
{
    DoubleVector x0;
    if (optimizeK)
    {
        x0 << -4.5400 << -7.4500; //k
        //        x0 << -4.9700 << -1.4900; //
    }
    else
    {
        K.clear();
        K << -8.5000 << -2.7000; //k
    }

    if (optimizeZ)
    {
        x0 << +12.1000 << +14.900; //z
        //        x0 << +10.7400 << +7.0600; //
    }
    else
    {
        Z.clear();
        Z << +2.1000 << +4.9000; //z
    }

    if (optimizeE)
    {
        x0 << 0.2500 << 0.7500; //e
        //        x0 << 0.0100 << 0.0600;
    }
    else
    {
        E.clear();
        E << +0.02000 << +0.08000; //e
    }

    FILE *file = fopen("sample2.txt", "w");
    fclose(file);
    DD = 1.0;
    R = 1.0;
    optimize(x0);
    while (R < 10000.0)
    {
        R *= 10.0;
        IPrinter::printSeperatorLine();
        optimize(x0);
    }

    DoubleMatrix u;
    DoubleVector k,z,e;
    getParameters(k,z,e,x0);

    IPrinter::printSeperatorLine();
    printf("Optimal k: %20.14f %20.14f\n", k[0], k[1]);
    printf("Optimal z: %20.14f %20.14f\n", z[0], z[1]);
    printf("Optimal e: %20.14f %20.14f\n", e[0], e[1]);
    IPrinter::printSeperatorLine();

    /*
     * Image 1 v1
     */

    // initial temperatures
    //vfi << 0.0 << 1.0 << 2.0;
    // environment temperatures
    //vtt << +19.0 << +20.0 << +21.0;

    //M *= 3;

    //    {

    //        FILE *file1 = fopen("image_1_v.txt", "w");

    //        fi = +0.0;
    //        tt = +19.0;
    //        DoubleMatrix u1;
    //        calculateU(u1,k,z,e);
    //        for (int unsigned m=0; m<=M; m++)
    //        {
    //            double v = vf(m,k,z,e,u1);
    //            fprintf(file1, "%.10f ",v);
    //            fflush(file1);
    //        }
    //        fprintf(file1, "\n");

    //        //    vfi.clear(); vfi << 5.0;
    //        //    vtt.clear(); vtt << 25.0;

    //        fi = -3.0;
    //        tt = +25.0;
    //        DoubleMatrix u2;
    //        calculateU(u2,k,z,e);
    //        for (int unsigned m=0; m<=M; m++)
    //        {
    //            double v = vf(m,k,z,e,u2);
    //            fprintf(file1, "%.10f ",v);
    //            fflush(file1);
    //        }

    //        fclose(file1);
    //    }

    //    {
    //        FILE *file1 = fopen("image_2_du.txt", "w");
    //        double MM = M;

    //        fi = +0.0;
    //        tt = +19.0;
    //        for (int unsigned m=0; m<=MM; m++)
    //        {
    //            M = m;
    //            DoubleMatrix u1;
    //            calculateU(u1,k,z,e);

    //            DoubleVector uT = u1.row(u1.rows()-1);
    //            double max = -1000000.0;
    //            for (unsigned int n=0; n<=N; n++)
    //            {
    //                if (max < fabs(uT[n]-V[n])) max = fabs(uT[n]-V[n]);
    //            }

    //            fprintf(file1, "%.10f ",max);
    //            fflush(file1);
    //        }
    //        fprintf(file1, "\n");

    //    vfi.clear(); vfi << 5.0;
    //    vtt.clear(); vtt << 25.0;

    //        fi = -3.0;
    //        tt = +25.0;
    //        for (int unsigned m=0; m<=MM; m++)
    //        {
    //            M = m;
    //            DoubleMatrix u2;
    //            calculateU(u2,k,z,e);

    //            DoubleVector uT = u2.row(u2.rows()-1);
    //            double max = -1000000.0;
    //            for (unsigned int n=0; n<=N; n++)
    //            {
    //                if (max < fabs(uT[n]-V[n])) max = fabs(uT[n]-V[n]);
    //            }

    //            fprintf(file1, "%.10f ",max);
    //            fflush(file1);
    //        }

    //        fclose(file1);
    //    }

    //    {
    //        //withError = true;
    //        persent = 0.01;
    //        FILE *file1 = fopen("image_3_du_1.txt", "w");
    //        double MM = M;
    //        fi = +0.0;
    //        tt = +19.0;
    //        for (int unsigned m=0; m<=MM; m++)
    //        {
    //            M = m;
    //            DoubleMatrix u1;
    //            calculateU(u1,k,z,e);

    //            DoubleVector uT = u1.row(u1.rows()-1);
    //            double max = -1000000.0;
    //            for (unsigned int n=0; n<=N; n++)
    //            {
    //                if (max < fabs(uT[n]-V[n])) max = fabs(uT[n]-V[n]);
    //            }

    //            fprintf(file1, "%.10f ",max);
    //            fflush(file1);
    //        }
    //        fprintf(file1, "\n");
    //        fclose(file1);
    //    }

}

void ArtProblem1L2::optimize(DoubleVector &x0) const
{
    ArtProblem1L2* p = const_cast<ArtProblem1L2*>(this);
    ConjugateGradient g;
    g.setFunction(p);
    g.setGradient(p);
    g.setPrinter(p);
    g.setProjection(p);
    g.setEpsilon1(0.0001);//0.00000001
    g.setEpsilon2(0.0001);//0.00000001
    g.setEpsilon3(0.0001);//0.00000001
    g.setR1MinimizeEpsilon(10.0, 0.0001); //0.00000001
    g.setNormalize(true);
    g.showEndMessage(true);
    g.setResetIteration(false);
    g.calculate(x0);
}

void ArtProblem1L2::imager2L()
{
    DoubleVector y0;
    optimizeK = true;
    optimizeZ = true;
    optimizeE = true;

    if (optimizeK)
    {
        //y0 << -1.70338950304666 << -3.14795310656662;
        //y0 << -1.67351371107545 << -2.60464287363871;
        y0 << -1.46038450118722  <<  -3.39347208053438;
    }
    if (optimizeZ)
    {
        //y0 << +11.20714602437165 << +12.99567893148132;
        //y0 << +11.29506147253105 << +13.16809490237042;
        y0 << 11.15104030097655  <<  12.88213971188811;
    }
    if (optimizeE)
    {
        //y0 << 0.04500000000000 << 0.09500000000000;
        //y0 << 0.04431860079769 << 0.08715972328684;
        y0 << 0.06000000000000  <<   0.09500000000000;
    }


    DoubleMatrix u(101,101);
    R = 100.0;
    double max = -100000000.0;
    double min = +100000000.0;
    for (unsigned int k0=2; k0<=98; k0++)
    {
        y0[4] = 0.001*k0;
        for (unsigned int k1=2; k1<=98; k1++)
        {
            y0[5] = 0.001*k1;

            double f = 0.0;
//            if (k0 == k1 || k0 == k1-1 || k0 == k1+1)
//            {
//                f = 0.0;
//            }
//            else
//            {
//                f = fx(y0);
//            }
            f = fx(y0);
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

void ArtProblem1L2::imager3L()
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

void ArtProblem1L2::project(DoubleVector &x UNUSED_PARAM, int i UNUSED_PARAM)
{
    int p = 0;
    if (optimizeK) p+=L;

    if (optimizeZ) p+=L;

    /* e lower/upper limits */
    if (optimizeE)
    {
        if (L==2)
        {
            if (i==4) { if (x.at(4) < 0.05) x.at(4) = 0.05; if (x.at(4) > 0.95) x.at(4) = 0.95; }
            if (i==5) { if (x.at(5) < 0.05) x.at(5) = 0.05; if (x.at(5) > 0.95) x.at(5) = 0.95; }
        }

        if (L==3)
        {
            if (i == 6) { if (x.at(6) < 0.005) x.at(6) = 0.005; if (x.at(6) > 0.095) x.at(6) = 0.095; }
            if (i == 7) { if (x.at(7) < 0.005) x.at(7) = 0.005; if (x.at(7) > 0.095) x.at(7) = 0.095; }
            if (i == 8) { if (x.at(8) < 0.005) x.at(8) = 0.005; if (x.at(8) > 0.095) x.at(8) = 0.095; }
        }
    }
}

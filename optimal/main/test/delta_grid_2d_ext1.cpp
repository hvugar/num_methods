#include "delta_grid_2d_ext1.h"

void DeltaGrid2DExt1::Main(int /*argc*/, char **/*argv*/)
{
    //example1();
    example2();
}

void DeltaGrid2DExt1::example1()
{
    DeltaGrid2DExt1 ex1;

    unsigned int N = 100; double hx = 0.01;
    unsigned int M = 100; double hy = 0.01;
    DoubleMatrix u(M+1, N+1);
    for (unsigned int m=0; m<=M; m++)
    {
        double y = m*hy;
        for (unsigned int n=0; n<=N; n++)
        {
            double x = n*hx;
            u[m][n] = ex1.fx(x, y);
        }
    }
    IPrinter::printSeperatorLine();
    IPrinter::printMatrix(u);
    IPrinter::printSeperatorLine();

    DeltaGrid2D dg;
    dg.initGrid(N, hx, M, hy);
    //dg.distributeGauss(SpacePoint(0.4200, 0.6500), 1, 1);
    //dg.distributeGauss(SpacePoint(0.4224, 0.658), 1, 1);
    dg.distributeGauss(SpacePoint(0.4201, 0.6509), 1, 1);
    printf("%10.6f %10.6f\n", dg.p().x, dg.p().y);

    double z, dx, dy;
    printf("%14.10f %14.10f %14.10f\n", ex1.fx(dg.p().x, dg.p().y), ex1.dx(dg.p().x, dg.p().y), ex1.dy(dg.p().x, dg.p().y));
    puts("---");
    z = dg.consentrateInPoint(u, dx, dy, 3); printf("%14.10f %14.10f %14.10f\n", z, dx, dy);
    z = dg.consentrateInPoint(u, dx, dy, 4); printf("%14.10f %14.10f %14.10f\n", z, dx, dy);
    z = dg.consentrateInPoint(u, dx, dy, 5); printf("%14.10f %14.10f %14.10f\n", z, dx, dy);
    z = dg.lumpPointGauss(u, dx, dy);        printf("%14.10f %14.10f %14.10f\n", z, dx, dy);
    z = dg.lumpPointGauss1(u, dx, dy);       printf("%14.10f %14.10f %14.10f\n", z, dx, dy);
}

void DeltaGrid2DExt1::example2()
{
    struct A
    {
        //static double fx(double x) { return exp(x*x); }
        static double gs(SpacePoint sp, SpacePoint mu, double sigmaX, double sigmaY)
        {
            return (1.0/(2.0*M_PI*sigmaX*sigmaY))*
                    exp(-(((sp.x-mu.x)*(sp.x-mu.x))/(2.0*sigmaX*sigmaX)+((sp.y-mu.y)*(sp.y-mu.y))/(2.0*sigmaY*sigmaY)));
        }
    };

    const SpacePoint sp(0.503, 0.508);
    const SpacePoint mu(0.501, 0.509);
    const SpacePoint sg(0.01, 0.01);
    const double h = 0.01;
    DeltaGrid2D dg;
    dg.initGrid(100, h, 100, h);

    printf_s("%14.8f %14.8f %14.8f %14.8f\n", dg.onFlyWeight(sp, mu, 1), A::gs(sp, mu, h, h), dg.gaussWeight(sp, mu, 0.01, 0.01));

    double sum1 = 0.0, sum2 = 0.0, sum3 = 0.0;
    sum1 += 0.25*dg.onFlyWeight(SpacePoint(0.48, 0.49), mu, 1);
    sum1 += 0.25*dg.onFlyWeight(SpacePoint(0.48, 0.53), mu, 1);
    sum1 += 0.25*dg.onFlyWeight(SpacePoint(0.52, 0.49), mu, 1);
    sum1 += 0.25*dg.onFlyWeight(SpacePoint(0.52, 0.53), mu, 1);

    sum2 += 0.25*dg.gaussWeight(SpacePoint(0.48, 0.49), mu, 0.01, 0.01);
    sum2 += 0.25*dg.gaussWeight(SpacePoint(0.48, 0.53), mu, 0.01, 0.01);
    sum2 += 0.25*dg.gaussWeight(SpacePoint(0.52, 0.49), mu, 0.01, 0.01);
    sum2 += 0.25*dg.gaussWeight(SpacePoint(0.52, 0.53), mu, 0.01, 0.01);

    sum3 += 0.25*A::gs(SpacePoint(0.48, 0.49), mu, 0.01, 0.01);
    sum3 += 0.25*A::gs(SpacePoint(0.48, 0.53), mu, 0.01, 0.01);
    sum3 += 0.25*A::gs(SpacePoint(0.52, 0.49), mu, 0.01, 0.01);
    sum3 += 0.25*A::gs(SpacePoint(0.52, 0.53), mu, 0.01, 0.01);

    for (unsigned int i=48; i<=52; i++)
    {
        sum1 += 0.5*dg.onFlyWeight(SpacePoint(i*0.01, 0.0), mu, 1);
        sum1 += 0.5*dg.onFlyWeight(SpacePoint(i*0.01, 1.0), mu, 1);
        sum2 += 0.5*dg.gaussWeight(SpacePoint(i*0.01, 0.0), mu, 0.01, 0.01);
        sum2 += 0.5*dg.gaussWeight(SpacePoint(i*0.01, 1.0), mu, 0.01, 0.01);
        sum3 += 0.5*A::gs(SpacePoint(i*0.01, 0.0), mu, 0.01, 0.01);
        sum3 += 0.5*A::gs(SpacePoint(i*0.01, 1.0), mu, 0.01, 0.01);
    }
    for (unsigned int i=49; i<=53; i++)
    {
        sum1 += 0.5*dg.onFlyWeight(SpacePoint(0.0, i*0.01), mu, 1);
        sum1 += 0.5*dg.onFlyWeight(SpacePoint(1.0, i*0.01), mu, 1);
        sum2 += 0.5*dg.gaussWeight(SpacePoint(0.0, i*0.01), mu, 0.01, 0.01);
        sum2 += 0.5*dg.gaussWeight(SpacePoint(1.0, i*0.01), mu, 0.01, 0.01);
        sum3 += 0.5*A::gs(SpacePoint(0.0, i*0.01), mu, 0.01, 0.01);
        sum3 += 0.5*A::gs(SpacePoint(1.0, i*0.01), mu, 0.01, 0.01);
    }
    for (unsigned int i=48; i<=52; i++)
    {
        for (unsigned int j=49; j<=53; j++)
        {
            sum1 += dg.onFlyWeight(SpacePoint(i*0.01, j*0.01), mu, 1);
            sum2 += dg.gaussWeight(SpacePoint(i*0.01, j*0.01), mu, 0.01, 0.01);
            sum3 += A::gs(SpacePoint(i*0.01, j*0.01), mu, 0.01, 0.01);
        }
    }
    sum1 *= 0.0001;
    sum2 *= 0.0001;
    sum3 *= 0.0001;

    printf_s("sum1: %20.14f sum2: %20.14f sum3: %20.14f\n", sum1, sum2, sum3);
}

double DeltaGrid2DExt1::fx(double x, double y) const
{
    //return 1000000.0;
    //return x+y;
    return x*x+y*y;
    //return sin(M_PI*x)*sin(M_PI*y);
    //return sin(2.0*M_PI*x)*sin(2.0*M_PI*y);
    //return sin(8.0*M_PI*x)*sin(8.0*M_PI*y);
    //return sin(50.0*x*x)*cos(5.0*x)+sin(20.0*y*y)*sin(5.0*y)-0.5;
}

double DeltaGrid2DExt1::dx(double x, double /*y*/) const
{
    //return 0.0;
    //return 1.0;
    return 2.0*x;
    //return M_PI*cos(M_PI*x)*sin(M_PI*y);
    //return 2.0*M_PI*cos(2.0*M_PI*x)*sin(2.0*M_PI*y);
    //return 8.0*M_PI*cos(8.0*M_PI*x)*sin(8.0*M_PI*y);
    //return 100.0*x*cos(50.0*x*x)*cos(5.0*x) - 5.0*sin(50.0*x*x)*sin(5.0*x);
}

double DeltaGrid2DExt1::dy(double /*x*/, double y) const
{
    //return 0.0;
    //return 1.0;
    return 2.0*y;
    //return M_PI*sin(M_PI*x)*cos(M_PI*y);
    //return 2.0*M_PI*sin(2.0*M_PI*x)*cos(2.0*M_PI*y);
    //return 8.0*M_PI*sin(8.0*M_PI*x)*cos(8.0*M_PI*y);
    //return 40.0*y*cos(20.0*y*y)*sin(5.0*y) + 5.0*sin(20.0*y*y)*cos(5.0*y);
}

void DeltaGrid1DExt1::Main(int /*argc*/, char **/*argv*/)
{
    //example1();
    //example2();
    example3();
}

void DeltaGrid1DExt1::example1()
{
    DeltaGrid1DExt1 ex1;

    unsigned int N = 100; double hx = 0.01;
    DoubleVector u(N+1, 0.0);
    for (unsigned int n=0; n<=N; n++)
    {
        u[n] = ex1.fx(hx*n);
    }

    DeltaGrid1D dg;
    dg.initGrid(N, hx);
    dg.distributeGauss(SpacePoint(0.42), 1);
    //dg.distributeGauss(SpacePoint(0.4224, 0.6538), 1, 1);
    printf("%10.6f\n", dg.p().x);

    double z, dx;
    printf("%14.10f %14.10f\n", ex1.fx(dg.p().x), ex1.dx(dg.p().x));
    puts("---");
    z = dg.consentrateInPoint(u, dx, 3); printf("%14.10f %14.10f\n", z, dx);
    z = dg.consentrateInPoint(u, dx, 4); printf("%14.10f %14.10f\n", z, dx);
    z = dg.consentrateInPoint(u, dx, 5); printf("%14.10f %14.10f\n", z, dx);
    z = dg.lumpPointGauss(u, dx);        printf("%14.10f %14.10f\n", z, dx);
    z = dg.lumpPointGauss1(u, dx);       printf("%14.10f %14.10f\n", z, dx);
}

void DeltaGrid1DExt1::example2()
{
    struct A
    {
        //static double fx(double x) { return exp(x*x); }
        //static double gs(double x, double mu, double sigma) { return (1.0/(sqrt(2.0*M_PI)*sigma))*exp(-((x-mu)*(x-mu))/(2.0*sigma*sigma)); }
    };

    //const double x = 0.00;
    //const double c = 0.00;
    const double h = 0.01;
    const double m = 0.004;
    DeltaGrid1D dg;
    dg.initGrid(100, h);
    double node[] = {0.00, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.10};
    for (int i=0; i<=10; i++)
    {
        double w = dg.onFlyWeight(SpaceNodePDE(0, 10.0*node[i]), SpacePoint(0.5+m), 1);

        double sum = 0.0;
        sum += 0.5*dg.onFlyWeight(SpaceNodePDE(0, 0.00), SpacePoint(0.5+m), 1);
        for (int j=1; j<=99; j++)
        {
            sum += dg.onFlyWeight(SpaceNodePDE(j, j*0.01), SpacePoint(0.5+m), 1);
        }
        sum += 0.5*dg.onFlyWeight(SpaceNodePDE(100, 1.00), SpacePoint(0.5+m), 1);
        sum *= h;

        printf_s("%14.10f %14.10f %14.10f\n", 0.00+node[i], w, sum);
    }

    //    double w0 = dg.onFlyWeight(SpaceNodePDE(0,0.990), SpacePoint(m), 1);
    //    double w1 = dg.onFlyWeight(SpaceNodePDE(0,0.991), SpacePoint(m), 1);
    //    double w2 = dg.onFlyWeight(SpaceNodePDE(0,0.992), SpacePoint(m), 1);
    //    double w3 = dg.onFlyWeight(SpaceNodePDE(0,0.993), SpacePoint(m), 1);
    //    double w4 = dg.onFlyWeight(SpaceNodePDE(0,0.994), SpacePoint(m), 1);
    //    double w5 = dg.onFlyWeight(SpaceNodePDE(0,0.995), SpacePoint(m), 1);
    //    double w6 = dg.onFlyWeight(SpaceNodePDE(0,0.996), SpacePoint(m), 1);
    //    double w7 = dg.onFlyWeight(SpaceNodePDE(0,0.997), SpacePoint(m), 1);
    //    double w8 = dg.onFlyWeight(SpaceNodePDE(0,0.998), SpacePoint(m), 1);
    //    double w9 = dg.onFlyWeight(SpaceNodePDE(0,0.999), SpacePoint(m), 1);
    //    double w10= dg.onFlyWeight(SpaceNodePDE(0,1.000), SpacePoint(m), 1);

    //    printf_s("%14.10f %14.10f\n", w0, A::gs(0.200, m, 0.01));
    //    printf_s("%14.10f %14.10f\n", w1, A::gs(0.190, m, 0.01));
    //    printf_s("%14.10f %14.10f\n", w2, A::gs(0.180, m, 0.01));
    //    printf_s("%14.10f %14.10f\n", w3, A::gs(0.170, m, 0.01));
    //    printf_s("%14.10f %14.10f\n", w4, A::gs(0.160, m, 0.01));
    //    printf_s("%14.10f %14.10f\n", w5, A::gs(0.150, m, 0.01));
    //    printf_s("%14.10f %14.10f\n", w6, A::gs(0.150, m, 0.01));
    //    printf_s("%14.10f %14.10f\n", w7, A::gs(0.150, m, 0.01));
    //    printf_s("%14.10f %14.10f\n", w8, A::gs(0.150, m, 0.01));
    //    printf_s("%14.10f %14.10f\n", w9, A::gs(0.150, m, 0.01));
    //    printf_s("%14.10f %14.10f\n", w10,A::gs(0.150, m, 0.01));
}

void DeltaGrid1DExt1::example3()
{
//    double a0 = DeltaFunction::gaussian(0.41, 0.4, 0.01);
//    double a1 = DeltaFunction::gaussian(0.41, 0.4, 0.01, 1);
//    double a2 = DeltaFunction::gaussian(0.41, 0.4, 0.01, 2);
//    double a3 = DeltaFunction::gaussian(0.41, 0.4, 0.01, 3);
//    double a4 = DeltaFunction::gaussian(0.41, 0.4, 0.01, 4);
//    double a5 = DeltaFunction::gaussian(0.41, 0.4, 0.01, 5);

    DoubleMatrix m(101, 101);
    for (unsigned int j=0; j<=100; j++)
        for (unsigned int i=0; i<=100; i++)
            m[j][i] = DeltaFunction::gaussian(SpacePoint(i*0.01, j*0.01),
                                              SpacePoint(0.50, 0.50),
                                              SpacePoint(0.01, 0.01));
    double sum1 = Integral::trapezoidal(m, 0.01, 0.01);
    double sum2 = Integral::simpsons(m, 0.01, 0.01);
    printf_s("%20.14f %20.14f\n", sum1, sum2);


//    double a0 = DeltaFunction::gaussian(SpacePoint(0.51, 0.51), SpacePoint(0.50, 0.50), SpacePoint(0.01, 0.01));
//    double a1 = DeltaFunction::gaussian(SpacePoint(0.51, 0.51), SpacePoint(0.50, 0.50), SpacePoint(0.01, 0.01), 1);
//    double a2 = DeltaFunction::gaussian(SpacePoint(0.51, 0.51), SpacePoint(0.50, 0.50), SpacePoint(0.01, 0.01), 2);
//    double a3 = DeltaFunction::gaussian(SpacePoint(0.51, 0.51), SpacePoint(0.50, 0.50), SpacePoint(0.01, 0.01), 3);
//    double a4 = DeltaFunction::gaussian(SpacePoint(0.51, 0.51), SpacePoint(0.50, 0.50), SpacePoint(0.01, 0.01), 4);
//    double a5 = DeltaFunction::gaussian(SpacePoint(0.51, 0.51), SpacePoint(0.50, 0.50), SpacePoint(0.01, 0.01), 5);

//    printf_s("%20.14f %20.14f %20.14f %20.14f %20.14f %20.14f\n", a0, a1, a2, a3, a4, a5);
}

double DeltaGrid1DExt1::fx(double x) const
{
    //return 1.0;
    //return x;
    //return x*x;
    return 0.5*sin(4.0*M_PI*x*x);
    //return sin(50.0*x*x)*cos(5.0*x)-0.5;
    //return sin(20.0*x*x)*sin(5.0*x)-0.5;
}

double DeltaGrid1DExt1::dx(double x) const
{
    //return 0.0;
    //return 1.0;
    //return 2.0*x;
    return 4.0*M_PI*x*cos(4.0*M_PI*x*x);
    //return 100.0*x*cos(50.0*x*x)*cos(5.0*x) - 5.0*sin(50.0*x*x)*sin(5.0*x);
    //return 40.0*x*cos(20.0*x*x)*sin(5.0*x) + 5.0*sin(20.0*x*x)*cos(5.0*x);
}

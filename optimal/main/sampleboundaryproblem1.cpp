#include "sampleboundaryproblem1.h"

void BoundaryValueProblem1::Main(int argc UNUSED_PARAM, char **argv UNUSED_PARAM)
{
    BoundaryValueProblem1 sbp;

    {
        sbp.h = 0.1;
        sbp.N = 10;
        DoubleVector rx(sbp.N+1);
        for (unsigned int i=0; i<=sbp.N; i++) rx.at(i) = sbp.fx(i);
        IPrinter::printVector(18,14,rx);
        rx.clear();
        IPrinter::printSeperatorLine();
    }

    {
        sbp.h = 0.001;
        sbp.N = 1000;
        DoubleVector x;
        sbp.calculate2N(x, sbp.h, sbp.N);
        IPrinter::printVector(18,14,x);

//        DoubleVector rx(sbp.N+1);
//        for (unsigned int i=0; i<=sbp.N; i++) rx.at(i) = sbp.fx(i);
//        DoubleVector x1(sbp.N+1);
//        for (unsigned int i=0; i<=sbp.N; i++) x1.at(i) = rx.at(i) - x.at(i);
//        IPrinter::printVector(18,14,x1);
//        printf("Norm: %14.10f\n", x1.EuclideanNorm());

        x.clear();
//        IPrinter::printSeperatorLine();
    }
    {
        sbp.h = 0.001;
        sbp.N = 1000;
        DoubleVector x;
        sbp.calculate4NL2R(x, sbp.h, sbp.N);
        IPrinter::printVector(18,14,x);

//        DoubleVector rx(sbp.N+1);
//        for (unsigned int i=0; i<=sbp.N; i++) rx.at(i) = sbp.fx(i);
//        DoubleVector x1(sbp.N+1);
//        for (unsigned int i=0; i<=sbp.N; i++) x1.at(i) = rx.at(i) - x.at(i);
//        IPrinter::printVector(18,14,x1);
//        printf("Norm: %14.10f\n", x1.EuclideanNorm());

        x.clear();
//        IPrinter::printSeperatorLine();
    }
    /*
    {
        sbp.h = 0.001;
        sbp.N = 1000;
        DoubleVector x;
        sbp.calculate4NR2L(x, sbp.h, sbp.N);
        printf("%14.10f %14.10f %14.10f %14.10f\n", sbp.fx(1), sbp.fx(2), sbp.fx(3), sbp.fx(4));
        IPrinter::printVector(14,10,x);

        DoubleVector rx(sbp.N+1);
        for (unsigned int i=0; i<=sbp.N; i++) rx.at(i) = sbp.fx(i);
//        DoubleVector x1(sbp.N+1);
//        for (unsigned int i=0; i<=sbp.N; i++) x1.at(i) = rx.at(i) - x.at(i);
//        IPrinter::printVector(14,10,x1);
//        printf("Norm: %14.10f\n", x1.EuclideanNorm());

        x.clear();

        //        double x0 = 0*sbp.h;
        //        double x1 = 1*sbp.h;
        //        double x2 = 2*sbp.h;
        //        double x3 = 3*sbp.h;
        //        double x4 = 4*sbp.h;
        //        double x5 = 5*sbp.h;
        //        double x6 = 6*sbp.h;

        //        double a0 = -96.0*exp(-2.0*x0)*sin(10.0*x0) - 40.0*exp(-2.0*x0)*cos(10.0*x0);
        //        double a1 = -96.0*exp(-2.0*x1)*sin(10.0*x1) - 40.0*exp(-2.0*x1)*cos(10.0*x1);
        //        double a2 = -96.0*exp(-2.0*x2)*sin(10.0*x2) - 40.0*exp(-2.0*x2)*cos(10.0*x2);
        //        double a3 = -96.0*exp(-2.0*x3)*sin(10.0*x3) - 40.0*exp(-2.0*x3)*cos(10.0*x3);
        //        double a4 = -96.0*exp(-2.0*x4)*sin(10.0*x4) - 40.0*exp(-2.0*x4)*cos(10.0*x4);
        //        double a5 = -96.0*exp(-2.0*x5)*sin(10.0*x5) - 40.0*exp(-2.0*x5)*cos(10.0*x5);
        //        double a6 = -96.0*exp(-2.0*x6)*sin(10.0*x6) - 40.0*exp(-2.0*x6)*cos(10.0*x6);

        //        double b0 = (+812.0*sbp.fx(0) - 3132.0*sbp.fx(1)+ 5265.0*sbp.fx(2) - 5080.0*sbp.fx(3) + 2970.0*sbp.fx(4) - 972.0*sbp.fx(5)  + 137.0*sbp.fx(6))/(180.0*sbp.h*sbp.h);
        //        double b1 = (+137.0*sbp.fx(0) - 147.0*sbp.fx(1) - 255.0*sbp.fx(2)  + 470.0*sbp.fx(3)  - 285.0*sbp.fx(4)  + 93.0*sbp.fx(5)   - 13.0*sbp.fx(6))/(180.0*sbp.h*sbp.h);
        //        double b2 = (-13.0*sbp.fx(0)  + 228.0*sbp.fx(1) - 420.0*sbp.fx(2)  + 200.0*sbp.fx(3)  + 15.0*sbp.fx(4)   - 12.0*sbp.fx(5)   + 2.0*sbp.fx(6))/(180.0*sbp.h*sbp.h);
        //        double b3 = (+2.0*sbp.fx(0)   - 27.0*sbp.fx(1)  + 270.0*sbp.fx(2)  - 490.0*sbp.fx(3)  + 270.0*sbp.fx(4)  - 27.0*sbp.fx(5)   + 2.0*sbp.fx(6))/(180.0*sbp.h*sbp.h);
        //        double b4 = (+2.0*sbp.fx(0)   - 12.0*sbp.fx(1)  + 15.0*sbp.fx(2)   + 200.0*sbp.fx(3)  - 420.0*sbp.fx(4)  + 228.0*sbp.fx(5)  - 13.0*sbp.fx(6))/(180.0*sbp.h*sbp.h);
        //        double b5 = (-13.0*sbp.fx(0)  + 93.0*sbp.fx(1)  - 285.0*sbp.fx(2)  + 470.0*sbp.fx(3)  - 255.0*sbp.fx(4)  - 147.0*sbp.fx(5)  + 137.0*sbp.fx(6))/(180.0*sbp.h*sbp.h);
        //        double b6 = (+137.0*sbp.fx(0) - 972.0*sbp.fx(1) + 2970.0*sbp.fx(2) - 5080.0*sbp.fx(3) + 5265.0*sbp.fx(4) - 3132.0*sbp.fx(5) + 812.0*sbp.fx(6))/(180.0*sbp.h*sbp.h);

        //        printf("%14.10f %14.10f %14.10f %14.10f %14.10f %14.10f %14.10f\n", a0, a1, a2, a3, a4, a5, a6);
        //        printf("%14.10f %14.10f %14.10f %14.10f %14.10f %14.10f %14.10f\n", b0, b1, b2, b3, b4, b5, b6);

        IPrinter::printSeperatorLine();
    }
    */
    {
        sbp.h = 0.01;
        sbp.N = 100;
        DoubleVector x;
        sbp.calculate6NL2R(x, sbp.h, sbp.N);
        //printf("%14.10f %14.10f %14.10f %14.10f %14.10f %14.10f %14.10f\n", sbp.fx(sbp.N-6), sbp.fx(sbp.N-5), sbp.fx(sbp.N-4), sbp.fx(sbp.N-3), sbp.fx(sbp.N-2), sbp.fx(sbp.N-1), sbp.fx(sbp.N));
        IPrinter::printVector(18,14,x);

//        DoubleVector rx(sbp.N+1);
//        for (unsigned int i=0; i<=sbp.N; i++) rx.at(i) = sbp.fx(i);
//        DoubleVector x1(sbp.N+1);
//        for (unsigned int i=0; i<=sbp.N; i++) x1.at(i) = rx.at(i) - x.at(i);
//        IPrinter::printVector(18,14,x1);
//        printf("Norm: %18.14f\n", x1.EuclideanNorm());

        x.clear();
//        IPrinter::printSeperatorLine();
    }
    {
        sbp.h = 0.01;
        sbp.N = 100;
        DoubleVector x;
        sbp.calculate6NR2L(x, sbp.h, sbp.N);
        //printf("%14.10f %14.10f %14.10f %14.10f %14.10f %14.10f %14.10f\n", sbp.fx(sbp.N-6), sbp.fx(sbp.N-5), sbp.fx(sbp.N-4), sbp.fx(sbp.N-3), sbp.fx(sbp.N-2), sbp.fx(sbp.N-1), sbp.fx(sbp.N));
        IPrinter::printVector(18,14,x);

//        DoubleVector rx(sbp.N+1);
//        for (unsigned int i=0; i<=sbp.N; i++) rx.at(i) = sbp.fx(i);
//        DoubleVector x1(sbp.N+1);
//        for (unsigned int i=0; i<=sbp.N; i++) x1.at(i) = rx.at(i) - x.at(i);
//        IPrinter::printVector(18,14,x1);
//        printf("Norm: %18.14f\n", x1.EuclideanNorm());

        x.clear();
//        IPrinter::printSeperatorLine();
    }
}


double BoundaryValueProblem1::r(unsigned int i UNUSED_PARAM) const
{
    return +1.0;
}

double BoundaryValueProblem1::p(unsigned int i UNUSED_PARAM) const
{
#ifdef SAMPLE_1
    return +2.0;
#endif
#ifdef SAMPLE_2
    return +0.0;
#endif
#ifdef SAMPLE_3
    return +0.0;
#endif
}

double BoundaryValueProblem1::q(unsigned int i UNUSED_PARAM) const
{
#ifdef SAMPLE_1
    return +100.0;
#endif
#ifdef SAMPLE_2
    return +100.0;
#endif
#ifdef SAMPLE_3
    return +0.0;
#endif
}

double BoundaryValueProblem1::f(unsigned int i UNUSED_PARAM) const
{
    double x UNUSED_PARAM = i*h;
#ifdef SAMPLE_1
    return -20.0*exp(-2.0*x)*cos(10.0*x);
#endif
#ifdef SAMPLE_2
    return +20.0*cos(10.0*x);
#endif
#ifdef SAMPLE_3
    return +6.0*x;
#endif
}

double BoundaryValueProblem1::boundary(Boundary bound) const
{
    if (bound==Left) return fx(0);
    if (bound==Right) return fx(N);
    return 0.0;
}

double BoundaryValueProblem1::fx(unsigned int i) const
{
    double x UNUSED_PARAM = i*h;
#ifdef SAMPLE_1
    return exp(-2.0*x)*sin(10.0*x);
#endif
#ifdef SAMPLE_2
    return x*sin(10.0*x);
#endif
#ifdef SAMPLE_3
    return x*x*x;
#endif
}



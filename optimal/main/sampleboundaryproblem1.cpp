#include "sampleboundaryproblem1.h"

void SampleBoundaryProblem1::Main(int argc UNUSED_PARAM, char **argv UNUSED_PARAM)
{
    SampleBoundaryProblem1 sbp;

    {
        sbp.h = 0.1;
        sbp.N = 10;
        DoubleVector rx(sbp.N+1);
        for (unsigned int i=0; i<=sbp.N; i++) rx.at(i) = sbp.fx(i);
        IPrinter::printVector(14,10,rx);
        rx.clear();
        IPrinter::printSeperatorLine();
    }
    {
        sbp.h = 0.001;
        sbp.N = 1000;
        DoubleVector x;
        sbp.calculate2N(x, sbp.h, sbp.N);
        IPrinter::printVector(14,10,x);

        DoubleVector rx(sbp.N+1);
        for (unsigned int i=0; i<=sbp.N; i++) rx.at(i) = sbp.fx(i);
        DoubleVector x1(sbp.N+1);
        for (unsigned int i=0; i<=sbp.N; i++) x1.at(i) = rx.at(i) - x.at(i);
        IPrinter::printVector(14,10,x1);
        printf("Norm: %14.10f\n", x1.EuclideanNorm());

        x.clear();
    }

    IPrinter::printSeperatorLine();
    {
        sbp.h = 0.01;
        sbp.N = 100;
        DoubleVector x;
        sbp.calculate4NL2R(x, sbp.h, sbp.N);
        IPrinter::printVector(14,10,x);

        DoubleVector rx(sbp.N+1);
        for (unsigned int i=0; i<=sbp.N; i++) rx.at(i) = sbp.fx(i);
        DoubleVector x1(sbp.N+1);
        for (unsigned int i=0; i<=sbp.N; i++) x1.at(i) = rx.at(i) - x.at(i);
        IPrinter::printVector(14,10,x1);
        printf("Norm: %14.10f\n", x1.EuclideanNorm());

        x.clear();
    }
    IPrinter::printSeperatorLine();
    {
        sbp.h = 0.01;
        sbp.N = 100;
        DoubleVector x;
        sbp.calculate4NR2L(x, sbp.h, sbp.N);
        IPrinter::printVector(14,10,x);

        DoubleVector rx(sbp.N+1);
        for (unsigned int i=0; i<=sbp.N; i++) rx.at(i) = sbp.fx(i);
        DoubleVector x1(sbp.N+1);
        for (unsigned int i=0; i<=sbp.N; i++) x1.at(i) = rx.at(i) - x.at(i);
        IPrinter::printVector(14,10,x1);
        printf("Norm: %14.10f\n", x1.EuclideanNorm());

        x.clear();
    }
    IPrinter::printSeperatorLine();
    {
        sbp.h = 0.001;
        sbp.N = 1000;
        DoubleVector x;
        sbp.calculate6NL2R(x, sbp.h, sbp.N);
        IPrinter::printVector(14,10,x);
        printf("%14.10f %14.10f %14.10f %14.10f %14.10f %14.10f %14.10f\n",
               sbp.fx(sbp.N-6), sbp.fx(sbp.N-5), sbp.fx(sbp.N-4), sbp.fx(sbp.N-3),
               sbp.fx(sbp.N-2), sbp.fx(sbp.N-1), sbp.fx(sbp.N));

//        DoubleVector rx(sbp.N+1);
//        for (unsigned int i=0; i<=sbp.N; i++) rx.at(i) = sbp.fx(i);
//        DoubleVector x1(sbp.N+1);
//        for (unsigned int i=0; i<=sbp.N; i++) x1.at(i) = rx.at(i) - x.at(i);
//        IPrinter::printVector(14,10,x1);
//        printf("Norm: %14.10f\n", x1.EuclideanNorm());

        x.clear();
    }
    IPrinter::printSeperatorLine();
}


double SampleBoundaryProblem1::r(unsigned int i UNUSED_PARAM) const
{
    return +1.0;
}

double SampleBoundaryProblem1::p(unsigned int i UNUSED_PARAM) const
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

double SampleBoundaryProblem1::q(unsigned int i UNUSED_PARAM) const
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

double SampleBoundaryProblem1::f(unsigned int i UNUSED_PARAM) const
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

double SampleBoundaryProblem1::boundary(Boundary bound) const
{
    if (bound==Left) return fx(0);
    if (bound==Right) return fx(N);
    return 0.0;
}

double SampleBoundaryProblem1::fx(unsigned int i) const
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



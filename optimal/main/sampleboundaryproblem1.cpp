#include "sampleboundaryproblem1.h"

void BoundaryValueProblem1::Main(int argc UNUSED_PARAM, char **argv UNUSED_PARAM)
{
    BoundaryValueProblem1 sbp;

    unsigned int w = 14;
    unsigned int p = 10;
    {
        sbp.h = 0.1;
        sbp.N = 10;
        DoubleVector rx(sbp.N+1);
        for (unsigned int i=0; i<=sbp.N; i++) rx.at(i) = sbp.fx(i);
        IPrinter::printVector(w,p,rx);
        rx.clear();
        IPrinter::printSeperatorLine();
    }

    sbp.h = 0.01;
    sbp.N = 100;
    {
        //sbp.h = 0.01;
        //sbp.N = 100;
        DoubleVector x;
        sbp.calculateX(x, sbp.h, sbp.N);
        IPrinter::printVector(w,p,x);
        x.clear();
        IPrinter::printSeperatorLine();
    }
    {
        //sbp.h = 0.01;
        //sbp.N = 100;
        DoubleVector x;
        sbp.calculateN2L2RD(x, sbp.h, sbp.N);
        IPrinter::printVector(w,p,x);
        x.clear();
        //IPrinter::printSeperatorLine();
    }
    {
        //sbp.h = 0.01;
        //sbp.N = 100;
        DoubleVector x;
        sbp.calculateN2R2LD(x, sbp.h, sbp.N);
        IPrinter::printVector(w,p,x);
        x.clear();
        IPrinter::printSeperatorLine();
    }
    {
        //sbp.h = 0.01;
        //sbp.N = 100;
        DoubleVector x;
        sbp.calculateN4L2RD(x, sbp.h, sbp.N);
        IPrinter::printVector(w,p,x);
        x.clear();
        //IPrinter::printSeperatorLine();
    }
    {
        //sbp.h = 0.01;
        //sbp.N = 100;
        DoubleVector x;
        sbp.calculateN4R2LD(x, sbp.h, sbp.N);
        IPrinter::printVector(w,p,x);
        x.clear();
        IPrinter::printSeperatorLine();
    }
    {
        //sbp.h = 0.01;
        //sbp.N = 100;
        DoubleVector x;
        sbp.calculateN6L2RD(x, sbp.h, sbp.N);
        IPrinter::printVector(w,p,x);
        x.clear();
        //IPrinter::printSeperatorLine();
    }
    {
        //sbp.h = 0.01;
        //sbp.N = 100;
        DoubleVector x;
        sbp.calculateN6R2LD(x, sbp.h, sbp.N);
        IPrinter::printVector(w,p,x);
        x.clear();
        //IPrinter::printSeperatorLine();
    }
}

double BoundaryValueProblem1::r(unsigned int i UNUSED_PARAM) const
{
    return +1.0;
}

double BoundaryValueProblem1::p(unsigned int i UNUSED_PARAM) const
{
    double x UNUSED_PARAM = i*h;
#ifdef SAMPLE_1
    return +2.0;
#endif
#ifdef SAMPLE_2
    return +0.0;
#endif
#ifdef SAMPLE_3
    return +0.0;
#endif
#ifdef SAMPLE_4
    return +2.0*x;
#endif
#ifdef SAMPLE_5
    return +3.0*x;
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
#ifdef SAMPLE_4
    return -6.0;
#endif
#ifdef SAMPLE_5
    return -6.0;
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
#ifdef SAMPLE_4
    return +6.0*x;
#endif
#ifdef SAMPLE_5
    return +2.0;
#endif
}

double BoundaryValueProblem1::boundary(BoundaryType bound) const
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
#ifdef SAMPLE_4
    return x*x*x;
#endif
#ifdef SAMPLE_5
    return x*x;
#endif
}



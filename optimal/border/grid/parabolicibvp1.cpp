#include "parabolicibvp1.h"

void ParabolicIBVP1::Main(int argc UNUSED_PARAM, char *argv[] UNUSED_PARAM)
{
    ParabolicIBVP1 p;
    p.setTimeDimension(Dimension(0.001, 0.0, 1.0, 0, 1000));
    p.addSpaceDimension(Dimension(0.001, 0.0, 1.0, 0, 1000));
    p.time = p.timeDimension();
    p.dim1 = p.spaceDimension(SpaceDimension::Dim1);
    {
        DoubleMatrix u0;
        clock_t t = clock();
        p.gridMethod(u0);
        t = clock() - t;
        //IPrinter::printMatrix(14,10,u0);
        printf ("It took me %ld clicks (%f seconds).\n",t,((float)t)/CLOCKS_PER_SEC);
        //IPrinter::printSeperatorLine();
    }
    {

        DoubleMatrix u01;
        clock_t t = clock();
        p.gridMethod1(u01);
        t = clock() - t;
        //IPrinter::printMatrix(14,10,u01);
        printf ("It took me %ld clicks (%f seconds).\n",t,((float)t)/CLOCKS_PER_SEC);
        //IPrinter::printSeperatorLine();
    }
    {
        DoubleMatrix u1;
        clock_t t = clock();
        p.calculateN2L2RD(u1);
        t = clock() - t;
        //IPrinter::printMatrix(14,10,u1);
        printf ("It took me %ld clicks (%f seconds).\n",t,((float)t)/CLOCKS_PER_SEC);
        //IPrinter::printSeperatorLine();
    }
    {
        DoubleMatrix u11;
        clock_t t = clock();
        p.calculateN2L2RD1(u11);
        t = clock() - t;
        //IPrinter::printMatrix(14,10,u11);
        printf ("It took me %ld clicks (%f seconds).\n",t,((float)t)/CLOCKS_PER_SEC);
        //IPrinter::printSeperatorLine();
    }
    {
        DoubleMatrix u2;
        clock_t t = clock();
        p.calculateN4L2RD(u2);
        t = clock() - t;
        //IPrinter::printMatrix(14,10,u2);
        printf ("It took me %ld clicks (%f seconds).\n",t,((float)t)/CLOCKS_PER_SEC);
        //IPrinter::printSeperatorLine();
    }
    {
        DoubleMatrix u21;
        clock_t t = clock();
        p.calculateN4L2RD1(u21);
        t = clock() - t;
        //IPrinter::printMatrix(14,10,u21);
        printf ("It took me %ld clicks (%f seconds).\n",t,((float)t)/CLOCKS_PER_SEC);
        IPrinter::printSeperatorLine();
    }
}

ParabolicIBVP1::ParabolicIBVP1()
{
}

double ParabolicIBVP1::initial(unsigned int n) const
{
    double x UNUSED_PARAM = n*dim1.step();
#ifdef SAMPLE_1
    return x*x;
#endif
#ifdef SAMPLE_2
    return x*x;
#endif
#ifdef SAMPLE_3
    return 0.0;
#endif
#ifdef SAMPLE_4
    return 0.0;
#endif
#ifdef SAMPLE_5
    return 0.0;
#endif
    //    return U(n,mtimeDimension.minN());
}

double ParabolicIBVP1::initial(const SpaceNode& sn) const
{
    double x UNUSED_PARAM = sn.x;
#ifdef SAMPLE_1
    return x*x;
#endif
#ifdef SAMPLE_2
    return x*x;
#endif
#ifdef SAMPLE_3
    return 0.0;
#endif
#ifdef SAMPLE_4
    return 1.0;
#endif
#ifdef SAMPLE_5
    return 0.0;
#endif
    //    return U(n.i,time.minN());
}

double ParabolicIBVP1::boundary(unsigned int m, BoundaryType boundary) const
{
    double t UNUSED_PARAM = m*time.step();
#ifdef SAMPLE_1
    if (boundary == Left)  return t;
    if (boundary == Right) return t+1.0;
#endif
#ifdef SAMPLE_2
    if (boundary == Left)  return t*t;
    if (boundary == Right) return t*t+1.0;
#endif
#ifdef SAMPLE_3
    if (boundary == Left)  return 0.0;
    if (boundary == Right) return t*sin(1.0);
#endif
#ifdef SAMPLE_4
    if (boundary == Left)  return 1.0;
    if (boundary == Right) return t*sin(10.0)+exp(2.0*t);
#endif
#ifdef SAMPLE_5
    double k = 20.0;
    if (boundary == Left)  return t;
    if (boundary == Right) return exp(k*1.0)*t;
#endif
    return 0.0;
}

double ParabolicIBVP1::boundary(const SpaceNode &sn, const TimeNode &tn) const
{
    double t UNUSED_PARAM = tn.t;
    double x UNUSED_PARAM = sn.x;
#ifdef SAMPLE_1
    return x*x+t;
#endif
#ifdef SAMPLE_2
    return x*x+t*t;
#endif
#ifdef SAMPLE_3
    return t*sin(x);
#endif
#ifdef SAMPLE_4
    return t*sin(10.0*x)+exp(2*x*t);
#endif
#ifdef SAMPLE_5
    double k = 20.0;
    return exp(k*x)*t;
#endif
    //return U(sn.i,tn.i);
}

double f1(double x, double t, double a)
{
#ifdef SAMPLE_1
    return 1.0 - 2.0*a;
#endif
#ifdef SAMPLE_2
    return 2.0*t - 2.0*a;
#endif
#ifdef SAMPLE_3
    return sin(x)*(1.0+a*t);
#endif
#ifdef SAMPLE_4
    return (sin(10.0*x) + 2.0*x*exp(2.0*x*t) - a*(4.0*t*t*exp(2.0*x*t) - 100.0*t*sin(10.0*x)));
#endif
#ifdef SAMPLE_5
    double k = 20.0;
    return exp(k*x)*(1.0 - a*k*k*t);
#endif
}

double ParabolicIBVP1::f(unsigned int n, unsigned int m) const
{
    double t = m*time.step();
    double x = n*dim1.step();
    return f1(x,t,a(n,m));
}

double ParabolicIBVP1::f(const SpaceNode &sn, const TimeNode &tn) const
{
    double x = sn.x;
    double t = tn.t;
    return f1(x,t,a(sn,tn));
}

double ParabolicIBVP1::a(unsigned int n UNUSED_PARAM, unsigned int m UNUSED_PARAM) const
{
    return 1.0;
}

double ParabolicIBVP1::a(const SpaceNode &sn UNUSED_PARAM, const TimeNode &tn UNUSED_PARAM) const
{
    return 1.0;
}

//double ParabolicIBVP1::U(unsigned int n, unsigned int m) const
//{
//    double t = m*time.step();
//    double x = n*dim1.step();
//#ifdef SAMPLE_1
//    return x*x + t;
//#endif
//#ifdef SAMPLE_2
//    return x*x + t*t;
//#endif
//#ifdef SAMPLE_3
//    return sin(x)*t;
//#endif
//#ifdef SAMPLE_4
//    return t*sin(10.0*x) + exp(2.0*x*t);
//#endif
//#ifdef SAMPLE_5
//    double k = 20.0;
//    return exp(k*x)*t;
//#endif
//    return 0.0;
//}

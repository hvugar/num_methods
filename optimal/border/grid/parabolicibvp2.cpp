#include "parabolicibvp2.h"

void ParabolicIBVP2::Main(int argc UNUSED_PARAM, char *argv[] UNUSED_PARAM)
{
    ParabolicIBVP2 p;
    p.setTimeDimension(Dimension(0.001, 1000, 0));
    p.addSpaceDimension(Dimension(0.01, 100, 0));
    p.addSpaceDimension(Dimension(0.01, 100, 0));
//    {
//        Dimension time = mtimeDimension;
//        Dimension dim1 = mspaceDimension.at(0);
//        Dimension dim2 = mspaceDimension.at(1);

//        double ht = time.step();
//        double h1 = dim1.step();
//        double h2 = dim2.step();
//        unsigned int M = time.maxN();
//        unsigned int N1 = dim1.maxN();
//        unsigned int N2 = dim2.maxN();

//        DoubleMatrix u(N2+1, N1+1);
//        clock_t t = clock();
//        for (unsigned int n2=0; n2<=N2; n2++)
//        {
//            for (unsigned int n1=0; n1<=N1; n1++)
//            {
//                u[n2][n1] = U(n1,n2,M);
//            }
//        }
//        t = clock() - t;
//        IPrinter::printMatrix(14,10,u);
//        printf ("It took me %ld clicks (%f seconds).\n",t,((float)t)/CLOCKS_PER_SEC);
//        IPrinter::printSeperatorLine();
//    }
    {
        DoubleMatrix u;
        clock_t t = clock();
        p.calculateMVD(u);
        t = clock() - t;
        IPrinter::printMatrix(14,10,u);
        printf ("It took me %ld clicks (%f seconds).\n",t,((float)t)/CLOCKS_PER_SEC);
        IPrinter::printSeperatorLine();
    }
    {
        DoubleMatrix u;
        clock_t t = clock();
        p.calculateMVD1(u);
        t = clock() - t;
        IPrinter::printMatrix(14,10,u);
        printf ("It took me %ld clicks (%f seconds).\n",t,((float)t)/CLOCKS_PER_SEC);
        IPrinter::printSeperatorLine();
    }
    {
        DoubleMatrix u;
        clock_t t = clock();
        p.calculateMVD2(u);
        t = clock() - t;
        IPrinter::printMatrix(14,10,u);
        printf ("It took me %ld clicks (%f seconds).\n",t,((float)t)/CLOCKS_PER_SEC);
        IPrinter::printSeperatorLine();
    }
}

ParabolicIBVP2::ParabolicIBVP2()
{}

double ParabolicIBVP2::initial(const SpaceNode &sn) const
{
    C_UNUSED(sn);
#ifdef SAMPLE_0
    return 0.0;
#endif
#ifdef SAMPLE_1
    double x = sn.x;
    double y = sn.y;
    return x*x + y*y;
#endif
#ifdef SAMPLE_2
    double x = sn.x;
    double y = sn.y;
    return x*x + y*y;
#endif
#ifdef SAMPLE_3
    return 0.0;
#endif
    return 0.0;
}

double ParabolicIBVP2::boundary(const SpaceNode &sn, const TimeNode &tn, BoundaryType) const
{
#ifdef SAMPLE_0
    return 0.0;
#endif
#ifdef SAMPLE_1
    double t = tn.t;
    double x = sn.x;
    double y = sn.y;
    return x*x + y*y + t;
#endif
#ifdef SAMPLE_2
    double t = tn.t;
    double x = sn.x;
    double y = sn.y;
    return x*x + y*y + t*t;
#endif
#ifdef SAMPLE_3
    double t = tn.t;
    double x = sn.x;
    double y = sn.y;
    return (sin(2.0*x)+cos(3.0*y))*t;
#endif
    return 0.0;
}

double ParabolicIBVP2::f(const SpaceNode &sn, const TimeNode &tn) const
{
#ifdef SAMPLE_0
    return 0.0;
#endif
#ifdef SAMPLE_1
    double a1 = a(sn,tn);
    return 1.0 - 4.0*a1;
#endif
#ifdef SAMPLE_2
    double a1 = a(sn,tn);
    return 2.0*t - 4.0*a1;
#endif
#ifdef SAMPLE_3
    double t = tn.t;
    double x = sn.x;
    double y = sn.y;
    double a1 = a(sn,tn);
    return (sin(2.0*x)*(1.0+4.0*a1*t))+(cos(3.0*y)*(1.0+9.0*a1*t));
#endif
    return 0.0;
}

double ParabolicIBVP2::a(const SpaceNode &sn UNUSED_PARAM, const TimeNode &tn UNUSED_PARAM) const
{
    return 1.0;
}

double ParabolicIBVP2::U(unsigned int i, unsigned int j, unsigned int k) const
{
    C_UNUSED(i);
    C_UNUSED(j);
    C_UNUSED(k);
    return 0.0;
}

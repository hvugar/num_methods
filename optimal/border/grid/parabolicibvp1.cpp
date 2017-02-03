#include "parabolicibvp1.h"

void ParabolicIBVP1::Main(int argc UNUSED_PARAM, char *argv[] UNUSED_PARAM)
{
    ParabolicIBVP1 p;
    p.setTimeDimension(Dimension(0.01, 300, 200));
    p.addSpaceDimension(Dimension(0.001, 1100, 100));
    p.time = p.timeDimension();
    p.dim1 = p.spaceDimension(SpaceDimension::Dim1);
    printf("%d %d %d %d\n", p.time.minN(), p.time.maxN(), p.dim1.minN(), p.dim1.maxN());
    {
        unsigned int minN = p.dim1.minN();
        unsigned int maxN = p.dim1.maxN();
        unsigned int N = p.dim1.sizeN();

        unsigned int minM = p.time.minN();
        unsigned int maxM = p.time.maxN();

        DoubleMatrix u(p.time.sizeN()+1, N+1);

        clock_t t = clock();
        for (unsigned int j=minM; j<=maxM; j++)
        {
            for (unsigned int i=minN; i<=maxN; i++)
            {
                u[j-minM][i-minN] = p.U(i,j);
            }
        }
        t = clock() - t;
        IPrinter::printMatrix(14,10,u);
        printf ("It took me %ld clicks (%f seconds).\n",t,((float)t)/CLOCKS_PER_SEC);
        IPrinter::printSeperatorLine();
    }
    {
        DoubleMatrix u;
        clock_t t = clock();
        p.gridMethod2(u);
        t = clock() - t;
        IPrinter::printMatrix(14,10,u);
        printf ("It took me %ld clicks (%f seconds).\n",t,((float)t)/CLOCKS_PER_SEC);
        IPrinter::printSeperatorLine();
    }
    {
        DoubleMatrix u;
        clock_t t = clock();
        p.gridMethod(u);
        t = clock() - t;
        IPrinter::printMatrix(14,10,u);
        printf ("It took me %ld clicks (%f seconds).\n",t,((float)t)/CLOCKS_PER_SEC);
        IPrinter::printSeperatorLine();
    }
    {
        DoubleMatrix u;
        clock_t t = clock();
        p.calculateN2L2RD(u);
        t = clock() - t;
        IPrinter::printMatrix(14,10,u);
        printf ("It took me %ld clicks (%f seconds).\n",t,((float)t)/CLOCKS_PER_SEC);
        IPrinter::printSeperatorLine();
    }
    return;
    {
        DoubleMatrix u;
        clock_t t = clock();
        p.calculateN2L2RD1(u);
        t = clock() - t;
        IPrinter::printMatrix(14,10,u);
        printf ("It took me %ld clicks (%f seconds).\n",t,((float)t)/CLOCKS_PER_SEC);
        IPrinter::printSeperatorLine();
    }
    {
        DoubleMatrix u;
        clock_t t = clock();
        p.calculateN4L2RD(u);
        t = clock() - t;
        IPrinter::printMatrix(14,10,u);
        printf ("It took me %ld clicks (%f seconds).\n",t,((float)t)/CLOCKS_PER_SEC);
        IPrinter::printSeperatorLine();
    }
    {
        DoubleMatrix u;
        clock_t t = clock();
        p.calculateN4L2RD1(u);
        t = clock() - t;
        IPrinter::printMatrix(14,10,u);
        printf ("It took me %ld clicks (%f seconds).\n",t,((float)t)/CLOCKS_PER_SEC);
        IPrinter::printSeperatorLine();
    }
}

ParabolicIBVP1::ParabolicIBVP1()
{}

double finitial(double x)
{
#ifdef SAMPLE_1
    return x*x+2.0;
#endif
#ifdef SAMPLE_2
    return x*x;
#endif
#ifdef SAMPLE_3
    return 0.0;
#endif
#ifdef SAMPLE_4
    return 2.0*sin(10.0*x)+exp(4.0*x);
#endif
#ifdef SAMPLE_5
    return 0.0;
#endif
    return 0.0;
}

double ParabolicIBVP1::initial(unsigned int n) const
{
    double x UNUSED_PARAM = n*dim1.step();
    return finitial(x);
}

double ParabolicIBVP1::initial(const SpaceNode& sn) const
{
    double x UNUSED_PARAM = sn.x;
    return finitial(x);
}

double fboundary(double t, BoundaryValueProblem::BoundaryType boundary)
{
#ifdef SAMPLE_1
    if (boundary == BoundaryValueProblem::Left)  return 0.01+t;
    if (boundary == BoundaryValueProblem::Right) return 1.21+t;
#endif
#ifdef SAMPLE_2
    return x*x + t*t;
#endif
#ifdef SAMPLE_3
    return sin(x)*t;
#endif
#ifdef SAMPLE_4
    if (boundary == BoundaryValueProblem::Left)  return t*sin(1.0)+exp(0.2*t);
    if (boundary == BoundaryValueProblem::Right) return t*sin(11.0)+exp(2.2*t);
#endif
#ifdef SAMPLE_5
    double k = 20.0;
    return exp(k*x)*t;
#endif
    return 0.0;
}

double ParabolicIBVP1::boundary(unsigned int m, BoundaryType boundary) const
{
    double t UNUSED_PARAM = m*time.step();
    return fboundary(t,boundary);
}

double ParabolicIBVP1::boundary(const SpaceNode &sn, const TimeNode &tn) const
{
    if (sn.i==dim1.minN()) return fboundary(tn.t, BoundaryValueProblem::Left);
    if (sn.i==dim1.maxN()) return fboundary(tn.t, BoundaryValueProblem::Right);
    return 0.0;
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

double ParabolicIBVP1::U(unsigned int n, unsigned int m) const
{
    double t = m*time.step();
    double x = n*dim1.step();
#ifdef SAMPLE_1
    return x*x + t;
#endif
#ifdef SAMPLE_2
    return x*x + t*t;
#endif
#ifdef SAMPLE_3
    return sin(x)*t;
#endif
#ifdef SAMPLE_4
    return t*sin(10.0*x) + exp(2.0*x*t);
#endif
#ifdef SAMPLE_5
    double k = 20.0;
    return exp(k*x)*t;
#endif
    return 0.0;
}

#include "parabolicibvp1.h"

void ParabolicIBVP1::Main(int argc UNUSED_PARAM, char *argv[] UNUSED_PARAM)
{
    ParabolicIBVP1 p;
    p.setTimeDimension(Dimension(0.001, 1000, 0));
    p.addSpaceDimension(Dimension(0.001, 1000, 0));
    {
        unsigned int minN = p.spaceDimension(Dimension::Dim1).minN();
        unsigned int maxN = p.spaceDimension(Dimension::Dim1).maxN();
        unsigned int N = p.spaceDimension(Dimension::Dim1).sizeN();

        unsigned int minM = p.timeDimension().minN();
        unsigned int maxM = p.timeDimension().maxN();
        unsigned int M = p.timeDimension().sizeN();

        DoubleMatrix u(M+1, N+1);

        for (unsigned int j=minM; j<=maxM; j++)
        {
            for (unsigned int i=minN; i<=maxN; i++)
            {
                u[j-minM][i-minN] = p.U(i,j);
            }
        }
        IPrinter::printMatrix(14,10,u);
        IPrinter::printSeperatorLine();
    }

    {
        DoubleMatrix u;
        //p.gridMethod(u);
        //IPrinter::printMatrix(14,10,u);
        IPrinter::printSeperatorLine();
    }

    {
        DoubleMatrix u;
        //p.gridMethod1L(u);
        //IPrinter::printMatrix(14,10,u);
        IPrinter::printSeperatorLine();
    }
    {
        //DoubleMatrix u;
        //p.gridMethod1R(u);
        //IPrinter::printMatrix(14,10,u);
        //IPrinter::printSeperatorLine();
    }
    {
        //DoubleMatrix u;
        //clock_t t = clock();
        //p.gridMethod11(u);
        //t = clock() - t;
        //IPrinter::printMatrix(14,10,u);
        //printf ("It took me %ld clicks (%f seconds).\n",t,((float)t)/CLOCKS_PER_SEC);
        //IPrinter::printSeperatorLine();
    }
    {
        //DoubleVector u;
        //clock_t t = clock();
        //p.gridMethod(u);
        //t = clock() - t;
        //IPrinter::printVector(14,10,u);
        //printf ("It took me %ld clicks (%f seconds).\n",t,((float)t)/CLOCKS_PER_SEC);
        //IPrinter::printSeperatorLine();
    }
    {
        //DoubleMatrix u;
        //clock_t t = clock();
        //p.calculateN2L2RD(u);
        //t = clock() - t;
        //IPrinter::printMatrix(14,10,u);
        //printf ("It took me %ld clicks (%f seconds).\n",t,((float)t)/CLOCKS_PER_SEC);
        //IPrinter::printSeperatorLine();
    }
    {
        DoubleMatrix u;
        p.calculateN4L2RD(u);
        IPrinter::printMatrix(14,10,u);
        IPrinter::printSeperatorLine();
    }
}

double ParabolicIBVP1::initial(const SpaceNode& sn) const
{
    double x UNUSED_PARAM = sn.x;

#ifdef SAMPLE_0
    return 0.0;
#endif
#ifdef SAMPLE_13
    return 0.0;
#endif
#ifdef SAMPLE_14
    return 0.0;
#endif
#ifdef SAMPLE_15
    return 0.0;
#endif

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
    return 2.0*sin(10.0*x)+exp(4.0*x);
#endif
#ifdef SAMPLE_5
    return 0.0;
#endif
#ifdef SAMPLE_6
    return 2.0*sin(20.0*x)*cos(10.0*x);
#endif
#ifdef SAMPLE_7
    return 0.0;
#endif
#ifdef SAMPLE_8
    return 0.0;
#endif
    return NAN;
}

double ParabolicIBVP1::boundary(const SpaceNode &sn, const TimeNode &tn, BoundaryType boundary) const
{
    C_UNUSED(sn);

#ifdef SAMPLE_0
    double t = tn.t;
    if (boundary == BoundaryValueProblem::Left)  return 0.0;
    if (boundary == BoundaryValueProblem::Right) return t;
#endif

#if defined(SAMPLE_13)
    double t = tn.t;
    if (boundary == BoundaryValueProblem::Left)  return 0.0;
    if (boundary == BoundaryValueProblem::Right) return t;
#endif
#ifdef SAMPLE_14
    double t = tn.t;
    if (boundary == BoundaryValueProblem::Left)  return 0.0;
    if (boundary == BoundaryValueProblem::Right) return t;
#endif
#ifdef SAMPLE_15
    double t = tn.t;
    if (boundary == BoundaryValueProblem::Left)  return 0.0;
    if (boundary == BoundaryValueProblem::Right) return t;
#endif

#ifdef SAMPLE_1
    double t = tn.t;
    if (boundary == BoundaryValueProblem::Left)  return t;
    if (boundary == BoundaryValueProblem::Right) return 1.0+t;
#endif
#ifdef SAMPLE_2
    double x = sn.x;
    double t = tn.t;
    return x*x + t*t;
#endif
#ifdef SAMPLE_3
    double x = sn.x;
    double t = tn.t;
    return sin(x)*t;
#endif
#ifdef SAMPLE_4
    double t = tn.t;
    if (boundary == BoundaryValueProblem::Left)  return t*sin(1.0)+exp(0.2*t);
    if (boundary == BoundaryValueProblem::Right) return t*sin(11.0)+exp(2.2*t);
#endif
#ifdef SAMPLE_5
    double x = sn.x;
    double t = tn.t;
    double k = 20.0;
    return exp(k*x)*t;
#endif
#ifdef SAMPLE_6
    double t = tn.t;
    if (boundary == BoundaryValueProblem::Left)  return t;
    if (boundary == BoundaryValueProblem::Right) return 2.0*sin(20.0)*cos(10.0)+t;
#endif
#ifdef SAMPLE_7
    double t = tn.t;
    if (boundary == BoundaryValueProblem::Left)  return 0.0;
    if (boundary == BoundaryValueProblem::Right) return t;
#endif
#ifdef SAMPLE_8
    double t = tn.t;
    if (boundary == BoundaryValueProblem::Left)  return 0.0;
    if (boundary == BoundaryValueProblem::Right) return sin(4.0)*t;
#endif
    return NAN;
}

double ParabolicIBVP1::f(const SpaceNode &sn, const TimeNode &tn) const
{
    double x = sn.x;
    double t = tn.t;

#ifdef SAMPLE_0
    return x*x - 2.0*a(sn,tn)*t;
#endif
#ifdef SAMPLE_13
    return x*x*x - 6.0*a(sn,tn)*x*t;
#endif
#ifdef SAMPLE_14
    return x*x*x*x - 12.0*a(sn,tn)*x*x*t;
#endif
#ifdef SAMPLE_15
    return x*x*x*x*x - 20.0*a(sn,tn)*x*x*x*t;
#endif

#ifdef SAMPLE_1
    return 1.0 - 2.0*a(sn,tn);
#endif
#ifdef SAMPLE_2
    return 2.0*t - 2.0*a(sn,tn);
#endif
#ifdef SAMPLE_3
    return sin(x)*(1.0+a(sn,tn)*t);
#endif
#ifdef SAMPLE_4
    return (sin(10.0*x) + 2.0*x*exp(2.0*x*t) - a(sn,tn)*(4.0*t*t*exp(2.0*x*t) - 100.0*t*sin(10.0*x)));
#endif
#ifdef SAMPLE_5
    double k = 20.0;
    return exp(k*x)*(1.0 - a(sn,tn)*k*k*t);
#endif
#ifdef SAMPLE_6
    return 1 + 200.0*a(sn,tn)*(5.0*sin(20.0*x)*cos(10.0*x)+4.0*cos(20.0*x)*sin(10.0*x));
#endif
#ifdef SAMPLE_7
    return x*x*x*x - 12.0*t*x*x*a(sn,tn);
#endif
#ifdef SAMPLE_8
    return sin(4.0*x)+16.0*a(sn,tn)*sin(4.0*x);
#endif

    return NAN;
}

double ParabolicIBVP1::a(const SpaceNode &sn UNUSED_PARAM, const TimeNode &tn UNUSED_PARAM) const
{
    return 1.0;
}

double ParabolicIBVP1::U(unsigned int n, unsigned int m) const
{
    double t = m*timeDimension().step();
    double x = n*spaceDimension(Dimension::Dim1).step();

#ifdef SAMPLE_0
    return x*x*t;
#endif
#ifdef SAMPLE_13
    return x*x*x*t;
#endif
#ifdef SAMPLE_14
    return x*x*x*x*t;
#endif
#ifdef SAMPLE_15
    return x*x*x*x*x*t;
#endif

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
#ifdef SAMPLE_6
    return 2.0*sin(20.0*x)*cos(10.0*x)+t;
#endif
#ifdef SAMPLE_7
    return x*x*x*x*t;
#endif
#ifdef SAMPLE_8
    return sin(4.0*x)*t;
#endif
    return NAN;
}

void ParabolicIBVP1::layerInfo(const DoubleVector &u, unsigned int i) const
{
    if (i%100==0)
    {
        printf("%4d| ", i);
        IPrinter::printVector(14,10,u);
    }
}

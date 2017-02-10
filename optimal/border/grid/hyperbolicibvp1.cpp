#include "hyperbolicibvp1.h"

void HyperbolicIBVP1::Main(int argc UNUSED_PARAM, char *argv[] UNUSED_PARAM)
{
    HyperbolicIBVP1 p;
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
        p.gridMethod0(u);
        t = clock() - t;
        IPrinter::printMatrix(14,10,u);
        printf ("It took me %ld clicks (%f seconds).\n",t,((float)t)/CLOCKS_PER_SEC);
        IPrinter::printSeperatorLine();
    }
    {
        DoubleMatrix u;
        clock_t t = clock();
        p.gridMethod1(u);
        t = clock() - t;
        IPrinter::printMatrix(14,10,u);
        printf ("It took me %ld clicks (%f seconds).\n",t,((float)t)/CLOCKS_PER_SEC);
        IPrinter::printSeperatorLine();

        FILE* file = fopen("data.txt", "w");
        IPrinter::printMatrix(14, 10, u, u.rows(), u.cols(), NULL, file);
        fclose(file);
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
}

HyperbolicIBVP1::HyperbolicIBVP1()
{
}

double HyperbolicIBVP1::initial1(const SpaceNode &sn) const
{
    return 0.0;
}

double HyperbolicIBVP1::initial2(const SpaceNode &sn) const
{
    return sn.x*sn.x;
}

double HyperbolicIBVP1::boundary(const SpaceNode &sn, const TimeNode &tn, BoundaryType boundary) const
{
    double t = tn.t;
    if (boundary == BoundaryValueProblem::Left)  return 0.0;
    if (boundary == BoundaryValueProblem::Right) return t;
    return 0.0;
}

double HyperbolicIBVP1::f(const SpaceNode &sn, const TimeNode &tn) const
{
    return /*2.0 */- 2.0 * a(sn,tn)*tn.t;
}

double HyperbolicIBVP1::a(const SpaceNode &sn, const TimeNode &tn) const
{
    return 1.0;
}

double HyperbolicIBVP1::U(unsigned int n, unsigned int m) const
{
    double t = m*timeDimension().step();
    double x = n*spaceDimension(Dimension::Dim1).step();
    return x*x*t;
}


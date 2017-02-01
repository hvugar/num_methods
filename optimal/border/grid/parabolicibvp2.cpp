#include "parabolicibvp2.h"

void ParabolicIBVP2::Main(int argc UNUSED_PARAM, char *argv[] UNUSED_PARAM)
{
    ParabolicIBVP2 p;
    p.setTimeDimension(Dimension(0.001, 0.0, 1.0, 0, 1000));
    p.addSpaceDimension(Dimension(0.01, 0.0, 1.0, 0, 100));
    p.addSpaceDimension(Dimension(0.001, 0.0, 1.0, 0, 1000));
    p.time = p.mtimeDimension;
    p.dim1 = p.spaceDimension(SpaceDimension::Dim1);
    p.dim2 = p.spaceDimension(SpaceDimension::Dim2);
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
}

ParabolicIBVP2::ParabolicIBVP2()
{
}

double ParabolicIBVP2::initial(unsigned int i UNUSED_PARAM, unsigned int j UNUSED_PARAM) const
{
    double x UNUSED_PARAM = i*dim1.step();
    double y UNUSED_PARAM = j*dim2.step();
    return x*x + y+y;
}

double ParabolicIBVP2::initial(const SpaceNode& sn) const
{
    double x UNUSED_PARAM = sn.x;
    double y UNUSED_PARAM = sn.y;
    return x*x + y+y;
}

double ParabolicIBVP2::boundary(unsigned int i UNUSED_PARAM, unsigned int j UNUSED_PARAM, unsigned int k UNUSED_PARAM) const
{
    double t UNUSED_PARAM = 0.5*k*time.step();
    double x UNUSED_PARAM = i*dim1.step();
    double y UNUSED_PARAM = j*dim2.step();
    return x*x + y*y + t;
}

double ParabolicIBVP2::boundary(const SpaceNode &sn, const TimeNode &tn) const
{
    double t UNUSED_PARAM = 0.5*tn.t;
    double x UNUSED_PARAM = sn.x;
    double y UNUSED_PARAM = sn.y;
    return x*x + y*y + t;
}

double f1(double x UNUSED_PARAM, double y UNUSED_PARAM, double t UNUSED_PARAM, double a UNUSED_PARAM)
{
    return 1.0 - 4.0*a;
}

double ParabolicIBVP2::f(unsigned int i UNUSED_PARAM, unsigned int j UNUSED_PARAM, unsigned int k UNUSED_PARAM) const
{
    double t = 0.5*k*time.step();
    double x = i*dim1.step();
    double y = j*dim2.step();
    return f1(x,y,t,a(i,j,k));
}

double ParabolicIBVP2::f(const SpaceNode &sn, const TimeNode &tn) const
{
    double t = 0.5*tn.t;
    double x = sn.x;
    double y = sn.y;
    return f1(x,y,t,a(sn,tn));
}

double ParabolicIBVP2::a(const SpaceNode &sn UNUSED_PARAM, const TimeNode &tn UNUSED_PARAM) const
{
    return 1.0;
}

double ParabolicIBVP2::a(unsigned int i UNUSED_PARAM, unsigned int j UNUSED_PARAM, unsigned int k UNUSED_PARAM) const
{
    return 1.0;
}

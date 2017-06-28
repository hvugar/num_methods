#include "loadedheatequation.h"
#include <time.h>

void LoadedHeatEquation::Main(int argc UNUSED_PARAM, char** argv UNUSED_PARAM)
{
    LoadedHeatEquation ex1;
    ex1.lambda0 = 0.001;
    ex1.lambda1 = 1000.0;
    ex1.lambda2 = 1.0;
    ex1.a = 1.0;
    ex1.theta = 3.0;

    ex1.L = 3;
    ex1.params = new Parameter[ex1.L];
    ex1.params[0].k = -1.1; ex1.params[0].z = 10.2; ex1.params[0].xi = 25; ex1.params[0].e = 0.25;
    ex1.params[1].k = -2.5; ex1.params[1].z = 12.5; ex1.params[1].xi = 50; ex1.params[1].e = 0.50;
    ex1.params[2].k = -0.1; ex1.params[2].z = 20.5; ex1.params[2].xi = 75; ex1.params[2].e = 0.75;

    Dimension time(0.0001, 10000, 0);
    ex1.setTimeDimension(time);
    Dimension dim1(0.01, 100, 0);
    ex1.addSpaceDimension(dim1);

    DoubleVector U(dim1.sizeN()+1);
    TimeNode tn;
    tn.i = time.sizeN();
    tn.t = tn.i*time.step();
    for (unsigned int n=dim1.minN(); n<=dim1.maxN(); n++)
    {
        SpaceNode sn;
        sn.i = n;
        sn.x = n*dim1.step();
        U[n] = ex1.U(sn, tn);
    }
    IPrinter::printVector(14, 10, U);
    IPrinter::printSeperatorLine();

    DoubleVector u;
    //clock_t t = clock();
    ex1.calculateM1(u);
    //t = clock() - t;
    //printf("It took me %ld clicks (%f seconds).\n",t,((float)t)/CLOCKS_PER_SEC);
    IPrinter::printSeperatorLine();
    IPrinter::printVector(14, 10, u);
    IPrinter::printSeperatorLine();
}

double LoadedHeatEquation::U(const SpaceNode &sn, const TimeNode &tn) const
{
    double t = tn.t; C_UNUSED(t);
    double x = sn.x; C_UNUSED(x);
    return x*x*t;
}

double LoadedHeatEquation::initial(const SpaceNode &sn) const
{
    C_UNUSED(sn);
    return 2.0;
}

double LoadedHeatEquation::boundary(const SpaceNode &sn UNUSED_PARAM, const TimeNode &tn UNUSED_PARAM,
                                    BoundaryType boundary UNUSED_PARAM) const
{
    return NAN;
}

double LoadedHeatEquation::f(const SpaceNode &sn, const TimeNode &tn) const
{
    double t = tn.t; C_UNUSED(t);
    double x = sn.x; C_UNUSED(x);
    //return x*x - 2.0*a*a*t + lambda0*(x*x*t-theta);
    return 0.0;
}

double LoadedHeatEquation::g(const TimeNode &tn) const
{
    double t = tn.t; C_UNUSED(t);
    //double v = 0.0;
    //for (unsigned int s=0; s<L; s++)
    //    v += params[s].k*(params[s].e*params[s].e*t - params[s].z);
    //return 0.0 - lambda1 * (0.0 - v);
    return 0.0;
}

double LoadedHeatEquation::h(const TimeNode &tn) const
{
    double t = tn.t; C_UNUSED(t);
    //return 2.0*t + lambda2*(t-theta);
    return 0.0;
}

void LoadedHeatEquation::layerInfo(const DoubleVector &u, unsigned int m) const
{
    //IPrinter::printVector(14, 10, u);
}

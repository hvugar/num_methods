#include "loadedheatequation.h"
#include <time.h>

void LoadedHeatEquation::Main(int argc UNUSED_PARAM, char** argv UNUSED_PARAM)
{
    LoadedHeatEquation ex1;

    Dimension time(0.1, 10, 0);
    ex1.setTimeDimension(time);
    Dimension dim1(0.001, 1000, 0);
    ex1.addSpaceDimension(dim1);

    ex1.lambda0 = 0.001;
    ex1.lambda1 = 1000.0;
    ex1.lambda2 = 1.0;
    ex1.a = 1.0;
    ex1.theta = 3.0;

    ex1.L = 3;
    ex1.params.resize(ex1.L);
    ex1.params[0].k = -1.1; ex1.params[0].z = 10.2; ex1.params[0].e = 0.7415;
    ex1.params[1].k = -2.5; ex1.params[1].z = 12.5; ex1.params[1].e = 0.5845;
    ex1.params[2].k = -0.1; ex1.params[2].z = 20.5; ex1.params[2].e = 0.7145;

    IPrinter::printSeperatorLine("Real process:");

    DoubleVector U(dim1.sizeN()+1);
    TimeNodePDE tn;
    tn.i = 1;//time.sizeN();
    tn.t = tn.i*time.step();
    for (unsigned int n=dim1.minN(); n<=dim1.maxN(); n++)
    {
        SpaceNodePDE sn;
        sn.i = n;
        sn.x = n*dim1.step();
        U[n] = ex1.U(sn, tn);
    }
    IPrinter::printVector(14, 10, U);

    IPrinter::printSeperatorLine();
    DoubleVector u1;
    ex1.transferringConditions2(u1);
    IPrinter::printVector(14, 10, u1);
    IPrinter::printSeperatorLine();

    IPrinter::printSeperatorLine();
    DoubleVector u2;
    ex1.calculateM2(u2);
    IPrinter::printVector(14, 10, u2);
}

double LoadedHeatEquation::U(const SpaceNodePDE &sn, const TimeNodePDE &tn) const
{
    double t = tn.t; C_UNUSED(t);
    double x = sn.x; C_UNUSED(x);
    return x*x*t;
}

double LoadedHeatEquation::initial(const SpaceNodePDE &sn) const
{
    C_UNUSED(sn);
    return 0.0;
    //return 2.0;
}

double LoadedHeatEquation::boundary(const SpaceNodePDE &sn UNUSED_PARAM, const TimeNodePDE &tn UNUSED_PARAM,
                                    BoundaryType boundary UNUSED_PARAM) const
{
    return NAN;
}

double LoadedHeatEquation::f(const SpaceNodePDE &sn UNUSED_PARAM, const TimeNodePDE &tn UNUSED_PARAM) const
{
//    double t = tn.t; C_UNUSED(t);
//    double x = sn.x; C_UNUSED(x);
//    return x*x - 2.0*a*a*t + lambda0*(x*x*t-theta);
    return 0.0;
}

double LoadedHeatEquation::g(const TimeNodePDE &tn UNUSED_PARAM) const
{
//    double t = tn.t; C_UNUSED(t);
//    double v = 0.0;

//    Dimension dim1 = mspaceDimension.at(0);
//    double hx = dim1.step();
//    unsigned int minN = dim1.minN();
//    unsigned int maxN = dim1.maxN();
//    unsigned int N = maxN-minN;

//    for (unsigned int s=0; s<L; s++)
//    {
//        double e = params[s].e;
//        double u = e*e*t;
//        double e1 = 0.0;
//        double e2 = 0.0;
//        double u1 = 0.0;
//        double u2 = 0.0;
//        for (unsigned int n=0; n<=N-1; n++)
//        {
//            if (n*hx <= e && e < (n+1)*hx)
//            {
//                double diff = fabs(n*hx - e);

//                e1 = n*hx;
//                e2 = (n+1)*hx;

//                u1 = (1.0 - diff/hx)*e1*e1*t;
//                u2 = (diff/hx)*e2*e2*t;
//                //printf("%.10f %.10f %.10f %.10f %.10f %.10f\n", u1,u2, u, e1, e2, n*hx);
//                break;
//            }
//        }

//        v += params[s].k*( (u1 + u2) - params[s].z);
//    }
//    return 0.0 - lambda1 * (0.0 - v);
    return 0.0;
}

double LoadedHeatEquation::h(const TimeNodePDE &tn UNUSED_PARAM) const
{
//    double t = tn.t; C_UNUSED(t);
//    return 2.0*t + lambda2*(t-theta);
    return 0.0;
}

void LoadedHeatEquation::layerInfo(const DoubleVector &u, unsigned int m) const
{
    C_UNUSED(u);
    C_UNUSED(m);
    //IPrinter::printVector(14, 10, u);
}

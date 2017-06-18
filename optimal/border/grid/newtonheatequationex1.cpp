#include "newtonheatequationex1.h"

void NewtonHeatEquationEx1::Main(int argc UNUSED_PARAM, char *argv[] UNUSED_PARAM)
{
    NewtonHeatEquationEx1 ex1;
    ex1.lambda0 = 0.1;
    ex1.lambda1 = 1.0;
    ex1.lambda2 = 1.0;

    Dimension time(0.1, 10, 0);
    ex1.setTimeDimension(time);
    Dimension dim1(0.001, 1000, 0);
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

    DoubleVector u(dim1.sizeN()+1);
    ex1.calculateGM1(u,ForwardSweep);

    IPrinter::printVector(14, 10, U);
    IPrinter::printVector(14, 10, u);
}

double NewtonHeatEquationEx1::U(const SpaceNode &sn, const TimeNode &tn) const
{
    double t = tn.t; C_UNUSED(t);
    double x = sn.x; C_UNUSED(x);
#ifdef SAMPLE_1
    return x*t;
#endif
#ifdef SAMPLE_2
    return x*x*t;
#endif
    return NAN;
}

double NewtonHeatEquationEx1::initial(const SpaceNode &sn UNUSED_PARAM) const
{
#ifdef SAMPLE_1
    return 0.0;
#endif
#ifdef SAMPLE_2
    return 0.0;
#endif
    return NAN;
}

double NewtonHeatEquationEx1::boundary(const SpaceNode &sn UNUSED_PARAM, const TimeNode &tn UNUSED_PARAM, BoundaryType boundary UNUSED_PARAM) const
{
    double t = tn.t; C_UNUSED(t);
    double x = sn.x; C_UNUSED(x);
#ifdef SAMPLE_1
#endif
#ifdef SAMPLE_2
#endif
    return NAN;
}

double NewtonHeatEquationEx1::f(const SpaceNode &sn, const TimeNode &tn) const
{
    double t = tn.t; C_UNUSED(t);
    double x = sn.x; C_UNUSED(x);
#ifdef SAMPLE_1
    return x - a(sn,tn)*0.0 + lambda0*(x*t - theta0(tn));
#endif
#ifdef SAMPLE_2
    return x*x - 2.0*t*a(sn,tn) + lambda0*(x*x*t - theta0(tn));
#endif
    return NAN;
}

double NewtonHeatEquationEx1::a(const SpaceNode &sn UNUSED_PARAM, const TimeNode &tn UNUSED_PARAM) const
{
    double x = sn.x; C_UNUSED(x);
    double t = tn.t; C_UNUSED(t);
#ifdef SAMPLE_1
    return 1.0;
#endif
#ifdef SAMPLE_2
    return 1.0;
#endif
    return NAN;
}

double NewtonHeatEquationEx1::theta0(const TimeNode &tn UNUSED_PARAM) const
{
    double t = tn.t; C_UNUSED(t);
#ifdef SAMPLE_1
    return 4.0;
#endif
#ifdef SAMPLE_2
    return 4.0;
#endif
    return NAN;
}

double NewtonHeatEquationEx1::theta1(const TimeNode &tn UNUSED_PARAM) const
{
    double t = tn.t; C_UNUSED(t);
#ifdef SAMPLE_1
    return -t/lambda1;
#endif
#ifdef SAMPLE_2
    return 0.0;
#endif
    return NAN;
}

double NewtonHeatEquationEx1::theta2(const TimeNode &tn UNUSED_PARAM) const
{
    double t = tn.t; C_UNUSED(t);
#ifdef SAMPLE_1
    return t + t/lambda2;
#endif
#ifdef SAMPLE_2
    return t + (2.0*t)/lambda2;
#endif
    return NAN;
}

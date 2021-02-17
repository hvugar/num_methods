#include "heat_equation_ibvp.h"

#include "test_function.h"

#define HEAT_DIMENSION_2
//#define HEAT_HOMOGENIOUS
#define HEAT_QUADRATIC
//#define HEAT_DELTA

#if defined(HEAT_DIMENSION_1)
#define HEAT_LEFT_DIRICHLET
#define HEAT_RGHT_DIRICHLET
//#define HEAT_LEFT_ROBIN
//#define HEAT_RGHT_ROBIN
#endif

#if defined(HEAT_DIMENSION_2)
//#define HEAT_NORMAL_DIRICHLET
//#define HEAT_NORMAL_NEUMANN
#define HEAT_NORMAL_ROBIN
#endif

#if defined(HEAT_QUADRATIC)
#define HEAT_X1
#define HEAT_Y1
#define HEAT_T1
#endif

double u_fx(const IParabolicIBVP *p, const SpaceNodePDE &sn, const TimeNodePDE &tn, int dt = 0, int dx = 0, int dy = 0);
double p_fx(const IParabolicFBVP *p, const SpaceNodePDE &sn, const TimeNodePDE &tn, int dt = 0, int dx = 0, int dy = 0);

const double fa = +1.0;   // must be plus for forward
const double fb = +0.0;   // must be minus or plus for forward -  some problems on high values
const double fc = -0.0;   // must be minus for forward

const double ba = -1.0;   // must be minus for backward
const double bb = -0.0;   // must be minus or plus for forward -  some problems on high values
const double bc = +0.0;   // must be plus for backward

void HeatEquationIBVP::Main(int argc, char *argv[])
{
    C_UNUSED(argc);
    C_UNUSED(argv);
#ifdef USE_LIB_IMAGING
    QGuiApplication app(argc, argv);
#endif

    HeatEquationIBVP h;
    h.setTimeDimension(Dimension(0.010, 0, 100));
    h.setSpaceDimensionX(Dimension(0.10, 10, 20));
#ifdef HEAT_DIMENSION_2
    h.setSpaceDimensionY(Dimension(0.10, 20, 30));
#if defined(HEAT_DELTA)
    h.deltaGrid.initGrid(h.spaceDimensionX(), h.spaceDimensionY());
    h.deltaGrid.reset();
#endif
#endif

    h.setThermalDiffusivity(fa);
    h.setThermalConductivity(fb);
    h.setThermalConvection(fc);

    Benchmark bm;
    bm.tick();
#ifdef HEAT_DIMENSION_1
    h.implicit_calculate_D1V1();
    //h.explicit_calculate_D1V1();
#endif
#ifdef HEAT_DIMENSION_2
    h.implicit_calculate_D2V1();
    //h.explicit_calculate_D2V1();
#endif
    bm.tock();
    bm.printDuration();
}

void HeatEquationFBVP::Main(int argc UNUSED_PARAM, char *argv[] UNUSED_PARAM)
{
    C_UNUSED(argc);
    C_UNUSED(argv);
#ifdef USE_LIB_IMAGING
    QGuiApplication app(argc, argv);
#endif

    HeatEquationFBVP h;
    h.setTimeDimension(Dimension(0.01, 0, 100));
    h.setSpaceDimensionX(Dimension(0.010, 100, 200));
#ifdef HEAT_DIMENSION_2
    h.setSpaceDimensionY(Dimension(0.010, 200, 300));
#endif

    h.setThermalDiffusivity(ba);
    h.setThermalConductivity(bb);
    h.setThermalConvection(bc);

    Benchmark bm;
    bm.tick();
#ifdef HEAT_DIMENSION_1
    h.implicit_calculate_D1V1();
    //h.explicit_calculate_D1V1();
#endif
#ifdef HEAT_DIMENSION_2
    h.implicit_calculate_D2V1();
#endif
    bm.tock();
    bm.printDuration();
}

HeatEquationIBVP::HeatEquationIBVP() : IHeatEquationIBVP() {}

double HeatEquationIBVP::initial(const SpaceNodePDE &sn, InitialCondition) const
{
#if defined(HEAT_QUADRATIC)
    TimeNodePDE tn; tn.t = 0.0;
    //return ::u_fx(this, sn, tn);
    return TestFunction::u(tn, sn, TestFunction::FunctionValue);
#endif
#if defined(HEAT_DELTA)
    return 0.0;
#endif
#if defined(HEAT_HOMOGENIOUS)
    C_UNUSED(sn);
    return _initialTemperature;
#endif
}

double HeatEquationIBVP::boundary(const SpaceNodePDE &sn, const TimeNodePDE &tn, BoundaryConditionPDE &condition) const
{
#if defined(HEAT_DIMENSION_1)
    if (sn.i == spaceDimensionX().min())
    {
#if defined(HEAT_LEFT_DIRICHLET)
        condition = BoundaryConditionPDE(BoundaryCondition::Dirichlet); return ::u_fx(this, sn, tn);
#endif
#if defined(HEAT_LEFT_ROBIN)
        condition = BoundaryConditionPDE(BoundaryCondition::Robin, +4.0, -2.0, +1.0);
        return (condition.alpha()*::u_fx(this, sn, tn)+condition.beta()*::u_fx(this, sn, tn, -1, 1, -1))/condition.gamma();
#endif
    }
    if (sn.i == spaceDimensionX().max())
    {
#if defined(HEAT_RGHT_DIRICHLET)
        condition = BoundaryConditionPDE(BoundaryCondition::Dirichlet); return ::u_fx(this, sn, tn);
#endif
#if defined(HEAT_RGHT_ROBIN)
        condition = BoundaryConditionPDE(BoundaryCondition::Robin, +4.0, +2.0, +1.0);
        return (condition.alpha()*::u_fx(this, sn, tn)+condition.beta()*::u_fx(this, sn, tn, -1, 1, -1))/condition.gamma();
#endif
    }
#endif

#if defined(HEAT_DIMENSION_2)

#if defined(HEAT_NORMAL_DIRICHLET)
    condition = BoundaryConditionPDE(BoundaryCondition::Dirichlet); return ::u_fx(this, sn, tn);
#endif
#if defined(HEAT_NORMAL_NEUMANN)

#if defined(HEAT_QUADRATIC)
    condition = BoundaryConditionPDE(BoundaryCondition::Neumann); return (::u_fx(this, sn, tn, -1, 3, 3));
#endif

#if defined(HEAT_DELTA)
    condition = BoundaryConditionPDE(BoundaryCondition::Neumann, +0.0, +1.0, +1.0);
    return 0.0;
#endif

#endif
#if defined(HEAT_NORMAL_ROBIN)
    const double alpha = 1.0, beta = 2.0;
    condition = BoundaryConditionPDE::Robin(alpha, beta);
    //return (condition.alpha()*::u_fx(this, sn, tn)+condition.beta()*::u_fx(this, sn, tn, -1, 3, 3));
    return alpha * TestFunction::u(tn, sn, TestFunction::FunctionValue)
          + beta * TestFunction::u(tn, sn, TestFunction::SpaceNorm, spaceDimensionX(), spaceDimensionY());
#endif
#endif

#if defined(HEAT_HOMOGENIOUS)
    //condition = BoundaryConditionPDE::Robin(0.01, -1.0, 0.01);
    //condition = BoundaryConditionPDE::Neumann(1.0, 0.0); return 0.0;
    //condition = BoundaryConditionPDE(BoundaryCondition::Robin, _lambda, +1.00, _lambda); return _environmentTemperature;
    //condition = BoundaryConditionPDE::Robin(_lambda, -1.00, _lambda); return _environmentTemperature;
    condition = BoundaryConditionPDE::Neumann(); return 0.0;
#endif
}

double HeatEquationIBVP::f(const SpaceNodePDE &sn, const TimeNodePDE &tn) const
{
#if defined(HEAT_DIMENSION_1)
    const double a = thermalDiffusivity();
    const double b = thermalConductivity();
    const double c = thermalConvection();

    if (false) {
        TimeNodePDE tn0; tn0.i = tn1.i-1; tn0.t = tn0.i*0.010;
        double f1 = ::u_fx(this,sn,tn0,+1,-1,-1) - ::u_fx(this,sn,tn0,-1,+2,-1)*a - ::u_fx(this,sn,tn0,-1,+1,-1)*b - ::u_fx(this,sn,tn0,0,0,0)*c;
        double f2 = ::u_fx(this,sn,tn1,+1,-1,-1) - ::u_fx(this,sn,tn1,-1,+2,-1)*a - ::u_fx(this,sn,tn1,-1,+1,-1)*b - ::u_fx(this,sn,tn1,0,0,0)*c;
        return (f1+f2)/2.0;
    }
    else {
        return ::u_fx(this,sn,tn1,+1,-1,-1) - ::u_fx(this,sn,tn1,-1,+2,-1)*a - ::u_fx(this,sn,tn1,-1,+1,-1)*b - ::u_fx(this,sn,tn1,0,0,0)*c;
    }
#endif

#if defined(HEAT_DIMENSION_2)
#if defined(HEAT_QUADRATIC)
    const double a1 = thermalDiffusivity();
    const double a2 = thermalDiffusivity();
    const double b1 = thermalConductivity();
    const double b2 = thermalConductivity();
    const double c  = thermalConvection();
    return TestFunction::u(tn, sn, TestFunction::TimeFirstDerivative)
            - TestFunction::u(tn, sn, TestFunction::SpaceSecondDerivativeX) * a1
            - TestFunction::u(tn, sn, TestFunction::SpaceSecondDerivativeY) * a2
            - TestFunction::u(tn, sn, TestFunction::SpaceFirstDerivativeX) * b1
            - TestFunction::u(tn, sn, TestFunction::SpaceFirstDerivativeY) * b2
            - TestFunction::u(tn, sn, TestFunction::FunctionValue) * c;
#endif
#if defined(HEAT_DELTA)
    double q = 0.1;
    double val = q*deltaGrid.gaussWeight(sn, SpacePoint(1.5, 2.5), 0.01, 0.01);
    //printf("%6d %6d %14.8f\n", sn.i, sn.j, val);
    return val;
#endif
#endif

#if defined(HEAT_HOMOGENIOUS)
    return -thermalConvection() * _environmentTemperature;
#endif
}

void HeatEquationIBVP::layerInfo(const DoubleVector& u, const TimeNodePDE& tn) const
{
    C_UNUSED(u);
    C_UNUSED(tn);

    IPrinter::printVector(17, 8, u);
    //    unsigned int time_size = timeDimension().size();

    //    if (tn.i % ((time_size-1)/10) == 0) IPrinter::printVector(16, 8, u); return;

    //IPrinter::printVector(u);
    //    if (tn.i==0 || tn.i==1 || tn.i==timeDimension().max()-1 || tn.i==timeDimension().max()) IPrinter::printVector(u);
    return;

    //    int min = spaceDimensionX().min();
    //    int max = spaceDimensionX().max();

    //    if (tn.i==0 || tn.i==1 || tn.i==99 || tn.i==100) IPrinter::printVector(u);

    //    if (tn.i==1001)
    //    {
    //        double L2Norm = 0.0;
    //        double EuNorm = 0.0;
    //        double L1Norm = 0.0;
    //        TimeNodePDE tn; tn.t = 1.0;
    //        SpaceNodePDE sn;
    //        unsigned int n = 0;
    //        for (int i=min; i<=max; i++, n++)
    //        {
    //            sn.i = static_cast<int>(i);
    //            sn.x = sn.i*0.01;

    //            double k = 1.0; if (i==min || i== max) k = 0.5;
    //            L2Norm += 0.01*k*(u[n]-::u_fx(this, sn, tn))*(u[n]-::u_fx(this, sn, tn));

    //            EuNorm += (u[n]-u_fx(this, sn, tn))*(u[n]-u_fx(this, sn, tn));

    //            if (L1Norm < fabs(u[n]-::u_fx(this, sn, tn))) L1Norm = fabs(u[n]-::u_fx(this, sn, tn));
    //        }
    //        printf("L2Norm: %.10f EuNorm: %.10f L1Norm: %.10f\n", sqrt(L2Norm), sqrt(EuNorm), L1Norm);
    //    }
}

void HeatEquationIBVP::layerInfo(const DoubleMatrix &u, const TimeNodePDE &tn) const
{
    C_UNUSED(u);
    C_UNUSED(tn);

    unsigned int delimet = (timeDimension().size()-1) / 10;
    if (tn.i % delimet == 0 || tn.i == 0 || tn.i == 1)
    {
        std::string msg = std::string("time: ") + std::to_string(tn.t) +
                std::string(", min: ") + std::to_string(u.min()) +
                std::string(", max: ") + std::to_string(u.max());
        IPrinter::printSeperatorLine(msg.data());
        IPrinter::printMatrix(16, 8, u, 10, 10, msg.data());
    }
    return;

    //    if (tn.i==0 || tn.i==1 || tn.i==2 || tn.i==198 || tn.i==199 || tn.i==200)
    //    {
    //        IPrinter::printMatrix(u);
    //        IPrinter::printSeperatorLine();
    //        if (tn.i%2==0) IPrinter::printSeperatorLine();
    //    }
    //    return;

    //    if (tn.i==40000)
    //    {
    //        double norm = 0.0;
    //        double max = 0.0;
    //        TimeNodePDE tn; tn.t = 1.0;
    //        SpaceNodePDE sn;
    //        for (unsigned int j=0; j<=100; j++)
    //        {
    //            for (unsigned int i=0; i<=100; i++)
    //            {
    //                sn.x = i*0.01;
    //                sn.y = j*0.01;
    //                double k1 = 1.0; if (j==0 || j== 100) k1 = 0.5;
    //                double k2 = 1.0; if (i==0 || i== 100) k2 = 0.5;
    //                norm += 0.01*0.01*k1*k2*(u[j][i]-::u_fx(this,sn, tn))*(u[j][i]-::u_fx(this,sn, tn));

    //                if (max < fabs(u[j][i]-::u_fx(this, sn, tn))) max = fabs(u[j][i]-::u_fx(this, sn, tn));
    //            }
    //        }
    //        printf("norm: %.10f max: %.10f\n", sqrt(norm), max);
    //    }
}

//---------------------------------------------------------------------------------------------//

HeatEquationFBVP::HeatEquationFBVP() : IHeatEquationFBVP()
{}

double HeatEquationFBVP::final(const SpaceNodePDE &sn, FinalCondition) const
{
    TimeNodePDE tn; tn.t = timeDimension().max()*timeDimension().step();
    return ::p_fx(this, sn, tn);
}

double HeatEquationFBVP::boundary(const SpaceNodePDE &sn, const TimeNodePDE &tn, BoundaryConditionPDE &condition) const
{
#if defined(HEAT_DIMENSION_1)
    if (sn.i == spaceDimensionX().min())
    {
#if defined(HEAT_LEFT_DIRICHLET)
        condition = BoundaryConditionPDE(BoundaryCondition::Dirichlet); return ::p_fx(this, sn, tn);
#endif
#if defined(HEAT_LEFT_ROBIN)
        condition = BoundaryConditionPDE(BoundaryCondition::Robin, +4.0, -2.0, +1.0);
        return (condition.alpha()*::p_fx(this, sn, tn)+condition.beta()*::p_fx(this, sn, tn, -1, 1, -1))/condition.gamma();
#endif
    }
    if (sn.i == spaceDimensionX().max())
    {
#if defined(HEAT_RGHT_DIRICHLET)
        condition = BoundaryConditionPDE(BoundaryCondition::Dirichlet); return ::p_fx(this, sn, tn);
#endif
#if defined(HEAT_RGHT_ROBIN)
        condition = BoundaryConditionPDE(BoundaryCondition::Robin, +4.0, +2.0, +1.0);
        return (condition.alpha()*::p_fx(this, sn, tn)+condition.beta()*::p_fx(this, sn, tn, -1, 1, -1))/condition.gamma();
#endif
    }

#if defined(HEAT_RGHT_DIRICHLET)
    condition = BoundaryConditionPDE(BoundaryCondition::Dirichlet); return ::p_fx(this, sn, tn);
#endif
#endif

#if defined(HEAT_DIMENSION_2)
#if defined(HEAT_NORMAL_DIRICHLET)
    condition = BoundaryConditionPDE(BoundaryCondition::Dirichlet);
    return ::p_fx(this, sn, tn);
#endif
#if defined(HEAT_NORMAL_NEUMANN)
    condition = BoundaryConditionPDE(BoundaryCondition::Neumann);
    return (::p_fx(this, sn, tn)+::p_fx(this, sn, tn, -1, 3, 3));
#endif
#if defined(HEAT_NORMAL_ROBIN)
    condition = BoundaryConditionPDE::Robin(+4.0, +2.0);
    return (condition.alpha()*::p_fx(this, sn, tn)+condition.beta()*::p_fx(this, sn, tn, -1, 3, 3));
#endif
#endif
}

double HeatEquationFBVP::f(const SpaceNodePDE &sn, const TimeNodePDE &tn) const
{
#if defined(HEAT_DIMENSION_1)
    const double a = thermalDiffusivity();
    const double b = thermalConductivity();
    const double c = thermalConvection();
    return ::p_fx(this,sn,tn,+1,-1,-1)-::p_fx(this,sn,tn,-1,+2,-1)*a-::p_fx(this,sn,tn,-1,+1,-1)*b-::p_fx(this,sn,tn,0,0,0)*c;
#endif

#if defined(HEAT_DIMENSION_2)
    const double a1 = thermalDiffusivity();
    const double a2 = thermalDiffusivity();
    const double b1 = thermalConductivity();
    const double b2 = thermalConductivity();
    const double c  = thermalConvection();
    return ::p_fx(this, sn,tn,+1,-1,-1)
            - ::p_fx(this, sn,tn,-1,+2,-1)*a1 - ::p_fx(this, sn,tn,-1,-1,+2)*a2
            - ::p_fx(this, sn,tn,-1,+1,-1)*b1 - ::p_fx(this, sn,tn,-1,-1,+1)*b2
            - ::p_fx(this, sn,tn,+0,+0,+0)*c;
#endif
}

void HeatEquationFBVP::layerInfo(const DoubleVector& u, const TimeNodePDE& tn) const
{
    C_UNUSED(u);
    C_UNUSED(tn);

    if (tn.i % (timeDimension().size() / 10) == 0) IPrinter::printVector(16, 8, u); return;

    //IPrinter::printVector(u);
    if (tn.i==0 || tn.i==1 || tn.i==timeDimension().max()-1 || tn.i==timeDimension().max()) IPrinter::printVector(u);
}

void HeatEquationFBVP::layerInfo(const DoubleMatrix &u, const TimeNodePDE &tn) const
{
    C_UNUSED(u);
    C_UNUSED(tn);

    if (tn.i % (timeDimension().size() / 5) == 0)
    {
        IPrinter::printMatrix(16, 8, u);
        IPrinter::printSeperatorLine();
    }
    return;

    if (tn.i==0 || tn.i==1 || tn.i==2 || tn.i==198 || tn.i==199 || tn.i==200)
    {
        IPrinter::printMatrix(u);
        IPrinter::printSeperatorLine();
        if (tn.i%2==0) IPrinter::printSeperatorLine();
    }
}

//---------------------------------------------------------------------------------------------//

double u_fx(const IParabolicIBVP *p, const SpaceNodePDE &sn, const TimeNodePDE &tn, int dt, int dx, int dy)
{
    double res = 0.0;

#if defined(HEAT_QUADRATIC)

#if defined(HEAT_X1)
    if (dx == 0) res += sn.x;
    if (dx == 1) res += 1.0;
    if (dx == 2) res += 0.0;
    if (dx == 3)
    {
        const int xmin = p->spaceDimensionX().min();
        const int xmax = p->spaceDimensionX().max();
        const int ymin = p->spaceDimensionY().min();
        const int ymax = p->spaceDimensionY().max();

        //        if (tn.i%2==1)
        //        {
        //            if (sn.i == xmin && sn.j != ymin && sn.j != ymax) res += -1.0;
        //            if (sn.i == xmax && sn.j != ymin && sn.j != ymax) res += +1.0;
        //        }
        //        else
        //        {
        //            if (sn.i == xmin) res += -1.0;
        //            if (sn.i == xmax) res += +1.0;
        //        }

        if (sn.i == xmin && sn.j != ymin && sn.j != ymax) res += -1.0;
        if (sn.i == xmax && sn.j != ymin && sn.j != ymax) res += +1.0;

        //if (sn.i == xmin && (sn.j == ymin || sn.j == ymax)) return -1.0;
        //if (sn.i == xmax && (sn.j == ymin || sn.j == ymax)) return +1.0;
    }
#endif
#if defined(HEAT_X2)
    if (dx == 0) res += sn.x*sn.x;
    if (dx == 1) res += 2.0*sn.x;
    if (dx == 2) res += 2.0;
    if (dx == 3)
    {
        const int xmin = p->spaceDimensionX().min();
        const int xmax = p->spaceDimensionX().max();
        const int ymin = p->spaceDimensionY().min();
        const int ymax = p->spaceDimensionY().max();

        if (tn.i%2==1)
        {
            if (sn.i == xmin && sn.j != ymin && sn.j != ymax) res += -2.0*sn.x;
            if (sn.i == xmax && sn.j != ymin && sn.j != ymax) res += +2.0*sn.x;
        }
        else
        {
            if (sn.i == xmin) res += -2.0*sn.x;
            if (sn.i == xmax) res += +2.0*sn.x;
        }
    }
#endif

#if defined(HEAT_X3)
    if (dx == 0) res += sn.x*sn.x*sn.x;
    if (dx == 1) res += 3.0*sn.x*sn.x;
    if (dx == 2) res += 6.0*sn.x;
    if (dx == 3)
    {
        const int xmin = p->spaceDimensionX().min();
        const int xmax = p->spaceDimensionX().max();
        const int ymin = p->spaceDimensionY().min();
        const int ymax = p->spaceDimensionY().max();

        if (tn.i%2==1)
        {
            if (sn.i == xmin && sn.j != ymin && sn.j != ymax) res += -3.0*sn.x*sn.x;
            if (sn.i == xmax && sn.j != ymin && sn.j != ymax) res += +3.0*sn.x*sn.x;
        }
        else
        {
            if (sn.i == xmin) res += -3.0*sn.x*sn.x;
            if (sn.i == xmax) res += +3.0*sn.x*sn.x;
        }
    }
#endif

#if defined(HEAT_Y1)
    if (dy == 0) res += sn.y;
    if (dy == 1) res += 1.0;
    if (dy == 2) res += 0.0;
    if (dy == 3)
    {
        const int xmin = p->spaceDimensionX().min();
        const int xmax = p->spaceDimensionX().max();
        const int ymin = p->spaceDimensionY().min();
        const int ymax = p->spaceDimensionY().max();

        //        if (tn.i%2==0)
        //        {
        //            if (sn.j == ymin && sn.i != xmin && sn.i != xmax) res += -1.0;
        //            if (sn.j == ymax && sn.i != xmin && sn.i != xmax) res += +1.0;
        //        }
        //        else
        //        {
        //            if (sn.j == ymin) res += -1.0;
        //            if (sn.j == ymax) res += +1.0;
        //        }

        if (sn.j == ymin && sn.i != xmin && sn.i != xmax) res += -1.0;
        if (sn.j == ymax && sn.i != xmin && sn.i != xmax) res += +1.0;

        if (sn.j == ymin && (sn.i == xmin || sn.i == xmax)) return -1.0;
        if (sn.j == ymax && (sn.i == xmin || sn.i == xmax)) return +1.0;
    }
#endif
#if defined(HEAT_Y2)
    if (dy == 0) res += sn.y*sn.y;
    if (dy == 1) res += 2.0*sn.y;
    if (dy == 2) res += 2.0;
    if (dy == 3)
    {
        const int xmin = p->spaceDimensionX().min();
        const int xmax = p->spaceDimensionX().max();
        const int ymin = p->spaceDimensionY().min();
        const int ymax = p->spaceDimensionY().max();

        if (tn.i%2==0)
        {
            if (sn.j == ymin && sn.i != xmin && sn.i != xmax) res += -2.0*sn.y;
            if (sn.j == ymax && sn.i != xmin && sn.i != xmax) res += +2.0*sn.y;
        }
        else
        {
            if (sn.j == ymin) res += -2.0*sn.y;
            if (sn.j == ymax) res += +2.0*sn.y;
        }
    }
#endif
#if defined(HEAT_Y3)
    if (dy == 0) res += sn.y*sn.y*sn.y;
    if (dy == 1) res += 3.0*sn.y*sn.y;
    if (dy == 2) res += 6.0*sn.y;
    if (dy == 3)
    {
        const int xmin = p->spaceDimensionX().min();
        const int xmax = p->spaceDimensionX().max();
        const int ymin = p->spaceDimensionY().min();
        const int ymax = p->spaceDimensionY().max();

        if (tn.i%2==0)
        {
            if (sn.j == ymin && sn.i != xmin && sn.i != xmax) res += -3.0*sn.y*sn.y;
            if (sn.j == ymax && sn.i != xmin && sn.i != xmax) res += +3.0*sn.y*sn.y;
        }
        else
        {
            if (sn.j == ymin) res += -3.0*sn.y*sn.y;
            if (sn.j == ymax) res += +3.0*sn.y*sn.y;
        }
    }
#endif
#if defined(HEAT_T1)
    if (dt == 0) res += tn.t;
    if (dt == 1) res += 1.0;
    if (dt == 2) res += 0.0;
#endif
#if defined(HEAT_T2)
    if (dt == 0) res += tn.t*tn.t;
    if (dt == 1) res += 2.0*tn.t;
    if (dt == 2) res += 2.0;
#endif
#if defined(HEAT_T3)
    if (dt == 0) res += tn.t*tn.t*tn.t;
    if (dt == 1) res += 3.0*tn.t*tn.t;
    if (dt == 2) res += 6.0*tn.t;
#endif
#endif

    return res;
}

double p_fx(const IParabolicFBVP *p, const SpaceNodePDE &sn, const TimeNodePDE &tn, int dt, int dx, int dy)
{
    double res = 0.0;

#if defined(HEAT_QUADRATIC)

#if defined(HEAT_X1)
    if (dx == 0) res += sn.x;
    if (dx == 1) res += 1.0;
    if (dx == 2) res += 0.0;
    if (dx == 3)
    {
        const int xmin = p->spaceDimensionX().min();
        const int xmax = p->spaceDimensionX().max();
        const int ymin = p->spaceDimensionY().min();
        const int ymax = p->spaceDimensionY().max();

        if (tn.i%2==1)
        {
            if (sn.i == xmin && sn.j != ymin && sn.j != ymax) res += -1.0;
            if (sn.i == xmax && sn.j != ymin && sn.j != ymax) res += +1.0;
        }
        else
        {
            if (sn.i == xmin) res += -1.0;
            if (sn.i == xmax) res += +1.0;
        }
    }
#endif

#if defined(HEAT_X2)
    if (dx == 0) res += sn.x*sn.x;
    if (dx == 1) res += 2.0*sn.x;
    if (dx == 2) res += 2.0;
    if (dx == 3)
    {
        const int xmin = p->spaceDimensionX().min();
        const int xmax = p->spaceDimensionX().max();
        const int ymin = p->spaceDimensionY().min();
        const int ymax = p->spaceDimensionY().max();

        if (tn.i%2==1)
        {
            if (sn.i == xmin && sn.j != ymin && sn.j != ymax) res += -2.0*sn.x;
            if (sn.i == xmax && sn.j != ymin && sn.j != ymax) res += +2.0*sn.x;
        }
        else
        {
            if (sn.i == xmin) res += -2.0*sn.x;
            if (sn.i == xmax) res += +2.0*sn.x;
        }
    }
#endif

#if defined(HEAT_X3)
    if (dx == 0) res += sn.x*sn.x*sn.x;
    if (dx == 1) res += 3.0*sn.x*sn.x;
    if (dx == 2) res += 6.0*sn.x;
    if (dx == 3)
    {
        const int xmin = p->spaceDimensionX().min();
        const int xmax = p->spaceDimensionX().max();
        const int ymin = p->spaceDimensionY().min();
        const int ymax = p->spaceDimensionY().max();

        if (tn.i%2==1)
        {
            if (sn.i == xmin && sn.j != ymin && sn.j != ymax) res += -3.0*sn.x*sn.x;
            if (sn.i == xmax && sn.j != ymin && sn.j != ymax) res += +3.0*sn.x*sn.x;
        }
        else
        {
            if (sn.i == xmin) res += -3.0*sn.x*sn.x;
            if (sn.i == xmax) res += +3.0*sn.x*sn.x;
        }
    }
#endif

#if defined(HEAT_Y1)
    if (dy == 0) res += sn.y;
    if (dy == 1) res += 1.0;
    if (dy == 2) res += 0.0;
    if (dy == 3)
    {
        const int xmin = p->spaceDimensionX().min();
        const int xmax = p->spaceDimensionX().max();
        const int ymin = p->spaceDimensionY().min();
        const int ymax = p->spaceDimensionY().max();

        if (tn.i%2==0)
        {
            if (sn.j == ymin && sn.i != xmin && sn.i != xmax) res += -1.0;
            if (sn.j == ymax && sn.i != xmin && sn.i != xmax) res += +1.0;
        }
        else
        {
            if (sn.j == ymin) res += -1.0;
            if (sn.j == ymax) res += +1.0;
        }
    }
#endif
#if defined(HEAT_Y2)
    if (dy == 0) res += sn.y*sn.y;
    if (dy == 1) res += 2.0*sn.y;
    if (dy == 2) res += 2.0;
    if (dy == 3)
    {
        const int xmin = p->spaceDimensionX().min();
        const int xmax = p->spaceDimensionX().max();
        const int ymin = p->spaceDimensionY().min();
        const int ymax = p->spaceDimensionY().max();

        if (tn.i%2==0)
        {
            if (sn.j == ymin && sn.i != xmin && sn.i != xmax) res += -2.0*sn.y;
            if (sn.j == ymax && sn.i != xmin && sn.i != xmax) res += +2.0*sn.y;
        }
        else
        {
            if (sn.j == ymin) res += -2.0*sn.y;
            if (sn.j == ymax) res += +2.0*sn.y;
        }
    }
#endif
#if defined(HEAT_Y3)
    if (dy == 0) res += sn.y*sn.y*sn.y;
    if (dy == 1) res += 3.0*sn.y*sn.y;
    if (dy == 2) res += 6.0*sn.y;
    if (dy == 3)
    {
        const int xmin = p->spaceDimensionX().min();
        const int xmax = p->spaceDimensionX().max();
        const int ymin = p->spaceDimensionY().min();
        const int ymax = p->spaceDimensionY().max();

        if (tn.i%2==0)
        {
            if (sn.j == ymin && sn.i != xmin && sn.i != xmax) res += -3.0*sn.y*sn.y;
            if (sn.j == ymax && sn.i != xmin && sn.i != xmax) res += +3.0*sn.y*sn.y;
        }
        else
        {
            if (sn.j == ymin) res += -3.0*sn.y*sn.y;
            if (sn.j == ymax) res += +3.0*sn.y*sn.y;
        }
    }
#endif
#if defined(HEAT_T1)
    if (dt == 0) res += tn.t;
    if (dt == 1) res += 1.0;
    if (dt == 2) res += 0.0;
#endif
#if defined(HEAT_T2)
    if (dt == 0) res += tn.t*tn.t;
    if (dt == 1) res += 2.0*tn.t;
    if (dt == 2) res += 2.0;
#endif
#if defined(HEAT_T3)
    if (dt == 0) res += tn.t*tn.t*tn.t;
    if (dt == 1) res += 3.0*tn.t*tn.t;
    if (dt == 2) res += 6.0*tn.t;
#endif
#endif

    return res;
}

/***********************************************************************************************************************************/

void LoadedHeatEquationIBVP::Main(int /*argc*/, char */*argv*/[])
{
    LoadedHeatEquationIBVP lheIBVP;
    lheIBVP.setThermalDiffusivity(1.0);
    lheIBVP.setThermalConductivity(0.0);
    lheIBVP.setThermalConvection(0.0);

    std::vector<LoadedSpacePoint> loadedPoints;
    loadedPoints.push_back(LoadedSpacePoint(0.20, 0.20, 0.00, 1.0));
    loadedPoints.push_back(LoadedSpacePoint(0.80, 0.80, 0.00, 1.0));
    lheIBVP.setLoadedPoints(loadedPoints);

    lheIBVP.implicit_calculate_D2V1();
}

double LoadedHeatEquationIBVP::initial(const SpaceNodePDE &sn, InitialCondition) const
{
    //return sn.x*sn.x + sn.y*sn.y;
    return sn.x + sn.y;
}

double LoadedHeatEquationIBVP::boundary(const SpaceNodePDE &sn, const TimeNodePDE &tn, BoundaryConditionPDE &bc) const
{
    //bc = BoundaryConditionPDE::Dirichlet(1.0, 1.0); return sn.x*sn.x + sn.y*sn.y + tn.t*tn.t;
    bc = BoundaryConditionPDE::Dirichlet(); return sn.x + sn.y + tn.t;
}

double LoadedHeatEquationIBVP::f(const SpaceNodePDE &sn, const TimeNodePDE &tn) const
{
    const double a =  thermalDiffusivity();
    const double b =  thermalConductivity();
    const double c =  thermalConvection();

    //return 2.0*tn.t - 4.0*a - 2.0*b*(sn.x+sn.y) - c*(sn.x*sn.x+sn.y*sn.y+tn.t*tn.t)
    //        - 0.1*(0.2*0.2 + 0.2*0.2 + tn.t*tn.t)
    //        - 0.1*(0.8*0.8 + 0.8*0.8 + tn.t*tn.t);

    return 1.0 - 2.0*b - c*(sn.x+sn.y+tn.t)
            - 1.0*(0.2 + 0.2 + tn.t)
            - 1.0*(0.8 + 0.8 + tn.t);
}

void LoadedHeatEquationIBVP::layerInfo(const DoubleVector&, const TimeNodePDE&) const
{}

void LoadedHeatEquationIBVP::layerInfo(const DoubleMatrix& u, const TimeNodePDE&) const
{
    IPrinter::printSeperatorLine();
    IPrinter::printMatrix(u);
}

/***********************************************************************************************************************************/

void LoadedHeatEquationFBVP::Main(int /*argc*/, char */*argv*/[])
{
    LoadedHeatEquationFBVP lheIBVP;
    lheIBVP.setThermalDiffusivity(-1.0);
    lheIBVP.setThermalConductivity(0.0);
    lheIBVP.setThermalConvection(0.0);

    std::vector<LoadedSpacePoint> loadedPoints;
    loadedPoints.push_back(LoadedSpacePoint(0.20, 0.20, 0.00, +0.1));
    loadedPoints.push_back(LoadedSpacePoint(0.80, 0.80, 0.00, +0.1));
    lheIBVP.setLoadedPoints(loadedPoints);

    lheIBVP.implicit_calculate_D2V1();
}

double LoadedHeatEquationFBVP::final(const SpaceNodePDE &sn, FinalCondition) const
{
    //return sn.x*sn.x + sn.y*sn.y;
    return sn.x + sn.y + 1.0;
}

double LoadedHeatEquationFBVP::boundary(const SpaceNodePDE &sn, const TimeNodePDE &tn, BoundaryConditionPDE &bc) const
{
    //bc = BoundaryConditionPDE::Dirichlet(1.0, 1.0); return sn.x*sn.x + sn.y*sn.y + tn.t*tn.t;
    bc = BoundaryConditionPDE::Dirichlet(); return sn.x + sn.y + tn.t;
}

double LoadedHeatEquationFBVP::f(const SpaceNodePDE &sn, const TimeNodePDE &tn) const
{
    const double a =  thermalDiffusivity();
    const double b =  thermalConductivity();
    const double c =  thermalConvection();

    //return 2.0*tn.t - 4.0*a - 2.0*b*(sn.x+sn.y) - c*(sn.x*sn.x+sn.y*sn.y+tn.t*tn.t)
    //        - 0.1*(0.2*0.2 + 0.2*0.2 + tn.t*tn.t)
    //        - 0.1*(0.8*0.8 + 0.8*0.8 + tn.t*tn.t);

    return 1.0 - 2.0*b - c*(sn.x+sn.y+tn.t)
            - 0.1*(0.2 + 0.2 + tn.t)
            - 0.1*(0.8 + 0.8 + tn.t);
}

void LoadedHeatEquationFBVP::layerInfo(const DoubleVector&, const TimeNodePDE&) const
{}

void LoadedHeatEquationFBVP::layerInfo(const DoubleMatrix& u, const TimeNodePDE&) const
{
    IPrinter::printSeperatorLine();
    IPrinter::printMatrix(u);
}

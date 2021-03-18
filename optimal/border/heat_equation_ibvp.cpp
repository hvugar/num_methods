#include "heat_equation_ibvp.h"

#include "test_function.h"

#define HEAT_DIMENSION_2
#define HEAT_QUADRATIC
//#define HEAT_HOMOGENIOUS
//#define HEAT_DELTA

#if defined(HEAT_DIMENSION_1)
#define HEAT_LEFT_DIRICHLET
#define HEAT_RGHT_DIRICHLET
//#define HEAT_LEFT_NEUMANN
//#define HEAT_RGHT_NEUMAN
//#define HEAT_LEFT_ROBIN
//#define HEAT_RGHT_ROBIN
#endif

#if defined(HEAT_DIMENSION_2)
//#define HEAT_NORMAL_DIRICHLET
//#define HEAT_NORMAL_NEUMANN
#define HEAT_NORMAL_ROBIN
#endif

//double p_fx(const IParabolicFBVP *p, const SpaceNodePDE &sn, const TimeNodePDE &tn, int dt = 0, int dx = 0, int dy = 0);

const double fa = +1.0;   // must be plus for forward
const double fb = +0.4;   // must be minus or plus for forward -  some problems on high values
const double fc = -0.5;   // must be minus for forward

const double ba = -1.0;   // must be minus for backward
const double bb = -0.4;   // must be minus or plus for forward -  some problems on high values
const double bc = +0.5;   // must be plus for backward

void HeatEquationIBVP::Main(int argc, char *argv[])
{
    C_UNUSED(argc);
    C_UNUSED(argv);
#ifdef USE_LIB_IMAGING
    QGuiApplication app(argc, argv);
#endif

    HeatEquationIBVP h;
    h.setTimeDimension(Dimension(0.010, 100, 200));
    h.setSpaceDimensionX(Dimension(0.010, 200, 300));
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

//---------------------------------------------------------------------------------------------//

HeatEquationIBVP::HeatEquationIBVP() : IHeatEquationIBVP() {}

double HeatEquationIBVP::initial(const SpaceNodePDE &sn, InitialCondition) const
{
#if defined(HEAT_QUADRATIC)
    TimeNodePDE tn; tn.i = timeDimension().min(); tn.t = tn.i*timeDimension().step();
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
        condition = BoundaryConditionPDE(BoundaryCondition::Dirichlet); return TestFunction::u(tn, sn, TestFunction::FunctionValue);
#endif
#if defined(HEAT_LEFT_NEUMAN)
        condition = BoundaryConditionPDE(BoundaryCondition::Neumann); return TestFunction::u(tn, sn, TestFunction::SpaceFirstDerivativeX);
#endif
#if defined(HEAT_LEFT_ROBIN)
        condition = BoundaryConditionPDE::Robin(+1.0, -1.0); return condition.alpha()*TestFunction::u(tn, sn, TestFunction::FunctionValue)+condition.beta()*TestFunction::u(tn, sn, TestFunction::SpaceFirstDerivativeX);
#endif
    }
    if (sn.i == spaceDimensionX().max())
    {
#if defined(HEAT_RGHT_DIRICHLET)
        condition = BoundaryConditionPDE(BoundaryCondition::Dirichlet); return TestFunction::u(tn, sn, TestFunction::FunctionValue);
#endif
#if defined(HEAT_RGHT_NEUMAN)
        condition = BoundaryConditionPDE(BoundaryCondition::Neumann); return TestFunction::u(tn, sn, TestFunction::SpaceFirstDerivativeX);
#endif
#if defined(HEAT_RGHT_ROBIN)
        condition = BoundaryConditionPDE::Robin(+1.0, +1.0); return condition.alpha()*TestFunction::u(tn, sn, TestFunction::FunctionValue) +condition.beta()*TestFunction::u(tn, sn, TestFunction::SpaceFirstDerivativeX);
#endif
    }
#endif

#if defined(HEAT_DIMENSION_2)

#if defined(HEAT_NORMAL_DIRICHLET)
    condition = BoundaryConditionPDE(BoundaryCondition::Dirichlet);
    return TestFunction::u(tn, sn, TestFunction::FunctionValue);
#endif
#if defined(HEAT_NORMAL_NEUMANN)

#if defined(HEAT_QUADRATIC)
    //condition = BoundaryConditionPDE(BoundaryCondition::Neumann); return (::u_fx(this, sn, tn, -1, 3, 3));
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

    return TestFunction::u(tn, sn, TestFunction::TimeFirstDerivative)
            - TestFunction::u(tn, sn, TestFunction::SpaceSecondDerivativeX) * a
            - TestFunction::u(tn, sn, TestFunction::SpaceFirstDerivativeX) * b
            - TestFunction::u(tn, sn, TestFunction::FunctionValue) * c;
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
    IPrinter::printVector(15, 6, u);
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
}

//---------------------------------------------------------------------------------------------//

HeatEquationFBVP::HeatEquationFBVP() : IHeatEquationFBVP()
{}

double HeatEquationFBVP::final(const SpaceNodePDE &sn, FinalCondition) const
{
#if defined(HEAT_QUADRATIC)
    TimeNodePDE tn; tn.i = timeDimension().max(); tn.t = tn.i*timeDimension().step();
    return TestFunction::u(tn, sn, TestFunction::FunctionValue);
#endif
}

double HeatEquationFBVP::boundary(const SpaceNodePDE &sn, const TimeNodePDE &tn, BoundaryConditionPDE &condition) const
{
#if defined(HEAT_DIMENSION_1)
    if (sn.i == spaceDimensionX().min())
    {
#if defined(HEAT_LEFT_DIRICHLET)
        condition = BoundaryConditionPDE(BoundaryCondition::Dirichlet); return TestFunction::u(tn, sn, TestFunction::FunctionValue);
#endif
#if defined(HEAT_LEFT_ROBIN)
        condition = BoundaryConditionPDE::Robin(+1.0, -1.0);
        return condition.alpha()*TestFunction::u(tn, sn, TestFunction::FunctionValue)
                +condition.beta()*TestFunction::u(tn, sn, TestFunction::SpaceFirstDerivativeX);
#endif
    }
    if (sn.i == spaceDimensionX().max())
    {
#if defined(HEAT_RGHT_DIRICHLET)
        condition = BoundaryConditionPDE(BoundaryCondition::Dirichlet); return TestFunction::u(tn, sn, TestFunction::FunctionValue);
#endif
#if defined(HEAT_RGHT_ROBIN)
        condition = BoundaryConditionPDE::Robin(+1.0, +1.0);
        return condition.alpha()*TestFunction::u(tn, sn, TestFunction::FunctionValue)
                +condition.beta()*TestFunction::u(tn, sn, TestFunction::SpaceFirstDerivativeX);
#endif
    }
#endif

#if defined(HEAT_DIMENSION_2)

#if defined(HEAT_NORMAL_DIRICHLET)
    //condition = BoundaryConditionPDE(BoundaryCondition::Dirichlet);
    //return ::p_fx(this, sn, tn);
#endif
#if defined(HEAT_NORMAL_NEUMANN)
    //condition = BoundaryConditionPDE(BoundaryCondition::Neumann);
    //return (::p_fx(this, sn, tn)+::p_fx(this, sn, tn, -1, 3, 3));
#endif
#if defined(HEAT_NORMAL_ROBIN)
    //condition = BoundaryConditionPDE::Robin(+4.0, +2.0);
    //return (condition.alpha()*::p_fx(this, sn, tn)+condition.beta()*::p_fx(this, sn, tn, -1, 3, 3));
#endif
#endif
}

double HeatEquationFBVP::f(const SpaceNodePDE &sn, const TimeNodePDE &tn) const
{
#if defined(HEAT_DIMENSION_1)
    const double a = thermalDiffusivity();
    const double b = thermalConductivity();
    const double c = thermalConvection();

    return TestFunction::u(tn, sn, TestFunction::TimeFirstDerivative)
            - TestFunction::u(tn, sn, TestFunction::SpaceSecondDerivativeX) * a
            - TestFunction::u(tn, sn, TestFunction::SpaceFirstDerivativeX) * b
            - TestFunction::u(tn, sn, TestFunction::FunctionValue) * c;
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
#endif
}

void HeatEquationFBVP::layerInfo(const DoubleVector& u, const TimeNodePDE& tn) const
{
    C_UNUSED(u);
    C_UNUSED(tn);
    IPrinter::printVector(15, 6, u);
}

void HeatEquationFBVP::layerInfo(const DoubleMatrix &u, const TimeNodePDE &tn) const
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
}

/***********************************************************************************************************************************/

void LoadedHeatEquationIBVP::Main(int /*argc*/, char */*argv*/[])
{
    LoadedHeatEquationIBVP lheIBVP;
    lheIBVP.setThermalDiffusivity(1.0);
    lheIBVP.setThermalConductivity(0.0);
    lheIBVP.setThermalConvection(0.0);

    std::vector<LoadedSpacePoint> loadedPoints;

#ifdef HEAT_DIMENSION_1
    //loadedPoints.push_back(LoadedSpacePoint(0.20, 0.00, 0.00, -0.4));
    //loadedPoints.push_back(LoadedSpacePoint(0.80, 0.00, 0.00, -0.5));
    //loadedPoints.push_back(LoadedSpacePoint(0.50, 0.00, 0.00, -0.3));
#endif

#ifdef HEAT_DIMENSION_2
    //    loadedPoints.push_back(LoadedSpacePoint(2.30, 3.20, 0.00, 0.4));
    //    loadedPoints.push_back(LoadedSpacePoint(2.80, 3.70, 0.00, 0.6));
    //    loadedPoints.push_back(LoadedSpacePoint(2.60, 3.50, 0.00, 0.5));
    //    loadedPoints.push_back(LoadedSpacePoint(2.50, 3.50, 0.00, 0.4));


    loadedPoints.push_back(LoadedSpacePoint(2.30, 3.30, 0.00, 0.1));
    loadedPoints.push_back(LoadedSpacePoint(2.30, 3.80, 0.00, 0.2));
    //loadedPoints.push_back(LoadedSpacePoint(2.50, 3.10, 0.00, 0.3));
    //loadedPoints.push_back(LoadedSpacePoint(2.60, 3.60, 0.00, 0.4));
    //loadedPoints.push_back(LoadedSpacePoint(2.80, 3.30, 0.00, 0.5));
    //loadedPoints.push_back(LoadedSpacePoint(2.80, 3.90, 0.00, 0.6));

#endif

    lheIBVP.setLoadedPoints(loadedPoints);

#ifdef HEAT_DIMENSION_1
    lheIBVP.implicit_calculate_D1V1();
#endif

#ifdef HEAT_DIMENSION_2
    lheIBVP.implicit_calculate_D2V1();
#endif
}

double LoadedHeatEquationIBVP::initial(const SpaceNodePDE &sn, InitialCondition) const
{
    TimeNodePDE tn; tn.t = 0.0;
    return TestFunction::u(tn, sn, TestFunction::FunctionValue);
}

double LoadedHeatEquationIBVP::boundary(const SpaceNodePDE &sn, const TimeNodePDE &tn, BoundaryConditionPDE &bc) const
{
    bc = BoundaryConditionPDE::Dirichlet(); return TestFunction::u(tn, sn, TestFunction::FunctionValue);
    //bc = BoundaryConditionPDE::Dirichlet(1.0, 1.0); return sn.x*sn.x + sn.y*sn.y + tn.t*tn.t;
    //bc = BoundaryConditionPDE::Dirichlet(); return sn.x + sn.y + tn.t;
}

double LoadedHeatEquationIBVP::f(const SpaceNodePDE &sn, const TimeNodePDE &tn) const
{
    const double a =  thermalDiffusivity();
    const double b =  thermalConductivity();
    const double c =  thermalConvection();

    double r = 0.0;

#ifdef HEAT_DIMENSION_1
    r = TestFunction::u(tn, sn, TestFunction::TimeFirstDerivative)
            - TestFunction::u(tn, sn, TestFunction::SpaceSecondDerivativeX) * a
            - TestFunction::u(tn, sn, TestFunction::SpaceFirstDerivativeX) * b
            - TestFunction::u(tn, sn, TestFunction::FunctionValue) * c;

    for (size_t i=0; i<loadedPoints().size(); i++)
    {
        SpaceNodePDE sn2; sn2.x = loadedPoints().at(i).x;
        r -= TestFunction::u(tn, sn2, TestFunction::FunctionValue) * loadedPoints().at(i).d;
    }
#endif

#ifdef HEAT_DIMENSION_2
    r = TestFunction::u(tn, sn, TestFunction::TimeFirstDerivative)
            - TestFunction::u(tn, sn, TestFunction::SpaceSecondDerivativeX) * a
            - TestFunction::u(tn, sn, TestFunction::SpaceSecondDerivativeY) * a
            - TestFunction::u(tn, sn, TestFunction::SpaceFirstDerivativeX) * b
            - TestFunction::u(tn, sn, TestFunction::SpaceFirstDerivativeY) * b
            - TestFunction::u(tn, sn, TestFunction::FunctionValue) * c;

    for (size_t i=0; i<loadedPoints().size(); i++)
    {
        SpaceNodePDE sn2;
        sn2.x = loadedPoints().at(i).x;
        sn2.y = loadedPoints().at(i).y;
        r -= TestFunction::u(tn, sn2, TestFunction::FunctionValue) * loadedPoints().at(i).d;
    }
#endif

    return r;
}

void LoadedHeatEquationIBVP::layerInfo(const DoubleVector& u, const TimeNodePDE&) const
{
    IPrinter::printVector(u);
}

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

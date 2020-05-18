#include "heat_equation_ibvp.h"

using namespace p3p;

void HeatEquationIBVP1::Main(int argc, char *argv[])
{
    C_UNUSED(argc);
    C_UNUSED(argv);
#ifdef USE_LIB_IMAGING
    QGuiApplication app(argc, argv);
#endif

    HeatEquationIBVP1 h;
    h.setTimeDimension(Dimension(1.0, 0, 500));
    h.setSpaceDimensionX(Dimension(0.005, 0, 100));
    h.setSpaceDimensionY(Dimension(0.005, 0, 100));
    h.deltaGrid.initGrid(h.spaceDimensionX(), h.spaceDimensionY());
    h.deltaGrid.reset();

    h.setThermalDiffusivity(0.0001);
    h.setThermalConductivity(0.0);
    h.setThermalConvection(0.0);

    h.implicit_calculate_D2V1();
}

void HeatEquationFBVP1::Main(int argc UNUSED_PARAM, char *argv[] UNUSED_PARAM)
{
    C_UNUSED(argc);
    C_UNUSED(argv);
#ifdef USE_LIB_IMAGING
    QGuiApplication app(argc, argv);
#endif

    HeatEquationFBVP1 h;
    h.setTimeDimension(Dimension(0.01, 0, 10000));
    h.setSpaceDimensionX(Dimension(0.010, 100, 200));
    h.setSpaceDimensionY(Dimension(0.010, 200, 300));

    h.setThermalDiffusivity(0.000111);
    h.setThermalConductivity(0.0);
    h.setThermalConvection(0.0);

    h.implicit_calculate_D2V1();
}

HeatEquationIBVP1::HeatEquationIBVP1() : IHeatEquationIBVP()
{}

double HeatEquationIBVP1::initial(const SpaceNodePDE &sn, InitialCondition) const
{
    return 0.0;
}

double HeatEquationIBVP1::boundary(const SpaceNodePDE &, const TimeNodePDE &, BoundaryConditionPDE &condition) const
{
    condition = BoundaryConditionPDE(BoundaryCondition::Neumann, +0.0, +1.0, +0.0);
    return 0.0;
}

double HeatEquationIBVP1::f(const SpaceNodePDE &sn, const TimeNodePDE &tn) const
{
    double q = 0.01;
    double val = q*deltaGrid.gaussWeight(sn, SpacePoint(0.25, 0.25), 0.005, 0.005);
    //printf("%6d %6d %14.8f\n", sn.i, sn.j, val);
    return val;
}

void HeatEquationIBVP1::layerInfo(const DoubleMatrix &u, const TimeNodePDE &tn) const
{
    C_UNUSED(u);
    C_UNUSED(tn);

    if (tn.i % ((timeDimension().size())/10) == 0)
    {
        std::string msg = std::string("time: ") + std::to_string(tn.t)
                        + std::string(", min: ") + std::to_string(u.min())
                        + std::string(", max: ") + std::to_string(u.max());
        IPrinter::printSeperatorLine(msg.data());
        IPrinter::printMatrix(9, 4, u, 10, 10, msg.data());
    }
}

//---------------------------------------------------------------------------------------------//

HeatEquationFBVP1::HeatEquationFBVP1() : IHeatEquationFBVP() {}

double HeatEquationFBVP1::final(const SpaceNodePDE &, FinalCondition) const { return 0.0; }

double HeatEquationFBVP1::boundary(const SpaceNodePDE &sn, const TimeNodePDE &tn, BoundaryConditionPDE &condition) const
{
    condition = BoundaryConditionPDE(BoundaryCondition::Robin, +0.0, +2.0, +1.0);
    return 0.0;
}

double HeatEquationFBVP1::f(const SpaceNodePDE &, const TimeNodePDE &) const { return 0.0; }

void HeatEquationFBVP1::layerInfo(const DoubleMatrix &, const TimeNodePDE &) const {}

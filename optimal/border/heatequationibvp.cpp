#include "heatequationibvp.h"

void HeatEquationIBVP::Main(int argc UNUSED_PARAM, char *argv[] UNUSED_PARAM)
{
#ifdef USE_LIB_IMAGING
    QGuiApplication app(argc, argv);
#endif

    HeatEquationIBVP w;
    w.setThermalDiffusivity(1.0);
    w.setTimeDimension(Dimension(0.005, 0, 200));
    w.setSpaceDimensionX(Dimension(0.01, 0, 100));
    w.setSpaceDimensionY(Dimension(0.01, 0, 100));

    Benchmark bm;
    bm.tick();
    IPrinter::printSeperatorLine();
    bm.tock();
    bm.printDuration();
}

double HeatEquationIBVP::initial(const SpaceNodePDE &sn, InitialCondition condition) const
{
    return 0.0;
}

double HeatEquationIBVP::boundary(const SpaceNodePDE &sn, const TimeNodePDE &tn, BoundaryConditionPDE &condition) const
{
    return 0.0;
}

double HeatEquationIBVP::f(const SpaceNodePDE &sn, const TimeNodePDE &tn) const
{
    const double a = thermalDiffusivity();
    return 0.0;
}

void HeatEquationIBVP::layerInfo(const DoubleVector& u, const TimeNodePDE& tn) const
{
    C_UNUSED(u);
    C_UNUSED(tn);
}

void HeatEquationIBVP::layerInfo(const DoubleMatrix &u, const TimeNodePDE &tn) const
{
    C_UNUSED(u);
    C_UNUSED(tn);
}

//---------------------------------------------------------------------------------------------//

void FinalHeatEquationIBVP::Main(int argc UNUSED_PARAM, char *argv[] UNUSED_PARAM)
{
#ifdef USE_LIB_IMAGING
    QGuiApplication app(argc, argv);
#endif

    FinalHeatEquationIBVP w;
    w.setThermalDiffusivity(1.0);
    w.setTimeDimension(Dimension(0.005, 0, 200));
    w.setSpaceDimensionX(Dimension(0.01, 0, 100));
    w.setSpaceDimensionY(Dimension(0.01, 0, 100));

    Benchmark bm;
    bm.tick();
    IPrinter::printSeperatorLine();
    bm.tock();
    bm.printDuration();
}

double FinalHeatEquationIBVP::final(const SpaceNodePDE &sn, FinalCondition condition) const
{
    return 0.0;
}

double FinalHeatEquationIBVP::boundary(const SpaceNodePDE &sn, const TimeNodePDE &tn, BoundaryConditionPDE &condition) const
{
    return 0.0;
}

double FinalHeatEquationIBVP::f(const SpaceNodePDE &sn, const TimeNodePDE &tn) const
{
    const double a = thermalDiffusivity();
    return 0.0;
}

void FinalHeatEquationIBVP::layerInfo(const DoubleVector& u, const TimeNodePDE& tn) const
{
    C_UNUSED(u);
    C_UNUSED(tn);
}

void FinalHeatEquationIBVP::layerInfo(const DoubleMatrix &u, const TimeNodePDE &tn) const
{
    C_UNUSED(u);
    C_UNUSED(tn);
}


#include "solver2.h"

using namespace p3p0;

void Solver2::Main(int /*argc*/, char **/*argv*/)
{
    HeatEquationIBVP he;
    he.setThermalDiffusivity(1.0);
    he.setThermalConvection(0.0);
    he.setThermalConductivity(0.0);
    he.implicit_calculate_D1V1();
}

HeatEquationIBVP::HeatEquationIBVP() {}

HeatEquationIBVP::HeatEquationIBVP(const HeatEquationIBVP &) {}

HeatEquationIBVP & HeatEquationIBVP::operator =(const HeatEquationIBVP &) { return *this; }

HeatEquationIBVP::~HeatEquationIBVP() {}

auto HeatEquationIBVP::q(size_t i, const TimeNodePDE &tn) const -> double
{
    if (i==0) return tn.t;
    if (i==1) return tn.t;
    throw std::exception();
}

auto HeatEquationIBVP::v(size_t /*i*/, const PointNodeODE &/*tn*/, SpacePoint &/*vl*/) const -> void { }

auto HeatEquationIBVP::initial(const SpaceNodePDE &/*sn*/, InitialCondition /*cn*/) const -> double { return 0.0; }

auto HeatEquationIBVP::boundary(const SpaceNodePDE &sn, const TimeNodePDE &/*tn*/, BoundaryConditionPDE &cn) const -> double
{
    if (sn.i == spaceDimensionX().min()) { cn = BoundaryConditionPDE::Robin(lambda1, -1.0, lambda1); return theta; }
    if (sn.i == spaceDimensionX().max()) { cn = BoundaryConditionPDE::Robin(lambda1, +1.0, lambda1); return theta; }
    throw std::exception();
}

auto HeatEquationIBVP::f(const SpaceNodePDE &sn, const TimeNodePDE &tn) const -> double
{
    auto fx = 0.0;
    fx += q(0, tn)*DeltaFunction::gaussian(sn, SpacePoint(0.25), SpacePoint(0.01));
    fx += q(1, tn)*DeltaFunction::gaussian(sn, SpacePoint(0.75), SpacePoint(0.01));
    return fx;
}

auto HeatEquationIBVP::layerInfo(const DoubleVector &u, const TimeNodePDE &tn) const -> void
{
    IPrinter::printVector(u, std::to_string(tn.i).data());
}

auto HeatEquationIBVP::timeDimension() const -> Dimension { return Dimension(0.001, 0, 1000); }

auto HeatEquationIBVP::spaceDimensionX() const -> Dimension { return Dimension(0.01, 0, 100); }

auto HeatEquationIBVP::spaceDimensionY() const -> Dimension { return 0.0; }

auto HeatEquationIBVP::spaceDimensionZ() const -> Dimension { return 0.0; }

Solver2::Solver2()
{}


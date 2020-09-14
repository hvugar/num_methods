#include "solver2.h"

using namespace p3p0;

void Solver2::Main(int /*argc*/, char **/*argv*/)
{
    HeatEquationIBVP he;
    he.s = new Solver2();
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

auto HeatEquationIBVP::z(size_t _i, const TimeNodePDE &/*tn*/) const -> double
{
    if (_i==0) return SpacePoint(0.30).x;
    if (_i==1) return SpacePoint(0.70).x;
    throw std::exception();
}

auto HeatEquationIBVP::initial(const SpaceNodePDE &/*sn*/, InitialCondition /*cn*/) const -> double { return 0.0; }

auto HeatEquationIBVP::boundary(const SpaceNodePDE &sn, const TimeNodePDE &/*tn*/, BoundaryConditionPDE &bc) const -> double
{
    if (sn.i == spaceDimensionX().min()) { bc = BoundaryConditionPDE::Robin(s->lambda1(), -1.0, s->lambda1()); return s->theta(); }
    if (sn.i == spaceDimensionX().max()) { bc = BoundaryConditionPDE::Robin(s->lambda1(), +1.0, s->lambda1()); return s->theta(); }
    throw std::exception();
}

auto HeatEquationIBVP::f(const SpaceNodePDE &sn, const TimeNodePDE &tn) const -> double
{
    auto fx = 0.0;
    fx += q(0, tn)*DeltaFunction::gaussian(sn.x, z(0, tn), SpacePoint(0.01).x);
    fx += q(1, tn)*DeltaFunction::gaussian(sn.x, z(1, tn), SpacePoint(0.01).x);
    return fx;
}

auto HeatEquationIBVP::layerInfo(const DoubleVector &u, const TimeNodePDE &tn) const -> void
{
    //IPrinter::printVector(u, std::to_string(tn.i).data());
    if (static_cast<int>(tn.i) == timeDimension().max()) { s->U = u; }
}

auto HeatEquationIBVP::timeDimension() const -> Dimension { return s->timeDimension(); }

auto HeatEquationIBVP::spaceDimensionX() const -> Dimension { return s->spaceDimensionX(); }

auto HeatEquationIBVP::spaceDimensionY() const -> Dimension { return spaceDimensionX(); }

auto HeatEquationIBVP::spaceDimensionZ() const -> Dimension { return spaceDimensionX(); }

/**************************************************************************************************************************/

HeatEquationFBVP::HeatEquationFBVP() {}

HeatEquationFBVP::HeatEquationFBVP(const HeatEquationFBVP &) {}

HeatEquationFBVP & HeatEquationFBVP::operator =(const HeatEquationFBVP &) { return *this; }

HeatEquationFBVP::~HeatEquationFBVP() {}

auto HeatEquationFBVP::final(const SpaceNodePDE &sn, FinalCondition /*ic*/) const -> double
{
    size_t n = static_cast<size_t>(sn.i);
    return 2.0*s->mu(sn) * (s->U.at(n) - s->V.at(n));
}

auto HeatEquationFBVP::boundary(const SpaceNodePDE &sn, const TimeNodePDE &/*tn*/, BoundaryConditionPDE &bc) const -> double
{
    if (sn.i == spaceDimensionX().min()) { bc = BoundaryConditionPDE::Robin(s->lambda1(), -1.0, s->lambda1()); return s->theta(); }
    if (sn.i == spaceDimensionX().max()) { bc = BoundaryConditionPDE::Robin(s->lambda1(), +1.0, s->lambda1()); return s->theta(); }
    throw std::exception();
}

auto HeatEquationFBVP::f(const SpaceNodePDE &/*sn*/, const TimeNodePDE &/*tn*/) const -> double
{
    return 0.0;
}

auto HeatEquationFBVP::HeatEquationFBVP::layerInfo(const DoubleVector &/*u*/, const TimeNodePDE &/*tn*/) const -> void {}

auto HeatEquationFBVP::timeDimension() const -> Dimension { return s->timeDimension(); }

auto HeatEquationFBVP::spaceDimensionX() const -> Dimension { return s->spaceDimensionX(); }

auto HeatEquationFBVP::spaceDimensionY() const -> Dimension { return spaceDimensionX(); }

auto HeatEquationFBVP::spaceDimensionZ() const -> Dimension { return spaceDimensionX(); }

/**************************************************************************************************************************/

auto Common::convert(const DoubleVector &x, size_t size, size_t length) -> void
{
    if (q != nullptr)
    {
        for (size_t i=0; i<size; i++)
        {
            delete [] q[i];
        }
        delete [] q;
    }

    q = new double* [size];

    for (size_t i=0; i<size; i++)
    {
        q[i] = new double [length];
        for (size_t n=0; n<length; n++)
        {
            q[i][n] = x[i*length+n];
        }
    }
}

auto Common::convert(size_t size, size_t length, DoubleVector &x) const -> void
{
    x.clear();
    x.resize(size*length);

    if (q != nullptr)
    {
        for (size_t i=0; i<size; i++)
        {
            for (size_t n=0; n<length; n++)
            {
                x[i*length+n] = q[i][n];
            }
        }
    }
}

/**************************************************************************************************************************/

Solver2::Solver2()
{
    V.resize(spaceDimensionX().size(), 3.0);
}

auto Solver2::gradient(const DoubleVector &x, DoubleVector &/*g*/) const -> void
{
    const unsigned int sizeX = spaceDimensionX().size();
    const_cast<Solver2*>(this)->convert(x, 2, sizeX);

    forward.implicit_calculate_D1V1();
    backward.implicit_calculate_D1V1();


}

auto Solver2::fx(const DoubleVector &x) const -> double
{
    const_cast<Solver2*>(this)->convert(x, 2, timeDimension().size());
    forward.implicit_calculate_D1V1();

    const unsigned int sizeX = spaceDimensionX().size()-1;
    const double stepX = spaceDimensionX().step();

    double sum = (U[0]-V[0])*(U[0]-V[0]);
    for (unsigned int i=1; i<sizeX; i++)
    {
        sum += 0.5*(U[i]-V[i])*(U[i]-V[i]);
    }
    sum += (U[sizeX]-V[sizeX])*(U[sizeX]-V[sizeX]);

    return sum*stepX;
}

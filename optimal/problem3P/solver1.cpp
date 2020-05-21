#include "solver1.h"

using namespace p3p;

void Solver1::Main(int argc, char **argv)
{}

Solver1::Solver1() {}

Solver1::~Solver1() {}

double Solver1::frw_initial(const SpaceNodePDE &, InitialCondition) const { return 0.0; }

double Solver1::frw_boundary(const SpaceNodePDE &, const TimeNodePDE &, BoundaryConditionPDE &cn) const
{
    cn = BoundaryConditionPDE(BoundaryCondition::Neumann, 0.0, 1.0, 1.0);
    return 0.0;
}

auto Solver1::frw_f(const SpaceNodePDE &sn, const TimeNodePDE &tn) const -> double
{
    unsigned int ln = static_cast<unsigned int>(tn.i) / 2;
    double q = externalSource[ln].q;
    const SpacePoint &z = externalSource[ln].z;
    return q * deltaZ.gaussWeight(sn, z, SIGMA_X, SIGMA_Y);

    //    const Dimension &dimensionX = spaceDimensionX();
    //    const Dimension &dimensionY = spaceDimensionY();
    //    const double hx = dimensionX.step();
    //    const double hy = dimensionY.step();

    //    if (sn.i==50 && sn.j==50)
    //        return mq[ln]*(1.0/(hx*hy));
    //    else
    //        return 0.0;

}

void Solver1::frw_layerInfo(const DoubleMatrix &u, const TimeNodePDE &tn) const
{
    unsigned int ln = tn.i;
    Solver1* solver = const_cast<Solver1*>(this);

    for (unsigned int j=0; j<measurePointNumber; j++)
    {
        const SpacePoint& measurePoint = measurePoints[j];
        solver->measurePointValues[j].z = DeltaFunction::lumpedPointG(u, measurePoint, _spaceDimensionX, _spaceDimensionY, 1, 4);
    }

    for (unsigned int i=0; i<heatSourceNumber; i++)
    {
        q[i] = v[i] = 0.0;
        SpacePoint* heatSourceRoute = heatSourceRoutes[i];
        const SpacePoint& zi = heatSourceRoute[ln];

        for (unsigned int j=0; j<measurePointNumber; j++)
        {
            double _alpha1 = alpha1[i*heatSourceNumber+j];
            double _alpha2 = alpha2[i*heatSourceNumber+j];
            double _alpha3 = alpha3[i*heatSourceNumber+j];
            double _betta1 = betta1[i*heatSourceNumber+j];
            double _betta2 = betta2[i*heatSourceNumber+j];
            double _betta3 = betta3[i*heatSourceNumber+j];
            double _uij = uij[i*heatSourceNumber+j];
            const SpacePoint& mp = measurePoints[j];

            double distance = (zi.x - mp.x)*(zi.x - mp.x) + (zi.y - mp.y)*(zi.y - mp.y);

            q[i] += (_alpha1*distance*distance + _alpha2*distance + _alpha3) * ( solver->measurePointValues[j].z - _uij );
            v[i] += (_betta1*distance*distance + _betta2*distance + _betta3) * ( solver->measurePointValues[j].z - _uij );
        }
    }
}

//--------------------------------------------------------------------------------------------------------------//

HeatEquationIBVP::HeatEquationIBVP(Solver1 *solver) { this->solver = solver; }

HeatEquationIBVP::~HeatEquationIBVP() {}

auto HeatEquationIBVP::initial(const SpaceNodePDE &sn, InitialCondition cn) const -> double { return solver1->frw_initial(sn, cn); }

auto HeatEquationIBVP::boundary(const SpaceNodePDE &sn, const TimeNodePDE &tn, BoundaryConditionPDE &cn) const -> double { return solver->frw_boundary(sn, tn, cn); }

auto HeatEquationIBVP::f(const SpaceNodePDE &sn, const TimeNodePDE &tn) const-> double { return solver1->frw_f(sn, tn); }

auto HeatEquationIBVP::layerInfo(const DoubleMatrix &u, const TimeNodePDE &tn) const -> void
{
    const_cast<HeatEquationIBVP*>(this)->lastLayerU = u;
    solver->frw_layerInfo(u, tn);
}

auto HeatEquationIBVP::weight() const -> double { throw std::runtime_error(""); }

auto HeatEquationIBVP::timeDimension() const -> Dimension { return solver1->timeDimension(); }

auto HeatEquationIBVP::spaceDimensionX() const -> Dimension { return solver1->spaceDimensionX(); }

auto HeatEquationIBVP::spaceDimensionY() const -> Dimension { return solver1->spaceDimensionY(); }

auto HeatEquationIBVP::spaceDimensionZ() const -> Dimension { throw std::runtime_error(""); }

//--------------------------------------------------------------------------------------------------------------//

auto HeatEquationIBVP::A(const PointNodeODE &, unsigned int, unsigned int) const -> double
{
    return 0.0;
}

auto HeatEquationIBVP::B(const PointNodeODE &, unsigned int, unsigned int) const -> double
{
    return 0.0;
}

auto HeatEquationIBVP::C(const PointNodeODE &node, unsigned int row) const -> double
{
    return 0.0 + D(node, row);
}

auto HeatEquationIBVP::D(const PointNodeODE &node, unsigned int row) const -> double
{
    const double vl[2] = {+1.0, +1.0};
    return vl[row-1]*node.x;
}

auto HeatEquationIBVP::initial(InitialCondition, unsigned int row) const -> double
{
    const double val[2] = { 0.10, 0.10 };
    return val[row-1];
}

auto HeatEquationIBVP::count() const -> unsigned int
{
    return 2;
}

auto HeatEquationIBVP::dimension() const -> Dimension { return solver->timeDimension(); }

auto HeatEquationIBVP::iterationInfo(const DoubleVector &z, const PointNodeODE &node) const -> void
{
    unsigned int ln = static_cast<unsigned int>(node.i);
    ExternalSource &es = const_cast<Solver*>(solver)->externalSource[ln];
    es.z.x = z[0];
    es.z.y = z[1];
}

//--------------------------------------------------------------------------------------------------------------//

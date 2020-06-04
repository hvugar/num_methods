#include "solver1.h"

using namespace p3p1;

void Solver1::Main(int argc, char **argv)
{

}

Solver1::Solver1() {}

Solver1::~Solver1() {}

void Solver1::setPointNumber(size_t heatSourceNumber, size_t measrPointNumber)
{
    this->heatSourceNumber = heatSourceNumber;
    this->measrPointNumber = measrPointNumber;

    alpha1.resize(heatSourceNumber, measrPointNumber, 0.13);
    alpha2.resize(heatSourceNumber, measrPointNumber, 0.12);
    alpha3.resize(heatSourceNumber, measrPointNumber, 0.15);
    betta1.resize(heatSourceNumber, measrPointNumber, 0.12);
    betta2.resize(heatSourceNumber, measrPointNumber, 0.18);
    betta3.resize(heatSourceNumber, measrPointNumber, 0.14);

    this->measurePointValues[0] = SpacePoint(0.35, 0.72);
    this->measurePointValues[1] = SpacePoint(0.68, 0.12);
    this->measurePointValues[2] = SpacePoint(0.28, 0.84);
    this->measurePointValues[3] = SpacePoint(0.55, 0.62);

    DoubleMatrix nominU;
}

double Solver1::frw_initial(const SpaceNodePDE &, InitialCondition) const { return _initialValue; }

double Solver1::frw_boundary(const SpaceNodePDE &, const TimeNodePDE &, BoundaryConditionPDE &cn) const
{
    cn = BoundaryConditionPDE(BoundaryCondition::Neumann, 0.0, 1.0, 1.0);
    return 0.0;
}

auto Solver1::frw_f(const SpaceNodePDE &sn, const TimeNodePDE &tn) const -> double
{
    double fx = 0.0;
    size_t ln = static_cast<size_t>(tn.i);

    for (size_t i=0; i<heatSourceNumber; i++)
    {
        fx += sourceParams[i][ln].q * DeltaFunction::gaussian(sn, sourceParams[i]->z, SpacePoint(spaceDimensionX().step(), spaceDimensionY().step()));
    }

    return fx;
}

void Solver1::frw_layerInfo(const DoubleMatrix &u, const TimeNodePDE &tn) const
{
}

//--------------------------------------------------------------------------------------------------------------//

HeatEquationIBVP::HeatEquationIBVP(Solver1 *solver) { this->solver = solver; }

HeatEquationIBVP::~HeatEquationIBVP() {}

auto HeatEquationIBVP::initial(const SpaceNodePDE &sn, InitialCondition cn) const -> double { return solver->frw_initial(sn, cn); }

auto HeatEquationIBVP::boundary(const SpaceNodePDE &sn, const TimeNodePDE &tn, BoundaryConditionPDE &cn) const -> double { return solver->frw_boundary(sn, tn, cn); }

auto HeatEquationIBVP::f(const SpaceNodePDE &sn, const TimeNodePDE &tn) const-> double { return solver->frw_f(sn, tn); }

auto HeatEquationIBVP::layerInfo(const DoubleMatrix &u, const TimeNodePDE &tn) const -> void
{
    HeatEquationIBVP *const_this = const_cast<HeatEquationIBVP*>(this);
    size_t ln = static_cast<size_t>(tn.i);
    DoubleVector z0(2), z1(2);
    PointNodeODE n0, n1;

    DeltaGrid2D **deltaGrid = new DeltaGrid2D* [solver->measrPointNumber];
    for (size_t j=0; j<solver->measrPointNumber; j++)
    {
        deltaGrid[j] = new DeltaGrid2D();
        deltaGrid[j]->initGrid(solver->spaceDimensionX(), solver->spaceDimensionY());
    }

    if (ln == 0)
    {
        n0.i = static_cast<int>(ln+0); n0.x = n0.i*timeDimension().step();
        n1.i = static_cast<int>(ln+1); n1.x = n1.i*timeDimension().step();

        for (size_t i=0; i<solver->heatSourceNumber; i++)
        {
            double qi = 0.0;
            double vi = 0.0;

            const_this->i = i+1;
            const_this->start( z0, n0 );
            solver->sourceParams[i][ln].z = SpacePoint(z0[0], z0[1]);

            for (size_t j=0; j<solver->measrPointNumber; j++)
            {
                double uj = DeltaFunction::lumpedPoint3(u, solver->measurePoints[j], spaceDimensionX(), spaceDimensionY());

                const SpacePoint &mp = solver->measurePoints[j];
                const SpacePoint &zi = solver->sourceParams[i][ln].z;
                double dist = sqrt((zi.x - mp.x)*(zi.x - mp.x) + (zi.y - mp.y)*(zi.y - mp.y));
                qi += (solver->alpha1[i][j]*dist*dist + solver->alpha2[i][j]*dist + solver->alpha3[i][j]) * (uj - solver->nominU[i][j]);
                vi += (solver->betta1[i][j]*dist*dist + solver->betta2[i][j]*dist + solver->betta3[i][j]) * (uj - solver->nominU[i][j]);
            }
            solver->sourceParams[i][ln].q = qi;
            solver->sourceParams[i][ln].v = vi;
        }
    }

    for (size_t i=0; i<solver->heatSourceNumber; i++)
    {
        double qi = 0.0;
        double vi = 0.0;

        const_this->i = i+1;
        z0[0] = solver->sourceParams[i][ln].z.x;
        z0[1] = solver->sourceParams[i][ln].z.y;
        const_this->next( z0, n0, z1, n1 );
        solver->sourceParams[i][ln+1].z = SpacePoint(z1[0], z1[1]);

        for (size_t j=0; j<solver->measrPointNumber; j++)
        {
            double uj = DeltaFunction::lumpedPoint3(u, solver->measurePoints[j], spaceDimensionX(), spaceDimensionY());

            const SpacePoint &mp = solver->measurePoints[j];
            const SpacePoint &zi = solver->sourceParams[i][ln+1].z;
            double dist = sqrt((zi.x - mp.x)*(zi.x - mp.x) + (zi.y - mp.y)*(zi.y - mp.y));
            qi += (solver->alpha1[i][j]*dist*dist + solver->alpha2[i][j]*dist + solver->alpha3[i][j]) * (uj - solver->nominU[i][j]);
            vi += (solver->betta1[i][j]*dist*dist + solver->betta2[i][j]*dist + solver->betta3[i][j]) * (uj - solver->nominU[i][j]);
        }

        solver->sourceParams[i][ln+1].q = qi;
        solver->sourceParams[i][ln+1].v = vi;
    }

    solver->frw_layerInfo(u, tn);
}

auto HeatEquationIBVP::timeDimension() const -> Dimension { return solver->timeDimension(); }

auto HeatEquationIBVP::spaceDimensionX() const -> Dimension { return solver->spaceDimensionX(); }

auto HeatEquationIBVP::spaceDimensionY() const -> Dimension { return solver->spaceDimensionY(); }

auto HeatEquationIBVP::spaceDimensionZ() const -> Dimension { throw std::runtime_error(""); }

//--------------------------------------------------------------------------------------------------------------//

auto HeatEquationIBVP::A(const PointNodeODE &, size_t r, size_t c) const -> double
{
    if (i==1) { const double data[2][2] = { {+0.1, +0.0}, {+0.0, -0.2} }; return data[r-1][c-1]; }
    if (i==2) { const double data[2][2] = { {+0.6, +0.0}, {+0.0, +0.4} }; return data[r-1][c-1]; }
    return 0.0;
}

auto HeatEquationIBVP::B(const PointNodeODE &, size_t r, size_t c) const -> double
{
    if (i==1) { const double data[2][2] = { {+0.0, +0.0}, {+0.0, +0.0} }; return data[r-1][c-1]; }
    if (i==2) { const double data[2][2] = { {+0.0, +0.0}, {+0.0, +0.0} }; return data[r-1][c-1]; }
    throw std::runtime_error(std::string("auto HeatEquationIBVP::B(const PointNodeODE &, size_t r, size_t c) const -> double"));
}

auto HeatEquationIBVP::C(const PointNodeODE &node, size_t r) const -> double
{
    size_t ln = static_cast<size_t>(node.i);
    if (i==1) { const double data[2] = { +0.1, -0.4 }; return data[r-1]*solver->sourceParams[i-1][ln].v; }
    if (i==2) { const double data[2] = { +0.2, -0.3 }; return data[r-1]*solver->sourceParams[i-1][ln].v; }
    throw std::runtime_error(std::string("auto HeatEquationIBVP::C(const PointNodeODE &node, size_t r) const -> double"));
}

auto HeatEquationIBVP::initial(InitialCondition c, size_t r) const -> double
{
    if (c == InitialCondition::InitialValue)
    {
        if (i==1) { const double data[2] = { 0.10, 0.90 }; return data[r-1]; }
        if (i==2) { const double data[2] = { 0.10, 0.10 }; return data[r-1]; }
    }
    if (c == InitialCondition::InitialFirstDerivative) return 0.0;
    throw std::runtime_error(std::string("auto HeatEquationIBVP::initial(InitialCondition, size_t row) const -> double"));
}

auto HeatEquationIBVP::count() const -> size_t
{
    return 2;
}

auto HeatEquationIBVP::dimension() const -> Dimension { return solver->timeDimension(); }

auto HeatEquationIBVP::iterationInfo(const DoubleVector &z, const PointNodeODE &node) const -> void
{
    //    unsigned int ln = static_cast<unsigned int>(node.i);
    //    ExternalSource &es = const_cast<Solver*>(solver)->externalSource[ln];
    //    es.z.x = z[0];
    //    es.z.y = z[1];
}

//--------------------------------------------------------------------------------------------------------------//

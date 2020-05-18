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

    for (unsigned int j=0; j<measurePointNumber; j++)
    {
        DeltaFunction::gaussian()
    }

    double q = externalSource[ln].q;
    const SpacePoint &z = externalSource[ln].z;
    return q * deltaZ.gaussWeight(sn, z, SIGMA_X, SIGMA_Y);

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

auto HeatEquationIBVP::A(const PointNodeODE &, unsigned int row, unsigned int col) const -> double
{
#ifdef VARIANT_1
    const double mx[2][2] = {
        {-3.0, -2.0}, {-4.0, -5.0}
    };
    return mx[row-1][col-1];
#endif
#ifdef VARIANT_2
    const double vl[4][4] = {
        {+0.0, +0.0, +1.0, +0.0},
        {+0.0, +0.0, +0.0, +1.0},
        {+0.0, +0.0, +0.0, +0.0},
        {+0.0, +0.0, +0.0, +0.0},
    };
    return vl[row-1][col-1];
#endif
#ifdef VARIANT_3
    const double vl[2][2] = { {+0.0, +0.0}, {+0.0, +0.0} };
    return vl[row-1][col-1];
#endif
}

auto HeatEquationIBVP::B(const PointNodeODE &node, unsigned int row) const -> double
{
    //return solver->zt(node, row)
    //        - (A(node, row, 1)*solver->z(node, 1)+A(node, row, 2)*solver->z(node, 2))
    //        - C(node, row)*solver->v(node)
    //        + C(node, row)*solver->mv[node.i];

    double t = node.x;
#ifdef VARIANT_1
    double v = solver->v(node);
    return (row == 1) ?
                0.4*M_PI*sin(2.0*M_PI*t) - 1.8*M_PI*t*sin(4.0*M_PI*t*t)
                + 3.0*(0.40*sin(1.0*M_PI*t)*sin(1.0*M_PI*t)+0.45*cos(2.0*M_PI*t*t)*cos(2.0*M_PI*t*t)+0.05)
                + 2.0*(0.30*sin(4.0*M_PI*t)*sin(4.0*M_PI*t)+0.45*cos(3.0*M_PI*t*t)*cos(3.0*M_PI*t*t)+0.05) - 5.0*sin(2.0*M_PI*t*t) + 5.0*v:
                1.2*M_PI*sin(8.0*M_PI*t) - 2.7*M_PI*t*sin(6.0*M_PI*t*t)
                + 4.0*(0.40*sin(1.0*M_PI*t)*sin(1.0*M_PI*t)+0.45*cos(2.0*M_PI*t*t)*cos(2.0*M_PI*t*t)+0.05)
                + 5.0*(0.30*sin(4.0*M_PI*t)*sin(4.0*M_PI*t)+0.45*cos(3.0*M_PI*t*t)*cos(3.0*M_PI*t*t)+0.05) - 4.0*sin(2.0*M_PI*t*t) + 4.0*v;
#endif
#ifdef VARIANT_2
    puts("B1");
    //const double mx[4] = {+0.0, +0.0, +1.0, +1.0};
    unsigned int ln = static_cast<unsigned int>(node.i);
    double aa = C(node, row-1);//*solver->externalSource[ln].v;
    puts("B1");
    return aa;
#endif
#ifdef VARIANT_3
    const double vl[2] = {+0.0, +0.0};
    unsigned int ln = static_cast<unsigned int>(node.i);
    double ret = C(node, row)*solver->externalSource[ln].v + vl[row-1];
    return ret;
#endif
}

auto HeatEquationIBVP::C(const PointNodeODE &node, unsigned int row) const -> double
{
#ifdef VARIANT_1
    const double c[2] = { +5.0, +4.0 };
    return c[row-1];
#endif
#ifdef VARIANT_2
    puts("C1");
    const double vl[4] = { +0.0, +0.0, +1.0, +1.0 };
    double bb = 0.0;//vl[row-1];
    printf("%f row %d\n", bb, row);
    return bb;
#endif
#ifdef VARIANT_3
    const double vl[2] = {+1.0, +1.0};
    return vl[row-1]*node.x;
#endif
}

auto HeatEquationIBVP::initial(InitialCondition, unsigned int row) const -> double
{
#ifdef VARIANT_1
    const double val[2] = { 0.50, 0.50 };
    return val[row-1];
#endif
#ifdef VARIANT_2
    const double val[4] = { 0.10, 0.10, 0.00, 0.00 };
    return val[row-1];
#endif
#ifdef VARIANT_3
    const double val[2] = { 0.10, 0.10 };
    return val[row-1];
#endif
}

auto HeatEquationIBVP::count() const -> unsigned int
{
#ifdef VARIANT_1
    return 2;
#endif
#ifdef VARIANT_2
    return 4;
#endif
#ifdef VARIANT_3
    return 2;
#endif
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

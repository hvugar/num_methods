#include "problem1p_solver.h"

using namespace p1p;

void ProblemSolver::Main(int argc, char* argv[])
{
    C_UNUSED(argc);
    C_UNUSED(argv);

    ProblemSolver solver(Dimension(0.001, 0, 100), Dimension(0.001, 0, 100));

    DoubleVector x(101, 3.0);
    DoubleVector ng(101, 0.0);
    DoubleVector ag(101, 0.0);

    solver.gradient(x, ag);
    IPrinter::printVector(ag.EuclideanNormalize());

    IGradient::Gradient(&solver, 0.0001, x, ng);
    IPrinter::printVector(ng.EuclideanNormalize());
}

ProblemSolver::ProblemSolver(const Dimension &timeDimension, const Dimension &spaceDimensionX)
{
    setTimeDimension(timeDimension);
    setSpaceDimensionX(spaceDimensionX);

    params.initialTemperature = 0.0;
    params.environmentTemperature = 0.5;
    params.thermalConductivity0 = 0.0;
    params.thermalConductivity1 = 0.1;
    params.thermalConductivity2 = 0.001;

    forward.setThermalDiffusivity(params.thermalDiffusivity);
    forward.setThermalConductivity(params.thermalConductivity0);
    forward.setThermalConvection(0.0);

    backward.setThermalDiffusivity(-params.thermalDiffusivity);
    backward.setThermalConductivity(params.thermalConductivity0);
    backward.setThermalConvection(0.0);

    forward.solver = this;
    backward.solver = this;
    const_this = const_cast<ProblemSolver*>(this);
}

void ProblemSolver::setTimeDimension(const Dimension &timeDimension)
{
    this->_timeDimension = timeDimension;
    forward.setTimeDimension(timeDimension);
    backward.setTimeDimension(timeDimension);

    uint32_t length = static_cast<uint32_t>(timeDimension.size()+1);
    p0.resize(length);
    params.v =  new double [length];
}

void ProblemSolver::setSpaceDimensionX(const Dimension &spaceDimensionX)
{
    this->_spaceDimensionX = spaceDimensionX;
    forward.setSpaceDimensionX(spaceDimensionX);
    backward.setSpaceDimensionX(spaceDimensionX);

    uint32_t length = static_cast<uint32_t>(spaceDimensionX.size()+1);
    u.resize(length);
    U.resize(length);
}

void ProblemSolver::gradient(const DoubleVector &x, DoubleVector &g) const
{
    C_UNUSED(x);
    C_UNUSED(g);

    const unsigned int length = x.length();
    for (unsigned int i=0; i<length; i++) const_this->params.v[i] = x[i];
    g.resize(length);

    forward.implicit_calculate_D1V1CN();
    backward.implicit_calculate_D1V1CN();

    for (unsigned int i=0; i<length; i++) g[i] = -params.thermalDiffusivity*p0[i];
    g[0] = 0.0;
}

double ProblemSolver::fx(const DoubleVector &x) const
{
    ProblemSolver *const_solver = const_cast<ProblemSolver *>(this);
    const_solver->params.v = const_cast<double *>(x.data());

    auto hx = _spaceDimensionX.step();
    auto N = static_cast<uint32_t>(_timeDimension.size());
    forward.implicit_calculate_D1V1CN();

    double sum = 0.0;
    sum += 0.5*(u[0]-U[0])*(u[0]-U[0]);
    for (unsigned int n=1; n<=N-1; n++)
    {
        sum += (u[n]-U[n])*(u[n]-U[n]);
    }
    sum += 0.5*(u[N]-U[N])*(u[N]-U[N]);
    sum *= hx;
    return sum;
}

/*********************************************************************************************************************************/

double HeatEquationIBVP::initial(const SpaceNodePDE &sn, InitialCondition condition) const
{
    C_UNUSED(sn);
    C_UNUSED(condition);
    return solver->params.initialTemperature;
}

double HeatEquationIBVP::boundary(const SpaceNodePDE &sn, const TimeNodePDE &tn, BoundaryConditionPDE &condition) const
{
    C_UNUSED(sn);
    C_UNUSED(tn);
    C_UNUSED(condition);

    if (sn.i == 0)
    {
        const EquationParameters &params = solver->params;
        condition = BoundaryConditionPDE(BoundaryCondition::Robin, params.thermalConductivity1, -1.0, params.thermalConductivity1);
        return params.v[tn.i];
    }
    else
    {
        const EquationParameters &params = solver->params;
        condition = BoundaryConditionPDE(BoundaryCondition::Robin, params.thermalConductivity2, +1.0, params.thermalConductivity2);
        return params.environmentTemperature;
    }
}

double HeatEquationIBVP::f(const SpaceNodePDE &sn, const TimeNodePDE &tn) const
{
    C_UNUSED(sn);
    C_UNUSED(tn);
    const EquationParameters &params = solver->params;
    return params.thermalConductivity0 * params.environmentTemperature;
}

void HeatEquationIBVP::layerInfo(const DoubleVector &u, const TimeNodePDE &tn) const
{
    C_UNUSED(u);
    C_UNUSED(tn);
    Dimension _timeDimension = timeDimension();
    if (_timeDimension.size() == static_cast<int>(tn.i)) const_cast<ProblemSolver*>(solver)->U = u;
}

/*********************************************************************************************************************************/

double HeatEquationFBVP::final(const SpaceNodePDE &sn, FinalCondition condition) const
{
    C_UNUSED(sn);
    C_UNUSED(condition);
    auto &u = solver->u;
    auto &V = solver->U;
    unsigned int i = static_cast<unsigned int>(sn.i);
    return -2.0*(u[i] - V[i]);
}

double HeatEquationFBVP::boundary(const SpaceNodePDE &sn, const TimeNodePDE &tn, BoundaryConditionPDE &condition) const
{
    C_UNUSED(sn);
    C_UNUSED(tn);
    C_UNUSED(condition);

    if (sn.i == 0) {
        const EquationParameters &params = solver->params;
        condition = BoundaryConditionPDE(BoundaryCondition::Robin, params.thermalConductivity1, -1.0, 0.0);
        return 0.0;
    } else {
        const EquationParameters &params = solver->params;
        condition = BoundaryConditionPDE(BoundaryCondition::Robin, params.thermalConductivity2, +1.0, 0.0);
        return 0.0;
    }
}

double HeatEquationFBVP::f(const SpaceNodePDE &sn, const TimeNodePDE &tn) const
{
    C_UNUSED(sn);
    C_UNUSED(tn);
    return 0.0;
}

void HeatEquationFBVP::layerInfo(const DoubleVector &p, const TimeNodePDE &tn) const
{
    C_UNUSED(p);
    C_UNUSED(tn);

    auto const_solver = const_cast<ProblemSolver *>(solver);
    const_solver->p0[tn.i] = p[0];
}

/*********************************************************************************************************************************/


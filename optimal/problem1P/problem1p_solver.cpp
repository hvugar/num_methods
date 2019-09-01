#include "problem1p_solver.h"

using namespace p1p;

void ProblemSolver::Main(int argc, char* argv[])
{
    C_UNUSED(argc);
    C_UNUSED(argv);

    unsigned int length = 101;
    ProblemSolver solver(Dimension(0.01, 0, 100), Dimension(0.01, 0, 100));

    DoubleVector x(length, 3.0);
    DoubleVector ng(length, 0.0);
    DoubleVector ag(length, 0.0);

    solver.gradient(x, ag);
    IGradient::Gradient(&solver, 0.0001, x, ng);
    IPrinter::printVector(10, 6, ag.EuclideanNormalize());
    IPrinter::printVector(10, 6, ng.EuclideanNormalize());
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

    uint32_t length = timeDimension.size();
    p0.resize(length);
    params.v =  new double [length];
    printf("_timeDimension.size: %d length: %d\n", _timeDimension.size(), length);
}

void ProblemSolver::setSpaceDimensionX(const Dimension &spaceDimensionX)
{
    this->_spaceDimensionX = spaceDimensionX;
    forward.setSpaceDimensionX(spaceDimensionX);
    backward.setSpaceDimensionX(spaceDimensionX);

    auto length = spaceDimensionX.size();
    U.resize(length, 0.0);
    V.resize(length, 10.0);
}

void ProblemSolver::gradient(const DoubleVector &x, DoubleVector &g) const
{
    C_UNUSED(x);
    C_UNUSED(g);

    const unsigned int length = x.length();
    for (unsigned int i=0; i<length; i++) const_this->params.v[i] = x[i];
    g.resize(length);
    printf("gradient length %d\n", length);

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
    auto N  = _spaceDimensionX.size() - 1;
    forward.implicit_calculate_D1V1CN();

    double sum = 0.0;
    sum += 0.5*(U[0]-V[0])*(U[0]-V[0]);
    for (unsigned int n=1; n<=N-1; n++)
    {
        sum += (U[n]-V[n])*(U[n]-V[n]);
    }
    sum += 0.5*(U[N]-V[N])*(U[N]-V[N]);
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
    if (_timeDimension.size() == static_cast<unsigned int>(tn.i))
        const_cast<ProblemSolver*>(solver)->U = u;
}

/*********************************************************************************************************************************/

double HeatEquationFBVP::final(const SpaceNodePDE &sn, FinalCondition condition) const
{
    C_UNUSED(sn);
    C_UNUSED(condition);
    auto &U = solver->U;
    auto &V = solver->V;
    unsigned int i = static_cast<unsigned int>(sn.i);
    return -2.0*(U[i] - V[i]);
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


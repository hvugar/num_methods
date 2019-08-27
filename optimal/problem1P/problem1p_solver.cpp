#include "problem1p_solver.h"

using namespace p1p;

void ProblemSolver::Main(int argc, char* argv[])
{
    C_UNUSED(argc);
    C_UNUSED(argv);




    HeatEquationIBVP hibvp;
    hibvp.setTimeDimension(Dimension(0.01, 0, 300));
    hibvp.setSpaceDimensionX(Dimension(0.01, 0, 100));

    EquationParameters param;
    param.v =  new double[301];
    param.initialTemperature = 0.0;
    param.environmentTemperature = 0.5;

    param.thermalConductivity0 = 0.0;
    param.thermalConductivity1 = 0.1;
    param.thermalConductivity2 = 0.01;

    for (unsigned int i=0; i<=300; i++) param.v[i] = 5.0;

    hibvp.params = &param;
    hibvp.implicit_calculate_D1V1CN();
}

ProblemSolver::ProblemSolver(const Dimension &timeDimension, const Dimension &spaceDimensionX)
{
    forward.setTimeDimension(timeDimension);
    forward.setSpaceDimensionX(spaceDimensionX);
    forward.setThermalDiffusivity(1.0);
    forward.setThermalConductivity(0.0);
    forward.setThermalConvection(0.0);

    backward.setTimeDimension(timeDimension);
    backward.setSpaceDimensionX(spaceDimensionX);
    backward.setThermalDiffusivity(-1.0);
    backward.setThermalConductivity(0.0);
    backward.setThermalConvection(0.0);
}

void ProblemSolver::setTimeDimension(const Dimension &timeDimension)
{
    this->_timeDimension = timeDimension;
    forward.setTimeDimension(timeDimension);
    backward.setTimeDimension(timeDimension);
}

void ProblemSolver::setSpaceDimensionX(const Dimension &spaceDimensionX)
{
    this->_spaceDimensionX = spaceDimensionX;
    forward.setSpaceDimensionX(spaceDimensionX);
    backward.setSpaceDimensionX(spaceDimensionX);
}

void ProblemSolver::gradient(const DoubleVector &x, DoubleVector &g) const
{
    ProblemSolver *const_solver = const_cast<ProblemSolver *>(this);

    if (const_solver->params.v != nullptr) const_solver->params.v;


    const_solver->forward.implicit_calculate_D1V1CN();
    const_solver->backward.implicit_calculate_D1V1CN();
}

double ProblemSolver::fx(const DoubleVector &x) const
{
    ProblemSolver *const_solver = const_cast<ProblemSolver *>(this);
    const_solver->params.v = const_cast<double *>(x.data());

    forward.implicit_calculate_D1V1CN();

    double sum = 0.0;
    sum += 0.5*(u[M][0]-V[0])*(u[M][0]-V[0]);
    for (unsigned int n=1; n<=N-1; n++)
    {
        sum += (u[M][n]-V[n])*(u[M][n]-V[n]);
    }
    sum += 0.5*(u[M][N]-V[N])*(u[M][N]-V[N]);
    sum *= hx;
    return sum;

    return 0.0;
}

/*********************************************************************************************************************************/

double HeatEquationIBVP::initial(const SpaceNodePDE &sn, InitialCondition condition) const
{
    C_UNUSED(sn);
    C_UNUSED(condition);

    const EquationParameters &params = solver->params;
    return params.initialTemperature;
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
    if (_timeDimension.size() == static_cast<int>(tn.i)) {
        const_cast<ProblemSolver*>(solver)->U = u;
    }
}

/*********************************************************************************************************************************/

double HeatEquationFBVP::final(const SpaceNodePDE &sn, FinalCondition condition) const
{
    C_UNUSED(sn);
    C_UNUSED(condition);

    DoubleVector &U = const_cast<ProblemSolver*>(solver)->U;
    DoubleVector &V = const_cast<ProblemSolver*>(solver)->V;
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

    return 0.0;
}

double HeatEquationFBVP::f(const SpaceNodePDE &sn, const TimeNodePDE &tn) const
{
    C_UNUSED(sn);
    C_UNUSED(tn);
    return 0.0;
}

void HeatEquationFBVP::layerInfo(const DoubleVector &u, const TimeNodePDE &tn) const
{
    C_UNUSED(u);
    C_UNUSED(tn);
}

/*********************************************************************************************************************************/


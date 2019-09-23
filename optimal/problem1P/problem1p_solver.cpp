#include "problem1p_solver.h"

#define EXAMPLE_LEFT_BORDER_ROBIN
//#define EXAMPLE_LEFT_BORDER_DIRICHLET
//#define EXAMPLE_FXT

using namespace p1p;

void ProblemSolver::Main(int argc, char* argv[])
{
    C_UNUSED(argc);
    C_UNUSED(argv);

    unsigned int length = 101;
    ProblemSolver solver(Dimension(0.01, 0, 100), Dimension(0.01, 0, 100));

#if defined (EXAMPLE_LEFT_BORDER_ROBIN) || defined(EXAMPLE_LEFT_BORDER_DIRICHLET)
    DoubleVector x(length, 3.0); for (unsigned int i=0; i<length; i++) x[i] = i*0.01;
#endif
    IPrinter::printVector(x);

    DoubleVector ng(length, 0.0);
    DoubleVector ag(length, 0.0);

    solver.gradient(x, ag);
    IGradient::Gradient(&solver, 0.01, x, ng);
    //IPrinter::printVector(ag);
    //IPrinter::printVector(ng);
    //ag[100] = ng[100] = 0.0;
    IPrinter::printVector(ag.EuclideanNormalize());
    IPrinter::printVector(ng.EuclideanNormalize());
}

ProblemSolver::ProblemSolver(const Dimension &timeDimension, const Dimension &spaceDimensionX)
{
    setTimeDimension(timeDimension);
    setSpaceDimensionX(spaceDimensionX);

    params.initialTemperature = 0.0;
    params.environmentTemperature = 0.5;
    params.thermalDiffusivity = 1.0;
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
}

void ProblemSolver::setSpaceDimensionX(const Dimension &spaceDimensionX)
{
    this->_spaceDimensionX = spaceDimensionX;
    forward.setSpaceDimensionX(spaceDimensionX);
    backward.setSpaceDimensionX(spaceDimensionX);

    auto length = spaceDimensionX.size();
    U.resize(length, 0.0);
    V.resize(length, 2.0);
}

void ProblemSolver::gradient(const DoubleVector &x, DoubleVector &g) const
{
    const unsigned int length = x.length();
    for (unsigned int i=0; i<length; i++) const_this->params.v[i] = x[i];

    forward.implicit_calculate_D1V1CN();
    backward.implicit_calculate_D1V1CN();

    g.resize(length);
#ifdef EXAMPLE_LEFT_BORDER_ROBIN
    for (unsigned int i=0; i<length; i++) g[i] = -params.thermalDiffusivity*params.thermalConductivity1*p0[i];
#endif
#ifdef EXAMPLE_LEFT_BORDER_DIRICHLET
    for (unsigned int i=0; i<length; i++) g[i] = -params.thermalDiffusivity*p0[i];
#endif
    g[0] = 0.0;
}

double ProblemSolver::fx(const DoubleVector &x) const
{
    const unsigned int length = x.length();
    for (unsigned int i=0; i<length; i++) const_this->params.v[i] = x[i];

    forward.implicit_calculate_D1V1CN();

    auto hx = _spaceDimensionX.step();
    auto N  = _spaceDimensionX.size() - 1;
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

    const EquationParameters &params = solver->params;
    if (sn.i == 0)
    {
#ifdef EXAMPLE_LEFT_BORDER_ROBIN
        condition = BoundaryConditionPDE(BoundaryCondition::Robin, params.thermalConductivity1, -1.0, params.thermalConductivity1);
#endif
#ifdef EXAMPLE_LEFT_BORDER_DIRICHLET
        condition = BoundaryConditionPDE(BoundaryCondition::Dirichlet, +1.0, +0.0, +1.0);
#endif
        return params.v[tn.i];
    }
    else
    {
#ifdef EXAMPLE_LEFT_BORDER_ROBIN
        condition = BoundaryConditionPDE(BoundaryCondition::Robin, params.thermalConductivity2, +1.0, params.thermalConductivity2);
        return params.environmentTemperature;
#endif
#ifdef EXAMPLE_LEFT_BORDER_DIRICHLET
        condition = BoundaryConditionPDE(BoundaryCondition::Dirichlet, +1.0, +0.0, +0.0);
        return 0.0;
#endif
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
    if (_timeDimension.size()-1 == static_cast<unsigned int>(tn.i))
    {
        const_cast<ProblemSolver*>(solver)->U = u;
    }
    //printf("%4d --- ", tn.i);
    //IPrinter::printVector(u);
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
#ifdef EXAMPLE_LEFT_BORDER_ROBIN
        const EquationParameters &params = solver->params;
        condition = BoundaryConditionPDE(BoundaryCondition::Robin, params.thermalConductivity1, -1.0, 0.0);
        return 0.0;
#endif
#ifdef EXAMPLE_LEFT_BORDER_DIRICHLET
        condition = BoundaryConditionPDE(BoundaryCondition::Dirichlet, +1.0, +0.0, +0.0);
        return 0.0;
#endif
    } else {
#ifdef EXAMPLE_LEFT_BORDER_ROBIN
        const EquationParameters &params = solver->params;
        condition = BoundaryConditionPDE(BoundaryCondition::Robin, params.thermalConductivity2, +1.0, 0.0);
        return 0.0;
#endif
#ifdef EXAMPLE_LEFT_BORDER_DIRICHLET
        condition = BoundaryConditionPDE(BoundaryCondition::Dirichlet, +1.0, +0.0, +0.0);
        return 0.0;
#endif
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
#ifdef EXAMPLE_LEFT_BORDER_ROBIN
    const_solver->p0[tn.i] = p[0];
#endif
#ifdef EXAMPLE_LEFT_BORDER_DIRICHLET
    const double hx = _spaceDimensionX.step();
    //const_solver->p0[tn.i] = (-3.0*p[0]+4.0*p[1]-p[2])/(2.0*hx);
    const_solver->p0[tn.i] = (p[1]-p[0])/hx;
#endif

    //printf("%4d --- ", tn.i);
    //IPrinter::printVector(p);

    //IPrinter::printVector(p);
    //printf("---- %4d %20.14f %20.14f %20.14f %20.14f %20.14f\n", tn.i, const_solver->params.v[tn.i], p[0], p[1], p[2], const_solver->p0[tn.i]);
}

/*********************************************************************************************************************************/


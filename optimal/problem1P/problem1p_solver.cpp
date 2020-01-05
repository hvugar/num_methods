#include "problem1p_solver.h"

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

    solver.forward.showLayers = true;
    solver.backward.showLayers = true;
    solver.gradient(x, ag);
    solver.forward.showLayers = false;
    solver.backward.showLayers = false;
    IGradient::Gradient(&solver, 0.01, x, ng);
    IPrinter::printVector(ag);
    IPrinter::printVector(ng);
    IPrinter::printSeperatorLine();
    ag[0]   = ng[0]   = 0.0;
    ag[100] = ng[100] = 0.0;
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
    params.thermalConductivity0 = 0.01;
    params.thermalConductivity1 = 0.1;
    params.thermalConductivity2 = 0.001;

    forward.setThermalDiffusivity(params.thermalDiffusivity);
    forward.setThermalConductivity(0.0);
    forward.setThermalConvection(-params.thermalConductivity0);

    backward.setThermalDiffusivity(-params.thermalDiffusivity);
    backward.setThermalConductivity(0.0);
    backward.setThermalConvection(params.thermalConductivity0);

    forward.solver = this;
    backward.solver = this;
    const_this = const_cast<ProblemSolver*>(this);
}

void ProblemSolver::setTimeDimension(const Dimension &timeDimension)
{
    this->_timeDimension = timeDimension;
    forward.setTimeDimension(timeDimension);
    backward.setTimeDimension(timeDimension);

    auto length = timeDimension.size();
    p0.resize(length);
    p0x.resize(length);
    params.v =  new double [length];
}

void ProblemSolver::setSpaceDimensionX(const Dimension &spaceDimensionX)
{
    this->_spaceDimensionX = spaceDimensionX;
    forward.setSpaceDimensionX(spaceDimensionX);
    backward.setSpaceDimensionX(spaceDimensionX);

    auto length = spaceDimensionX.size();
    U.resize(length, 0.0);
    V.resize(length, 0.5);
}

void ProblemSolver::gradient(const DoubleVector &x, DoubleVector &g) const
{
    const unsigned int length = x.length();
    for (unsigned int i=0; i<length; i++) const_this->params.v[i] = x[i];

    forward.implicit_calculate_D1V1();
    backward.implicit_calculate_D1V1();

    g.resize(length);
#ifdef EXAMPLE_LEFT_BORDER_ROBIN
    for (unsigned int i=0; i<length; i++) g[i] = -params.thermalDiffusivity*params.thermalConductivity1*p0[i];
#endif
#ifdef EXAMPLE_LEFT_BORDER_DIRICHLET
    for (unsigned int i=0; i<length; i++) g[i] = -params.thermalDiffusivity*p0x[i];
#endif
    g[0] = 0.0;
}

double ProblemSolver::fx(const DoubleVector &x) const
{
    const unsigned int length = x.length();
    for (unsigned int i=0; i<length; i++) const_this->params.v[i] = x[i];

    forward.implicit_calculate_D1V1();

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
{ return solver->frw_initial(sn, condition); }

double HeatEquationIBVP::boundary(const SpaceNodePDE &sn, const TimeNodePDE &tn, BoundaryConditionPDE &condition) const
{ return solver->frw_boundary(sn, tn, condition); }

double HeatEquationIBVP::f(const SpaceNodePDE &sn, const TimeNodePDE &tn) const
{ return solver->frw_f(sn, tn); }

void HeatEquationIBVP::layerInfo(const DoubleVector &u, const TimeNodePDE &tn) const
{ return solver->frw_layerInfo(u, tn); }

double ProblemSolver::frw_initial(const SpaceNodePDE &sn, InitialCondition condition) const
{
    C_UNUSED(sn);
    C_UNUSED(condition);
    return params.initialTemperature;
}

double ProblemSolver::frw_boundary(const SpaceNodePDE &sn, const TimeNodePDE &tn, BoundaryConditionPDE &condition) const
{
    C_UNUSED(sn);
    C_UNUSED(tn);
    C_UNUSED(condition);

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
        condition = BoundaryConditionPDE(BoundaryCondition::Dirichlet, +1.0, +0.0, +1.0);
        return 0.0;
#endif
    }
}

double ProblemSolver::frw_f(const SpaceNodePDE &sn, const TimeNodePDE &tn) const
{
    C_UNUSED(sn);
    C_UNUSED(tn);
    return params.thermalConductivity0 * params.environmentTemperature;
}

void ProblemSolver::frw_layerInfo(const DoubleVector &u, const TimeNodePDE &tn) const
{
    C_UNUSED(u);
    C_UNUSED(tn);
    //Dimension _timeDimension = timeDimension();
    if (_timeDimension.size()-1 == static_cast<unsigned int>(tn.i))
    {
        //const_cast<ProblemSolver*>(solver)->U = u;
        const_cast<ProblemSolver*>(this)->U = u;
    }

//    if (showLayers)
//    {
//        printf("%4d --- ", tn.i);
//        IPrinter::printVector(u);
//    }
}

/*********************************************************************************************************************************/

double HeatEquationFBVP::final(const SpaceNodePDE &sn, FinalCondition condition) const
{ return solver->bcw_final(sn, condition); }

double HeatEquationFBVP::boundary(const SpaceNodePDE &sn, const TimeNodePDE &tn, BoundaryConditionPDE &condition) const
{ return solver->bcw_boundary(sn, tn, condition); }

double HeatEquationFBVP::f(const SpaceNodePDE &sn, const TimeNodePDE &tn) const
{ return solver->bcw_f(sn, tn); }

void HeatEquationFBVP::layerInfo(const DoubleVector &p, const TimeNodePDE &tn) const
{ return solver->bcw_layerInfo(p, tn); }

double ProblemSolver::bcw_final(const SpaceNodePDE &sn, FinalCondition condition) const
{
    C_UNUSED(sn);
    C_UNUSED(condition);
    unsigned int i = static_cast<unsigned int>(sn.i);
    return -2.0*(U[i] - V[i]);
}

double ProblemSolver::bcw_boundary(const SpaceNodePDE &sn, const TimeNodePDE &tn, BoundaryConditionPDE &condition) const
{
    C_UNUSED(sn);
    C_UNUSED(tn);
    C_UNUSED(condition);

    if (sn.i == 0) {
#ifdef EXAMPLE_LEFT_BORDER_ROBIN
        //const EquationParameters &params = solver->params;
        condition = BoundaryConditionPDE(BoundaryCondition::Robin, params.thermalConductivity1, -1.0, 0.0);
        return 0.0;
#endif
#ifdef EXAMPLE_LEFT_BORDER_DIRICHLET
        condition = BoundaryConditionPDE(BoundaryCondition::Dirichlet, +1.0, +0.0, +1.0);
        return 0.0;
#endif
    } else {
#ifdef EXAMPLE_LEFT_BORDER_ROBIN
        //const EquationParameters &params = solver->params;
        condition = BoundaryConditionPDE(BoundaryCondition::Robin, params.thermalConductivity2, +1.0, 0.0);
        return 0.0;
#endif
#ifdef EXAMPLE_LEFT_BORDER_DIRICHLET
        condition = BoundaryConditionPDE(BoundaryCondition::Dirichlet, +1.0, +0.0, +1.0);
        return 0.0;
#endif
    }
}

double ProblemSolver::bcw_f(const SpaceNodePDE &sn, const TimeNodePDE &tn) const
{
    C_UNUSED(sn);
    C_UNUSED(tn);
    return 0.0;
}

void ProblemSolver::bcw_layerInfo(const DoubleVector &p, const TimeNodePDE &tn) const
{
    C_UNUSED(p);
    C_UNUSED(tn);

    auto const_solver = const_cast<ProblemSolver *>(this);
#ifdef EXAMPLE_LEFT_BORDER_ROBIN
    const_solver->p0[tn.i] = p[0];
#endif
#ifdef EXAMPLE_LEFT_BORDER_DIRICHLET
    const double hx = _spaceDimensionX.step();
    const_solver->p0x[tn.i] = (-3.0*p[0]+4.0*p[1]-p[2])/(2.0*hx);
    //const_solver->p0x[tn.i] = (p[1]-p[0])/hx;
#endif

//    if (showLayers)
//    {
//        printf("%4d --- ", tn.i);
//        IPrinter::printVector(p);
//    }


    //printf("%4d --- ", tn.i);
    //IPrinter::printVector(p);

    //IPrinter::printVector(p);
    //printf("---- %4d %20.14f %20.14f %20.14f %20.14f %20.14f\n", tn.i, const_solver->params.v[tn.i], p[0], p[1], p[2], const_solver->p0[tn.i]);
}

/*********************************************************************************************************************************/


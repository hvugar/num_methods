#include "problem1p_solver.h"

using namespace p1p;

void ProblemSolver::Main(int argc, char* argv[])
{
    C_UNUSED(argc);
    C_UNUSED(argv);

    unsigned int length = 101;
    ProblemSolver solver(Dimension(0.01, 0, 100), Dimension(0.01, 0, 100));

#if defined (EXAMPLE_LEFT_BORDER_ROBIN) || defined(EXAMPLE_LEFT_BORDER_DIRICHLET)
    DoubleVector x(length, 3.0);// for (unsigned int i=0; i<length; i++) x[i] = (i+1)*(1.0);
#endif

    IPrinter::printSeperatorLine("x");
    IPrinter::printVector(x);

    DoubleVector ng(length, 0.0);
    DoubleVector ag(length, 0.0);

    solver.gradient(x, ag);
    IGradient::Gradient(&solver, 0.01, x, ng);

    ag[0]   = ng[0]   = 0.0;
    //    ag[100] = ng[100] = 0.0;

    IPrinter::printSeperatorLine("gradients", '=');
    IPrinter::printVector(ag, "ag: ");
    IPrinter::printVector(ng, "ng: ");
    IPrinter::printSeperatorLine("normolized", '=');
    IPrinter::printVector(ag.EuclideanNormalize(), "ag: ");
    IPrinter::printVector(ng.EuclideanNormalize(), "ng: ");
    IPrinter::printSeperatorLine(nullptr, '=');
}

ProblemSolver::ProblemSolver(const Dimension &timeDimension, const Dimension &spaceDimensionX)
{
    forward.solver = backward.solver = this;

    setTimeDimension(timeDimension);
    setSpaceDimensionX(spaceDimensionX);

    initialTemperature = 0.0;
    environmentTemperature = 0.5;

    forward.setThermalDiffusivity(thermalDiffusivity);
    forward.setThermalConvection(-thermalConvection);
    forward.setThermalConductivity(0.0);

    backward.setThermalDiffusivity(-thermalDiffusivity);
    backward.setThermalConvection(thermalConvection);
    backward.setThermalConductivity(0.0);

    const_this = const_cast<ProblemSolver*>(this);
}

void ProblemSolver::setTimeDimension(const Dimension &timeDimension)
{
    _timeDimension = timeDimension;

    unsigned int size = timeDimension.size();

    if (heat_power != nullptr) delete [] heat_power;
    heat_power = new double[size];

    p0.resize(size);
    p1.resize(size);
    p2.resize(size);
    p0x.resize(size);
}

void ProblemSolver::setSpaceDimensionX(const Dimension &spaceDimensionX)
{
    _spaceDimensionX = spaceDimensionX;

    unsigned int size = spaceDimensionX.size();
    U.resize(size, 0.0);
    V.resize(size, 0.5);
}

void ProblemSolver::gradient(const DoubleVector &x, DoubleVector &g) const
{
    const unsigned int length = x.length();
    for (unsigned int i=0; i<length; i++) const_this->heat_power[i] = x[i];

    forward.implicit_calculate_D1V1();
    backward.implicit_calculate_D1V1();

    g.resize(length);
#ifdef EXAMPLE_LEFT_BORDER_ROBIN
    for (unsigned int i=0; i<length; i++) g[i] = -thermalDiffusivity*p0[i];
#endif
#ifdef EXAMPLE_LEFT_BORDER_DIRICHLET
    for (unsigned int i=0; i<length; i++) g[i] = -thermalDiffusivity*p0x[i];
#endif
    //    g[0] = 0.0;
    //    g[length-1] = 0.0;
}

double ProblemSolver::fx(const DoubleVector &x) const
{
    double sum = integral(x);
    return sum;
}

double ProblemSolver::integral(const DoubleVector &x) const
{
    const unsigned int length = x.length();
    const double hx = _spaceDimensionX.step();
    const unsigned int N = _spaceDimensionX.size() - 1;

    for (unsigned int i=0; i<length; i++) const_this->heat_power[i] = x[i];
    forward.implicit_calculate_D1V1();

    double integ_sum = 0.0;
    integ_sum += 0.5*(U[0]-V[0])*(U[0]-V[0]);
    for (unsigned int n=1; n<=N-1; n++)
    {
        integ_sum += (U[n]-V[n])*(U[n]-V[n]);
    }
    integ_sum += 0.5*(U[N]-V[N])*(U[N]-V[N]);
    integ_sum *= hx;

    return integ_sum;
}

/*********************************************************************************************************************************/

double ProblemSolver::frw_initial(const SpaceNodePDE &, InitialCondition) const { return initialTemperature; }

double ProblemSolver::frw_boundary(const SpaceNodePDE &sn, const TimeNodePDE &tn, BoundaryConditionPDE &c) const
{
    if (sn.i == 0)
    {
#ifdef EXAMPLE_LEFT_BORDER_ROBIN
        c = BoundaryConditionPDE(BoundaryCondition::Robin, thermalConductivity1, -1.0, thermalConductivity1);
#endif
#ifdef EXAMPLE_LEFT_BORDER_DIRICHLET
        c = BoundaryConditionPDE(BoundaryCondition::Dirichlet, +1.0, +0.0, +1.0);
#endif
        return heat_power[tn.i];
    }
    else
    {
#ifdef EXAMPLE_LEFT_BORDER_ROBIN
        c = BoundaryConditionPDE(BoundaryCondition::Robin, thermalConductivity2, +1.0, thermalConductivity2);
        return environmentTemperature;
#endif
#ifdef EXAMPLE_LEFT_BORDER_DIRICHLET
        c = BoundaryConditionPDE(BoundaryCondition::Dirichlet, +1.0, +0.0, +1.0);
        return 0.0;
#endif
    }
}

double ProblemSolver::frw_f(const SpaceNodePDE &, const TimeNodePDE &) const { return thermalConvection * environmentTemperature; }

void ProblemSolver::frw_layerInfo(const DoubleVector &u, const TimeNodePDE &tn) const
{
    if (static_cast<int>(tn.i) == timeDimension().max()) const_cast<ProblemSolver*>(this)->U = u;
}

/*********************************************************************************************************************************/

double ProblemSolver::bcw_final(const SpaceNodePDE &sn, FinalCondition) const
{
    unsigned int i = static_cast<unsigned int>(sn.i);
    return -2.0*(U[i] - V[i]);
}

double ProblemSolver::bcw_boundary(const SpaceNodePDE &sn, const TimeNodePDE &, BoundaryConditionPDE &c) const
{
    if (sn.i == 0) {
#ifdef EXAMPLE_LEFT_BORDER_ROBIN
        c = BoundaryConditionPDE(BoundaryCondition::Robin, thermalConductivity1, -1.0, 0.0);
        return 0.0;
#endif
#ifdef EXAMPLE_LEFT_BORDER_DIRICHLET
        c = BoundaryConditionPDE(BoundaryCondition::Dirichlet, +1.0, +0.0, +1.0);
        return 0.0;
#endif
    } else {
#ifdef EXAMPLE_LEFT_BORDER_ROBIN
        c = BoundaryConditionPDE(BoundaryCondition::Robin, thermalConductivity2, +1.0, 0.0);
        return 0.0;
#endif
#ifdef EXAMPLE_LEFT_BORDER_DIRICHLET
        c = BoundaryConditionPDE(BoundaryCondition::Dirichlet, +1.0, +0.0, +1.0);
        return 0.0;
#endif
    }
}

double ProblemSolver::bcw_f(const SpaceNodePDE &, const TimeNodePDE &) const { return 0.0; }

void ProblemSolver::bcw_layerInfo(const DoubleVector &p, const TimeNodePDE &tn) const
{
    C_UNUSED(p);
    C_UNUSED(tn);

    auto const_solver = const_cast<ProblemSolver *>(this);
#ifdef EXAMPLE_LEFT_BORDER_ROBIN
    if (tn.i==timeDimension().max())
    {
        const_solver->p0[tn.i] = 0.5*p[0];
    }
    else
    {
        const_solver->p0[tn.i+1] += 0.5*p[0];
        const_solver->p0[tn.i+0]  = 0.5*p[0];
    }
    //const_solver->p0[tn.i] = p[0];
#endif
#ifdef EXAMPLE_LEFT_BORDER_DIRICHLET
    const double hx = _spaceDimensionX.step();

    if (tn.i==timeDimension().max())
    {
        const_solver->p0[tn.i] = 0.5*p[0];
        const_solver->p1[tn.i] = 0.5*p[1];
        const_solver->p2[tn.i] = 0.5*p[2];
    }
    else
    {
        const_solver->p0[tn.i+1] += 0.5*p[0];
        const_solver->p1[tn.i+1] += 0.5*p[1];
        const_solver->p2[tn.i+1] += 0.5*p[2];

        const_solver->p0[tn.i] = 0.5*p[0];
        const_solver->p1[tn.i] = 0.5*p[1];
        const_solver->p2[tn.i] = 0.5*p[2];

        const_solver->p0x[tn.i+1] = (-3.0*p0[tn.i+1]+4.0*p1[tn.i+1]-p2[tn.i+1])/(2.0*hx);
    }
    //const_solver->p0x[tn.i] = (-3.0*p[0]+4.0*p[1]-p[2])/(2.0*hx);
#endif
}

/*********************************************************************************************************************************/

double HeatEquationIBVP::initial(const SpaceNodePDE &sn, InitialCondition c) const { return solver->frw_initial(sn, c); }

double HeatEquationIBVP::boundary(const SpaceNodePDE &sn, const TimeNodePDE &tn, BoundaryConditionPDE &c) const { return solver->frw_boundary(sn, tn, c); }

double HeatEquationIBVP::f(const SpaceNodePDE &sn, const TimeNodePDE &tn) const { return solver->frw_f(sn, tn); }

void HeatEquationIBVP::layerInfo(const DoubleVector &u, const TimeNodePDE &tn) const { return solver->frw_layerInfo(u, tn); }

Dimension HeatEquationIBVP::timeDimension() const { return solver->timeDimension(); }

Dimension HeatEquationIBVP::spaceDimensionX() const { return solver->spaceDimensionX(); }

Dimension HeatEquationIBVP::spaceDimensionY() const { return solver->spaceDimensionX(); }

Dimension HeatEquationIBVP::spaceDimensionZ() const { return solver->spaceDimensionX(); }

/*********************************************************************************************************************************/

double HeatEquationFBVP::final(const SpaceNodePDE &sn, FinalCondition c) const { return solver->bcw_final(sn, c); }

double HeatEquationFBVP::boundary(const SpaceNodePDE &sn, const TimeNodePDE &tn, BoundaryConditionPDE &c) const { return solver->bcw_boundary(sn, tn, c); }

double HeatEquationFBVP::f(const SpaceNodePDE &sn, const TimeNodePDE &tn) const { return solver->bcw_f(sn, tn); }

void HeatEquationFBVP::layerInfo(const DoubleVector &p, const TimeNodePDE &tn) const { return solver->bcw_layerInfo(p, tn); }

Dimension HeatEquationFBVP::timeDimension() const { return solver->timeDimension(); }

Dimension HeatEquationFBVP::spaceDimensionX() const { return solver->spaceDimensionX(); }

Dimension HeatEquationFBVP::spaceDimensionY() const { return solver->spaceDimensionX(); }

Dimension HeatEquationFBVP::spaceDimensionZ() const { return solver->spaceDimensionX(); }


#include "problem1h_solver.h"

using namespace h1p;

void ProblemSolver::Main(int argc, char* argv[])
{
    C_UNUSED(argc);
    C_UNUSED(argv);

    unsigned int length = 101;
    ProblemSolver solver(Dimension(0.01, 0, 100), Dimension(0.01, 0, 100));

#if defined (_LEFT_BORDER_DIRICHLET) || defined(_LEFT_BORDER_ROBIN)
    DoubleVector x(length, 0.0); for (unsigned int i=0; i<length; i++) x[i] = sin(i*1.01);
#endif
    IPrinter::printVector(x);
    IPrinter::printSeperatorLine();

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

    forward.setWaveSpeed(params.waveSpeed);
    forward.setWaveDissipation(params.waveDissipation);

    backward.setWaveSpeed(params.waveSpeed);
    backward.setWaveDissipation(-params.waveDissipation);

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

    params.v0 =  new double [length]; for (unsigned int i=0; i<length; i++) params.v0[i] = 0.0;
    params.v1 =  new double [length]; for (unsigned int i=0; i<length; i++) params.v1[i] = 0.0;
}

void ProblemSolver::setSpaceDimensionX(const Dimension &spaceDimensionX)
{
    this->_spaceDimensionX = spaceDimensionX;
    forward.setSpaceDimensionX(spaceDimensionX);
    backward.setSpaceDimensionX(spaceDimensionX);

    auto length = spaceDimensionX.size();
    U1.resize(length, 0.0);
    V1.resize(length, 0.5);
    U2.resize(length, 0.0);
    V2.resize(length, 0.3);
}

void ProblemSolver::gradient(const DoubleVector &x, DoubleVector &g) const
{
    const unsigned int length = x.length();
    for (unsigned int i=0; i<length; i++) const_this->params.v0[i] = x[i];

    forward.implicit_calculate_D1V1();
    backward.implicit_calculate_D1V1();

    g.resize(length);
#ifdef _LEFT_BORDER_ROBIN
    for (unsigned int i=0; i<length; i++) g[i] = -params.waveSpeed*params.waveSpeed*params.lambda1*p0[i];
#endif
#ifdef _LEFT_BORDER_DIRICHLET
    for (unsigned int i=0; i<length; i++) g[i] = -params.waveSpeed*params.waveSpeed*p0x[i];
#endif
    g[0] = 0.0;
}

double ProblemSolver::fx(const DoubleVector &x) const
{
    const unsigned int length = x.length();

    for (unsigned int i=0; i<length; i++) const_this->params.v0[i] = x[i];
    //for (unsigned int i=0; i<length/2; i++) const_this->params.v1[i] = x[i+length/2];

    forward.implicit_calculate_D1V1();

    auto hx = _spaceDimensionX.step();
    auto N  = _spaceDimensionX.size() - 1;

    double sum1 = 0.0;
    sum1 += 0.5*(U1[0]-V1[0])*(U1[0]-V1[0]);
    for (unsigned int n=1; n<=N-1; n++)
    {
        sum1 += (U1[n]-V1[n])*(U1[n]-V1[n]);
    }
    sum1 += 0.5*(U1[N]-V1[N])*(U1[N]-V1[N]);
    sum1 *= hx;

    double sum2 = 0.0;
    sum2 += 0.5*(U2[0]-V2[0])*(U2[0]-V2[0]);
    for (unsigned int n=1; n<=N-1; n++)
    {
        sum2 += (U2[n]-V2[n])*(U2[n]-V2[n]);
    }
    sum2 += 0.5*(U2[N]-V2[N])*(U2[N]-V2[N]);
    sum2 *= hx;

    return sum1+sum2;
}

/*********************************************************************************************************************************/

double WaveEquationIBVP::initial(const SpaceNodePDE &sn, InitialCondition condition) const
{ return solver->frw_initial(sn, condition); }

double WaveEquationIBVP::boundary(const SpaceNodePDE &sn, const TimeNodePDE &tn, BoundaryConditionPDE &condition) const
{ return solver->frw_boundary(sn, tn, condition); }

double WaveEquationIBVP::f(const SpaceNodePDE &sn, const TimeNodePDE &tn) const
{ return solver->frw_f(sn, tn); }

void WaveEquationIBVP::layerInfo(const DoubleVector &u, const TimeNodePDE &tn) const
{ return solver->frw_layerInfo(u, tn); }

double ProblemSolver::frw_initial(const SpaceNodePDE &sn, InitialCondition condition) const
{
    C_UNUSED(sn);
    C_UNUSED(condition);
    return 0.0;
}

double ProblemSolver::frw_boundary(const SpaceNodePDE &sn, const TimeNodePDE &tn, BoundaryConditionPDE &condition) const
{
    C_UNUSED(sn);
    C_UNUSED(tn);
    C_UNUSED(condition);

    if (sn.i == 0)
    {
#ifdef _LEFT_BORDER_DIRICHLET
        condition = BoundaryConditionPDE(BoundaryCondition::Dirichlet, +1.0, +0.0, +1.0);
#endif
#ifdef _LEFT_BORDER_ROBIN
        condition = BoundaryConditionPDE(BoundaryCondition::Robin, params.lambda1, -1.0, params.lambda1);
#endif
        return params.v0[tn.i];
    }
    else
    {
#ifdef _RIGHT_BORDER_DIRICHLET
        condition = BoundaryConditionPDE(BoundaryCondition::Dirichlet, +1.0, +0.0, +1.0);
#endif
#ifdef _RIGHT_BORDER_ROBIN
        condition = BoundaryConditionPDE(BoundaryCondition::Robin, params.lambda2, -1.0, params.lambda2);
#endif
        return params.v1[tn.i];
    }
}

double ProblemSolver::frw_f(const SpaceNodePDE &sn, const TimeNodePDE &tn) const
{
    C_UNUSED(sn);
    C_UNUSED(tn);
    return 0.0;
}

void ProblemSolver::frw_layerInfo(const DoubleVector &u, const TimeNodePDE &tn) const
{
    C_UNUSED(u);
    C_UNUSED(tn);
    if (_timeDimension.max() == tn.i)
    {
        const_cast<ProblemSolver*>(this)->U1 = u;
    }
    //IPrinter::printVector(u);
}

/*********************************************************************************************************************************/

double WaveEquationFBVP::final(const SpaceNodePDE &sn, FinalCondition condition) const
{ return solver->bcw_final(sn, condition); }

double WaveEquationFBVP::boundary(const SpaceNodePDE &sn, const TimeNodePDE &tn, BoundaryConditionPDE &condition) const
{ return solver->bcw_boundary(sn, tn, condition); }

double WaveEquationFBVP::f(const SpaceNodePDE &sn, const TimeNodePDE &tn) const
{ return solver->bcw_f(sn, tn); }

void WaveEquationFBVP::layerInfo(const DoubleVector &p, const TimeNodePDE &tn) const
{ return solver->bcw_layerInfo(p, tn); }

double ProblemSolver::bcw_final(const SpaceNodePDE &sn, FinalCondition condition) const
{
    unsigned int i = static_cast<unsigned int>(sn.i);
    if (condition == FinalCondition::FinalValue)
        return -2.0*(U2[i] - V2[i]);
    else
        return +2.0*(U1[i] - V1[i]) + params.waveDissipation*bcw_final(sn, FinalCondition::FinalValue);
}

double ProblemSolver::bcw_boundary(const SpaceNodePDE &sn, const TimeNodePDE &, BoundaryConditionPDE &condition) const
{
    if (sn.i == 0)
    {
#ifdef _LEFT_BORDER_DIRICHLET
        condition = BoundaryConditionPDE(BoundaryCondition::Dirichlet, +1.0, +0.0, +1.0);
#endif
#ifdef _LEFT_BORDER_ROBIN
        condition = BoundaryConditionPDE(BoundaryCondition::Robin, params.lambda1, -1.0, 0.0);
#endif
        return 0.0;
    }
    else
    {
#ifdef _RIGHT_BORDER_DIRICHLET
        condition = BoundaryConditionPDE(BoundaryCondition::Robin, params.lambda2, -1.0, 0.0);
#endif
#ifdef _RIGHT_BORDER_DIRICHLET
        condition = BoundaryConditionPDE(BoundaryCondition::Dirichlet, +1.0, +0.0, +1.0);
#endif
        return 0.0;
    }
}

double ProblemSolver::bcw_f(const SpaceNodePDE &, const TimeNodePDE &) const
{
    return 0.0;
}

void ProblemSolver::bcw_layerInfo(const DoubleVector &p, const TimeNodePDE &tn) const
{
    C_UNUSED(p);
    C_UNUSED(tn);

    auto const_solver = const_cast<ProblemSolver *>(this);
#ifdef _LEFT_BORDER_ROBIN
    const_solver->p0[tn.i] = p[0];
#endif
#ifdef _LEFT_BORDER_DIRICHLET
    const double hx = _spaceDimensionX.step();
    const_solver->p0x[tn.i] = (-3.0*p[0]+4.0*p[1]-p[2])/(2.0*hx);
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





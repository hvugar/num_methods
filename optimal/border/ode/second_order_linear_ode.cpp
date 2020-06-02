#include "second_order_linear_ode.h"

#define SO_TIME_STP 0.01
#define SO_TIME_MIN 0
#define SO_TIME_MAX 100

void SecondOrderLinearODEIBVP::Main(int argc, char **argv)
{
    C_UNUSED(argc);
    C_UNUSED(argv);
    CauchyProblemExample();
}

void SecondOrderLinearODEIBVP::CauchyProblemExample()
{
    SecondOrderLinearODEIBVP fnl;
    fnl.solveInitialValueProblem(ODESolverMethod::RUNGE_KUTTA_4);
    IPrinter::printSeperatorLine();
}

double SecondOrderLinearODEIBVP::A(const PointNodeODE &, size_t, size_t) const
{
    return -2.0;
}

double SecondOrderLinearODEIBVP::B(const PointNodeODE &, size_t, size_t) const
{
    return -3.0;
}

double SecondOrderLinearODEIBVP::C(const PointNodeODE &node, size_t r) const
{
    return 2.0 - A(node, 1, 1) * (2.0*node.x+1) - B(node, 1, 1) * (node.x*node.x+node.x);
}

size_t SecondOrderLinearODEIBVP::count() const
{
    return 1;
}

auto SecondOrderLinearODEIBVP::dimension() const -> Dimension { return Dimension(SO_TIME_STP, SO_TIME_MIN, SO_TIME_MAX); }

auto SecondOrderLinearODEIBVP::initial(InitialCondition c, size_t r) const -> double
{
    return c == InitialCondition::InitialValue ? 0.0 : 1.0;
}

void SecondOrderLinearODEIBVP::iterationInfo(const DoubleVector &v, const PointNodeODE &node) const
{
    if (node.i%((dimension().size()-1)/10)==0) { printf("%6d: ", node.i); IPrinter::print(v, v.length()); }
}

auto SecondOrderLinearODEIBVP::boundary(const PointNodeODE &, BoundaryConditionPDE &, size_t) const -> double
{
    throw runtime_error("");
}

/****************************************************************************************************************************/

void SecondOrderLinearODEFBVP::Main(int argc, char **argv)
{
    C_UNUSED(argc);
    C_UNUSED(argv);
    CauchyProblemExample();
}

void SecondOrderLinearODEFBVP::CauchyProblemExample()
{
    SecondOrderLinearODEFBVP fnl;
    fnl.solveInitialValueProblem2(ODESolverMethod::RUNGE_KUTTA_4);
    IPrinter::printSeperatorLine();
}

double SecondOrderLinearODEFBVP::A(const PointNodeODE &, size_t, size_t) const
{
    return +2.0;
}

double SecondOrderLinearODEFBVP::B(const PointNodeODE &, size_t, size_t) const
{
    return +3.0;
}

double SecondOrderLinearODEFBVP::C(const PointNodeODE &node, size_t r) const
{
    return 2.0 - A(node, 1, 1) * (2.0*node.x+1) - B(node, 1, 1) * (node.x*node.x+node.x);
}

size_t SecondOrderLinearODEFBVP::count() const
{
    return 1;
}

auto SecondOrderLinearODEFBVP::dimension() const -> Dimension { return Dimension(SO_TIME_STP, SO_TIME_MIN, SO_TIME_MAX); }

auto SecondOrderLinearODEFBVP::final(FinalCondition c, size_t r) const -> double
{
    return c == FinalCondition::FinalValue ? 2.0 : 3.0;
}

void SecondOrderLinearODEFBVP::iterationInfo(const DoubleVector &v, const PointNodeODE &node) const
{
    if (node.i%((dimension().size()-1)/10)==0) { printf("%6d: ", node.i); IPrinter::print(v, v.length()); }
}

auto SecondOrderLinearODEFBVP::boundary(const PointNodeODE &, BoundaryConditionPDE &, size_t) const -> double
{
    throw runtime_error("");
}


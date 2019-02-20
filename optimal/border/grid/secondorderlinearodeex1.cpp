#include "secondorderlinearodeex1.h"

void SecondOrderLinearODEEx1::Main(int argc, char **argv)
{
    SecondOrderLinearODEEx1 sol;
    DoubleVector rv;
    sol.solveInitialValueProblem(rv);
    IPrinter::printVector(rv);
}

SecondOrderLinearODEEx1::SecondOrderLinearODEEx1()
{}

SecondOrderLinearODEEx1::~SecondOrderLinearODEEx1()
{}

double SecondOrderLinearODEEx1::A(const PointNodeODE &node, unsigned int row, unsigned int col) const
{
    return 1.0;
}

double SecondOrderLinearODEEx1::B(const PointNodeODE &node, unsigned int row, unsigned int col) const
{
    return 2.0;
}

double SecondOrderLinearODEEx1::C(const PointNodeODE &node, unsigned int row) const
{
    double t = node.x;

    return 2.0 - 2.0*t*A(node) - t*t*B(node);
}

void SecondOrderLinearODEEx1::boundary(const PointNodeODE &node, BoundaryConditionODE &) const
{

}

void SecondOrderLinearODEEx1::initial(const PointNodeODE &node, InitialConditionODE &condition) const
{
    if (condition.initialConditionType == InitialConditionType::InitialValue)
        condition.value = 0.0;
    else
        condition.value = 0.0;
}

unsigned int SecondOrderLinearODEEx1::count() const
{
    return 1;
}

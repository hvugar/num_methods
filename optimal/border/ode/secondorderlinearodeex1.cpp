#include "secondorderlinearodeex1.h"

void SecondOrderLinearODEEx1::Main(int argc, char **argv)
{
    SecondOrderLinearODEEx1 sol;
    sol.setDimension(Dimension(0.1,0,10));
    DoubleVector rv;
    //sol.solveInitialValueProblem(rv);
    sol.solveBoundaryValueProblem(rv);
    IPrinter::printVector(rv);
}

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

double SecondOrderLinearODEEx1::boundary(const PointNodeODE &node, BoundaryConditionODE &condition, unsigned int row) const
{
    const Dimension &dim = dimension();
    if (node.i == dim.min())
    {
//        condition.boundaryConditionType = BoundaryConditionTypeODE::Dirichlet;
//        condition.lambda = 0.0;
//        return 0.0;
        condition.boundaryConditionType = BoundaryConditionTypeODE::Neumann;
        condition.lambda = 1.0;
        return 0.0;
    }
    else if (node.i == dim.max())
    {
//        condition.boundaryConditionType = BoundaryConditionTypeODE::Dirichlet;
//        condition.lambda = 0.0;
//        return 1.0;
        condition.boundaryConditionType = BoundaryConditionTypeODE::Neumann;
        condition.lambda = 1.0;
        return 2.0 - condition.lambda;
    }
    return NAN;
}

double SecondOrderLinearODEEx1::initial(InitialConditionTypeODE condition, unsigned int) const
{
    if (condition == InitialConditionTypeODE::InitialValue) return 0.0;
    if (condition == InitialConditionTypeODE::FirstDerivative) return 0.0;
    return NAN;
}

unsigned int SecondOrderLinearODEEx1::count() const
{
    return 1;
}

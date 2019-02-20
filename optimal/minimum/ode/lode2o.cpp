#include "lode2o.h"

void SecondOrderLinearODE::solveBoundaryValueProblem(DoubleVector &rv) const
{
    if (count() != 1) throw ExceptionODE(5);

    const Dimension &dim = dimension();
    int min = dim.min();
    int max = dim.max();
    double h = dim.step();
    unsigned int size = static_cast<unsigned int>(max-min);

    PointNodeODE leftNode(static_cast<double>(min*h), min);
    BoundaryConditionODE leftCondition;
    boundary(leftNode, leftCondition);
    double leftAlpha = leftCondition.a.at(0,0);
    double leftBetta = leftCondition.b.at(0,0);
    double leftGamma = leftCondition.c.at(0);

    PointNodeODE rightNode(static_cast<double>(max*h), max);
    BoundaryConditionODE rightCondition;
    boundary(rightNode, rightCondition);
    double rightAlpha = rightCondition.a.at(0,0);
    double rightBetta = rightCondition.b.at(0,0);
    double rightGamma = rightCondition.c.at(0);

    unsigned int N = size + 1;
    if (leftCondition.boundaryConditionType == BoundaryConditionType::Dirichlet) N--;
    if (rightCondition.boundaryConditionType == BoundaryConditionType::Dirichlet) N--;

    double *a = static_cast<double*>(malloc(sizeof(double)*N));
    double *b = static_cast<double*>(malloc(sizeof(double)*N));
    double *c = static_cast<double*>(malloc(sizeof(double)*N));
    double *d = static_cast<double*>(malloc(sizeof(double)*N));

    unsigned int start=0;
    unsigned int end=N;
    if (leftCondition.boundaryConditionType != BoundaryConditionType::Dirichlet)
    {
        PointNodeODE node(min*h, min);
        b[0] = -2.0 - h*h*B(node) + 2.0*h*leftBetta*(1.0+h*A(node));
        c[0] = leftAlpha;
        d[0] = leftAlpha*h*h*C(node) + h*leftGamma*(1.0+h*A(node));

        start = 1;
    }

    if (rightCondition.boundaryConditionType != BoundaryConditionType::Dirichlet)
    {
        PointNodeODE node(max*h, max);
        b[N-1] = -2.0 - h*h*B(node) - 2.0*h*rightBetta*(1.0+h*A(node));
        c[N-1] = rightAlpha;
        d[N-1] = rightAlpha*h*h*C(node)+h*rightGamma*(1.0+h*A(node));

        end = N-1;
    }

    for (unsigned int i=start; i<end; i++)
    {
        unsigned int mi = static_cast<unsigned int>(min)+i;
        PointNodeODE node(mi*h, static_cast<int>(mi));

        a[i] = +1.0 - 0.5*h*A(node);
        b[i] = -2.0 - h*h*B(node);
        c[i] = +1.0 + 0.5*h*A(node);
        d[i] = h*h*C(node);
    }
    a[0] = c[N-1] = 0.0;

    rv.resize(N);
    tomasAlgorithm(a, b, c, d, rv.data(), N);
}

void SecondOrderLinearODE::solveInitialValueProblem(DoubleVector &rv) const
{
    if (count() != 1) throw ExceptionODE(5);

    const Dimension &dim = dimension();
    int min = dim.min();
    int max = dim.max();
    double h = dim.step();

    unsigned int size = static_cast<unsigned int>(max-min);

    rv.resize(size+1);

    PointNodeODE node0(min*h, min);
    InitialConditionODE ic;
    ic.initialConditionType = InitialConditionType::InitialValue;
    initial(node0, ic);
    double value = ic.value;
    ic.initialConditionType = InitialConditionType::InitialDerivative;
    initial(node0, ic);
    double derivative = ic.value;

    rv[0] = value;
    rv[1] = value + h*derivative + 0.5*h*h*( A(node0)*derivative + B(node0)*value + C(node0) );

    for (unsigned int i=2; i<=size; i++)
    {
        PointNodeODE node(static_cast<double>((i-1)*h), static_cast<int>(i-1));
        double an = A(node);
        double bn = B(node);
        double cn = C(node);
        double mx = (1.0 - 0.5*h*an);

        rv[i] = (2.0+h*h*bn)*rv[i-1] - (1.0+0.5*h*an)*rv[i-2] + h*h*cn;
        rv[i] /= mx;
    }
}

void SecondOrderLinearODE::solveBoundaryValueProblem(const DoubleVector &left, const DoubleVector &right, std::vector<DoubleVector> &ry) const
{

}

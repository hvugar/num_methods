#include "lode2o.h"

void SecondOrderLinearODE::solveInitialValueProblem(DoubleVector &rv) const
{
    if (count() != 1) throw ExceptionODE(5);

    const Dimension &dim = dimension();
    int min = dim.min();
    int max = dim.max();
    double h = dim.step();

    unsigned int size = static_cast<unsigned int>(max-min);

    rv.resize(size+1);

    double value = initial(InitialCondition::InitialValue);
    double derivative = initial(InitialCondition::FirstDerivative);

    PointNodeODE node0(min*h, min);
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

void SecondOrderLinearODE::solveBoundaryValueProblem(DoubleVector &rv) const
{
    if (count() != 1) throw ExceptionODE(5);

    const Dimension &dim = dimension();
    int min = dim.min();
    int max = dim.max();
    double h = dim.step();
    unsigned int size = static_cast<unsigned int>(max-min);
    rv.resize(size+1);

    PointNodeODE leftNode(static_cast<double>(min*h), min);
    BoundaryConditionODE leftCondition;
    double leftValue = boundary(leftNode, leftCondition);
    BoundaryCondition leftConditionType = leftCondition.boundaryConditionType;
    double leftLambda = leftCondition.lambda;

    PointNodeODE rightNode(static_cast<double>(max*h), max);
    BoundaryConditionODE rightCondition;
    double rightValue = boundary(rightNode, rightCondition);
    BoundaryCondition rightConditionType = rightCondition.boundaryConditionType;
    double rightLambda = rightCondition.lambda;

    unsigned int N = size - 1;
    unsigned int start=0;
    unsigned int end=N;

    if (leftConditionType == BoundaryCondition::Dirichlet)
    {
        start = 1;
    }
    else
    {
        N++;
        start = 0;
    }

    if (rightConditionType == BoundaryCondition::Dirichlet)
    {
        end = N-1;
    }
    else
    {
        N++;
        end = N-2;
    }

    double *a = static_cast<double*>(malloc(sizeof(double)*N));
    double *b = static_cast<double*>(malloc(sizeof(double)*N));
    double *c = static_cast<double*>(malloc(sizeof(double)*N));
    double *d = static_cast<double*>(malloc(sizeof(double)*N));
    double *x = static_cast<double*>(malloc(sizeof(double)*N));

    a[0] = 0;
    if (leftCondition.boundaryConditionType == BoundaryCondition::Dirichlet)
    {
        rv[0] = leftValue;
        PointNodeODE node((min+1)*h, min+1);

        b[0] = -2.0 - h*h*B(node);
        c[0] = +1.0 - 0.5*h*A(node);
        d[0] = h*h*C(node) - (1.0 + 0.5*h*A(node))*rv[0];
    }
    else
    {
        PointNodeODE node(min*h, min);
        b[0] = -2.0 - h*h*B(node) - 2.0*h*leftLambda*(1.0+0.5*h*A(node));
        c[0] = +2.0;
        d[0] = h*h*C(node) + 2.0*h*leftValue*(1.0+0.5*h*A(node));
    }

    c[N-1] = 0;
    if (rightCondition.boundaryConditionType == BoundaryCondition::Dirichlet)
    {
        rv[size] = rightValue;
        PointNodeODE node((max-1)*h, static_cast<int>(max-1));

        a[N-1] = +1.0 + 0.5*h*A(node);
        b[N-1] = -2.0 - h*h*B(node);
        d[N-1] = h*h*C(node) - (1.0 - 0.5*h*A(node))*rv[size];
    }
    else
    {
        PointNodeODE node(max*h, max);
        a[N-1] = +2.0;
        b[N-1] = -2.0 - h*h*B(node) + 2.0*h*rightLambda*(1.0-0.5*h*A(node));
        d[N-1] = h*h*C(node) - 2.0*h*rightValue*(1.0-0.5*h*A(node));
    }

    unsigned int j=1;
    for (unsigned int i=start+1; i<=end; i++)
    {
        unsigned int mi = static_cast<unsigned int>(min)+i;
        PointNodeODE node(mi*h, static_cast<int>(mi));

        a[j] = +1.0 + 0.5*h*A(node);
        b[j] = -2.0 - h*h*B(node);
        c[j] = +1.0 - 0.5*h*A(node);
        d[j] = h*h*C(node);
        j++;
    }

    printf("%d\n", N);
    IPrinter::print(a, N);
    IPrinter::print(b, N);
    IPrinter::print(c, N);
    IPrinter::print(d, N);

    tomasAlgorithm(a, b, c, d, x, N);

    j=1;
    for (unsigned int i=0; i<N; i++) rv[start+i] = x[i];

    IPrinter::print(x, N);

    free(x);
    free(d);
    free(c);
    free(b);
    free(a);
}


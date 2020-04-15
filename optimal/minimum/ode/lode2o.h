#ifndef SECOND_ORDER_LINEAR_ODE_H
#define SECOND_ORDER_LINEAR_ODE_H

#include "diffequ.h"

/**
 * @brief Линейное дифференциальное уравнение второго порядка с переменными коэффициентами
 * The Linear ODE 2nd order in canonical (normal) form y"(x) = A(x)y'(x) + B(x)y(x) + C(x);
 */
class MINIMUMSHARED_EXPORT ISecondOrderLinearODEIBVP : virtual public ISecondOrderLinearODE, virtual public InitialValueProblemODE,
        virtual public BoundaryValueProblemODE
{
public:
    PUBLIC_CONSTRUCTORS_VIRTUAL_DESTRUCTOR(ISecondOrderLinearODEIBVP);

    void solveInitialValueProblem(DoubleVector &rv) const;
    void solveBoundaryValueProblem(DoubleVector &rv) const;
    void solveInitialValueProblem(std::vector<DoubleVector> &rv) const;
    void solveBoundaryValueProblem(std::vector<DoubleVector> &rv) const;
    void solveInitialValueProblem(ODESolverMethod method = ODESolverMethod::EULER) const;

    void solveInitialValueProblem2(ODESolverMethod method = ODESolverMethod::EULER) const;
};

class MINIMUMSHARED_EXPORT ISecondOrderLinearODEFBVP : virtual public ISecondOrderLinearODE, virtual public FinalValueProblemODE,
        virtual public BoundaryValueProblemODE
{
public:
    PUBLIC_CONSTRUCTORS_VIRTUAL_DESTRUCTOR(ISecondOrderLinearODEFBVP);

    void solveInitialValueProblem2(ODESolverMethod method = ODESolverMethod::EULER) const;
};

/**
 * @brief The LinearBoundaryValueProblemODE class
 * r(x)y"(x)+p(x)y'(x)+q(x)y(x)=f(x);
 * y(0) = boundary(Left);
 * y(N) = boundart(Right);
 * @see BoundaryValueProblemODE
 * @see BoundaryValueProblem
 */
class MINIMUMSHARED_EXPORT LinearBoundaryValueProblemODE : protected BoundaryValueProblemODE
{
protected:
    virtual double r(unsigned int n) const = 0;
    virtual double p(unsigned int n) const = 0;
    virtual double q(unsigned int n) const = 0;
    virtual double f(unsigned int n) const = 0;

public:
    void calculateX(DoubleVector &x, double h, unsigned int N) const;

    void calculateN2L2RD(DoubleVector &x, double h, unsigned int N) const;
    void calculateN2R2LD(DoubleVector &x, double h, unsigned int N) const;

    void calculateN4L2RD(DoubleVector &x, double h, unsigned int N) const;
    void calculateN4R2LD(DoubleVector &x, double h, unsigned int N) const;
    void calculateN4CNTR(DoubleVector &x, double h, unsigned int N) const;

    void calculateN6L2RD(DoubleVector &x, double h, unsigned int N) const;
    void calculateN6R2LD(DoubleVector &x, double h, unsigned int N) const;
};

#endif // SECOND_ORDER_LINEAR_ODE_H

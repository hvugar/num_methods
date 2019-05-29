#ifndef BOUNDARY_VALUE_PROBLEM_H
#define BOUNDARY_VALUE_PROBLEM_H

#include "grid.h"
#include "../matrix2d.h"

enum class BoundaryCondition
{
    Dirichlet = 1,
    Neumann = 2,
    Robin = 3
};

/**
 * @brief The BoundaryConditionODE class
 */
class MINIMUMSHARED_EXPORT BoundaryConditionODE
{
public:
    BoundaryCondition boundaryConditionType;
    DoubleMatrix a;
    DoubleMatrix b;
    double lambda;
};

/**
 * @brief The BoundaryConditionPDE class
 */
class MINIMUMSHARED_EXPORT BoundaryConditionPDE
{
public:
    BoundaryCondition type;
    double coefficientValue;
    double coefficientDerivative;
};

/**
 * @brief The BoundaryValueProblem class
 */
class MINIMUMSHARED_EXPORT BoundaryValueProblem {};

/**
 * @brief The BoundaryValueProblemODE class
 */
class MINIMUMSHARED_EXPORT BoundaryValueProblemODE : protected BoundaryValueProblem
{
public:
    virtual ~BoundaryValueProblemODE();

protected:
    virtual auto boundary(const PointNodeODE &node, BoundaryConditionODE &condition, unsigned int row = 1) const -> double = 0;
};

/**
 * @brief The BoundaryValueProblemPDE class
 */
class MINIMUMSHARED_EXPORT BoundaryValueProblemPDE : protected BoundaryValueProblem
{
public:
    virtual ~BoundaryValueProblemPDE();

protected:
    virtual auto boundary(const SpaceNodePDE &sn, const TimeNodePDE &tn) const -> double = 0;
};

#endif // BOUNDARYVALUEPROBLEM_H

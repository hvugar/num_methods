#ifndef BOUNDARY_VALUE_PROBLEM_H
#define BOUNDARY_VALUE_PROBLEM_H

#include "grid.h"

enum class BoundaryCondition
{
    Dirichlet = 1,
    Neumann = 2,
    Robin = 3
};

/**
 * @brief The BoundaryConditionPDE class
 */
class MINIMUMSHARED_EXPORT BoundaryConditionPDE
{
public:
    BoundaryConditionPDE(BoundaryCondition condition = BoundaryCondition::Dirichlet, double alpha = 1.0, double beta = 0.0, double gamma = 1.0);

    BoundaryCondition boundaryCondition() const;
    double alpha() const;
    double beta() const;
    double gamma() const;

private:
    BoundaryCondition _boundaryCondition;
    double _alpha;
    double _beta;
    double _gamma;
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
    PUBLIC_CONSTRUCTORS_VIRTUAL_DESTRUCTOR(BoundaryValueProblemODE);

protected:
    virtual auto boundary(const PointNodeODE &node, BoundaryConditionPDE &condition, unsigned int row = 1) const -> double = 0;
};

/**
 * @brief The BoundaryValueProblemPDE class
 */
class MINIMUMSHARED_EXPORT BoundaryValueProblemPDE : protected BoundaryValueProblem
{
    PUBLIC_CONSTRUCTORS_VIRTUAL_DESTRUCTOR(BoundaryValueProblemPDE);

protected:
    virtual auto boundary(const SpaceNodePDE &sn, const TimeNodePDE &tn, BoundaryConditionPDE &condition) const -> double = 0;
};

#endif // BOUNDARYVALUEPROBLEM_H

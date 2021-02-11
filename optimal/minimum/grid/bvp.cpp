#include "bvp.h"

BoundaryValueProblemODE::BoundaryValueProblemODE() {}

BoundaryValueProblemODE::BoundaryValueProblemODE(const BoundaryValueProblemODE &) {}

BoundaryValueProblemODE& BoundaryValueProblemODE::operator=(const BoundaryValueProblemODE &other)
{
    if (this == &other) { return *this; }
    return *this;
}

BoundaryValueProblemODE::~BoundaryValueProblemODE() {}

//**********************************************************************************************//

BoundaryValueProblemPDE::BoundaryValueProblemPDE() {}

BoundaryValueProblemPDE::BoundaryValueProblemPDE(const BoundaryValueProblemPDE &) {}

BoundaryValueProblemPDE& BoundaryValueProblemPDE::operator=(const BoundaryValueProblemPDE &other)
{
    if (this == &other) { return *this; }
    return *this;
}

BoundaryValueProblemPDE::~BoundaryValueProblemPDE() {}

//**********************************************************************************************//

BoundaryConditionPDE::BoundaryConditionPDE(BoundaryCondition condition)
{
    this->_boundaryCondition = condition;
    this->_alpha = 0.0;
    this->_beta = 0.0;
}

BoundaryConditionPDE BoundaryConditionPDE::Dirichlet()
{
    return BoundaryConditionPDE(BoundaryCondition::Dirichlet);
}

BoundaryConditionPDE BoundaryConditionPDE::Neumann()
{
    return BoundaryConditionPDE(BoundaryCondition::Neumann);
}

BoundaryConditionPDE BoundaryConditionPDE::Robin(double alpha, double beta)
{
    BoundaryConditionPDE condition(BoundaryCondition::Robin);
    condition._alpha = alpha;
    condition._beta = beta;
    return condition;
}

BoundaryCondition BoundaryConditionPDE::boundaryCondition() const { return _boundaryCondition; }

double BoundaryConditionPDE::alpha() const { return _alpha; }

double BoundaryConditionPDE::beta() const { return _beta; }

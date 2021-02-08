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

BoundaryConditionPDE::BoundaryConditionPDE(BoundaryCondition condition, double alpha, double beta)
{
    if (condition == BoundaryCondition::Robin)  { if (fabs(alpha) <= DBL_EPSILON || fabs(beta) <= DBL_EPSILON) { throw std::exception(); } }

    this->_boundaryCondition = condition;
    this->_alpha = alpha;
    this->_beta = beta;
}

BoundaryConditionPDE BoundaryConditionPDE::Dirichlet() { return BoundaryConditionPDE(BoundaryCondition::Dirichlet); }

BoundaryConditionPDE BoundaryConditionPDE::Neumann() { return BoundaryConditionPDE(BoundaryCondition::Neumann); }

BoundaryConditionPDE BoundaryConditionPDE::Robin(double alpha, double beta) { return BoundaryConditionPDE(BoundaryCondition::Robin, alpha, beta); }

BoundaryCondition BoundaryConditionPDE::boundaryCondition() const { return _boundaryCondition; }

double BoundaryConditionPDE::alpha() const { return _alpha; }

double BoundaryConditionPDE::beta() const { return _beta; }

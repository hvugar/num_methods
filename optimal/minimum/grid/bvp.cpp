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

BoundaryConditionPDE::BoundaryConditionPDE(BoundaryCondition condition, double alpha, double beta, double gamma)
{
    if (condition == BoundaryCondition::Dirichlet) { if (fabs(alpha) <= DBL_EPSILON || fabs(beta) >= DBL_EPSILON) { throw std::exception(); } } else
    if (condition == BoundaryCondition::Neumann)   { if (fabs(alpha) >= DBL_EPSILON || fabs(beta) <= DBL_EPSILON) { throw std::exception(); } } else
    if (condition == BoundaryCondition::Robin)     { if (fabs(alpha) <= DBL_EPSILON || fabs(beta) <= DBL_EPSILON) { throw std::exception(); } }

    this->_boundaryCondition = condition;
    this->_alpha = alpha;
    this->_beta = beta;
    this->_gamma = gamma;
}

BoundaryConditionPDE BoundaryConditionPDE::Dirichlet(double alpha, double beta, double gamma)
{
    return BoundaryConditionPDE(BoundaryCondition::Dirichlet, alpha, beta, gamma);
}

BoundaryConditionPDE BoundaryConditionPDE::Neumann(double alpha, double beta, double gamma)
{
    return BoundaryConditionPDE(BoundaryCondition::Neumann, alpha, beta, gamma);
}

BoundaryConditionPDE BoundaryConditionPDE::Robin(double alpha, double beta, double gamma)
{
    return BoundaryConditionPDE(BoundaryCondition::Robin, alpha, beta, gamma);
}

BoundaryCondition BoundaryConditionPDE::boundaryCondition() const { return _boundaryCondition; }

double BoundaryConditionPDE::alpha() const { return _alpha; }

double BoundaryConditionPDE::beta() const { return _beta; }

double BoundaryConditionPDE::gamma() const { return _gamma; }


#include "bvp.h"

BoundaryValueProblemODE::~BoundaryValueProblemODE() {}

BoundaryValueProblemPDE::~BoundaryValueProblemPDE() {}

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

BoundaryCondition BoundaryConditionPDE::boundaryCondition() const { return _boundaryCondition; }

double BoundaryConditionPDE::alpha() const { return _alpha; }

double BoundaryConditionPDE::beta() const { return _beta; }

double BoundaryConditionPDE::gamma() const { return _gamma; }


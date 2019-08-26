#include "problem1p_solver.h"

using namespace p1p;

double HeatEquationIBVP::initial(const SpaceNodePDE &sn, InitialCondition condition) const
{}

double HeatEquationIBVP::boundary(const SpaceNodePDE &sn, const TimeNodePDE &tn, BoundaryConditionPDE &condition) const
{}

double HeatEquationIBVP::f(const SpaceNodePDE &sn, const TimeNodePDE &tn) const
{}

void HeatEquationIBVP::layerInfo(const DoubleVector &u, const TimeNodePDE &tn) const
{}


double ConjugateHeatEquationIBVP::final(const SpaceNodePDE &sn, FinalCondition condition) const
{}

double ConjugateHeatEquationIBVP::boundary(const SpaceNodePDE &sn, const TimeNodePDE &tn, BoundaryConditionPDE &condition) const
{}

double ConjugateHeatEquationIBVP::f(const SpaceNodePDE &sn, const TimeNodePDE &tn) const
{}

void ConjugateHeatEquationIBVP::layerInfo(const DoubleVector &u, const TimeNodePDE &tn) const
{}

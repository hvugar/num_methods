#include "diffequ.h"

SystemDifferentialEquation::SystemDifferentialEquation(const ODEGrid &grid) : mgrid(grid) {}

const ODEGrid& SystemDifferentialEquation::grid() const { return mgrid; }

SystemLinearODE1stOrder::SystemLinearODE1stOrder(const ODEGrid &grid) : SystemLinearODE(grid) {}

#ifndef PROBLEM1H_SOLVER1_H
#define PROBLEM1H_SOLVER1_H

#include "problem1h_solver_base.h"

class PROBLEM1HSHARED_EXPORT Problem1HDirichlet1 : public Problem1HDirichletBase
{
public:
    Problem1HDirichlet1();
    virtual ~Problem1HDirichlet1();

    virtual auto gradient(const DoubleVector &, DoubleVector &) const -> void;

    /** Integral part of functional */
    virtual auto integral(const std::vector<DoubleVector> &u) const -> double;

    /** Penalty part of functional */
    virtual auto penalty(const spif_vector1H &info, const OptimizeParameter1H &o_prm) const -> double;

public:
    auto solveForwardIBVP(std::vector<DoubleVector> &u, spif_vector1H &u_info, bool use, double lambda=0.25) const -> void;
    auto solveBackwardIBVP(const std::vector<DoubleVector> &u, spif_vector1H &p_info, bool use, const spif_vector1H &u_info, double lambda=0.25) const -> void;
};

#endif // PROBLEM1H_SOLVER1_H

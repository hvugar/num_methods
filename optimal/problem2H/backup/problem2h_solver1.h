#ifndef PROBLEM2H_SOLVER1_H
#define PROBLEM2H_SOLVER1_H

#include "problem2h_solver_base.h"

class PROBLEM2HSHARED_EXPORT Problem2HDirichlet1 : public Problem2HDirichletBase
{
public:
    Problem2HDirichlet1();
    virtual ~Problem2HDirichlet1();

    virtual auto gradient(const DoubleVector &, DoubleVector &) const -> void;

    /** Integral part of functional */
    virtual auto integral(const std::vector<DoubleMatrix> &u) const -> double;

    /** Penalty part of functional */
    virtual auto penalty(const spif_vectorH &info, const OptimizeParameterH &o_prm) const -> double;
    virtual auto gpi(unsigned int i, unsigned int ln, const spif_vectorH &u_info, const OptimizeParameterH &o_prm) const -> double;
    virtual auto g0i(unsigned int i, unsigned int ln, const spif_vectorH &u_info, const OptimizeParameterH &o_prm) const -> double;
    virtual auto v(unsigned int i, OptimizeParameterH o_prm, EquationParameterH e_prm, const spif_vectorH &u_info, unsigned int s) const -> double;

    auto solveForwardIBVP(std::vector<DoubleMatrix> &u, spif_vectorH &u_info, bool use, const DoubleVector &pv, double lambda=0.25) const -> void;
    auto solveBackwardIBVP(const std::vector<DoubleMatrix> &u, spif_vectorH &p_info, bool use, const spif_vectorH &u_info, const DoubleVector &pv, double lambda=0.25) const -> void;
};

#endif // PROBLEM2H_SOLVER1_H

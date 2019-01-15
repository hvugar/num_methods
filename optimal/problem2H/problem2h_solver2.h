#ifndef PROBLEM2H_SOLVER2_H
#define PROBLEM2H_SOLVER2_H

#include "problem2h_solver_base.h"

class PROBLEM2HSHARED_EXPORT Problem2HDirichlet2 : public Problem2hDirichletBase
{
public:
    Problem2HDirichlet2();
    virtual ~Problem2HDirichlet2();

    virtual auto gradient(const DoubleVector &, DoubleVector &) const -> void;

    /** Integral part of functional */
    virtual auto integral(const std::vector<DoubleMatrix> &u) const -> double;

    /** Penalty part of functional */
    virtual auto penalty(const spif_vectorH &info, const OptimizeParameterH &o_prm) const -> double;

public:
    auto solveForwardIBVP(std::vector<DoubleMatrix> &u, spif_vectorH &u_info, bool use, const DoubleVector &pv, double lambda=0.25) const -> void;
    auto solveBackwardIBVP(const std::vector<DoubleMatrix> &u, spif_vectorH &p_info, bool use, const spif_vectorH &u_info, const DoubleVector &pv, double lambda=0.25) const -> void;

    virtual auto initConjtWeightMatrix(const std::vector<DoubleMatrix> &u) const -> void;

    virtual auto b_initial1(const SpaceNodePDE &sn) const -> double;
    virtual auto b_initial2(const SpaceNodePDE &sn) const -> double;

private:
    DoubleMatrix b_init1LayerMatrix;
    DoubleMatrix b_init2LayerMatrix;
};

#endif // PROBLEM2H_SOLVER2_H

#ifndef PROBLEM2H_IBVP_H
#define PROBLEM2H_IBVP_H

#include "problem2h_common.h"
#include <grid/hibvp.h>

class PROBLEM2HSHARED_EXPORT Problem2HNDirichletForward1 : public CC1IHyperbolicIBVP
{
public:
    const EquationParameterH *e_prm = nullptr;
    const OptimizeParameterH *o_prm = nullptr;
    const DoubleMatrix *u10 = nullptr;

protected:
    virtual double initial1(const SpaceNodePDE &sn) const;
    virtual double initial2(const SpaceNodePDE &sn) const;
    virtual double boundary(const SpaceNodePDE &sn, const TimeNodePDE &tn) const;
    virtual double f(const SpaceNodePDE &sn, const TimeNodePDE &) const;

public:
    void setU10(const DoubleMatrix &u10);
    void setEquationParameters(const EquationParameterH &e_prm, const OptimizeParameterH& o_prm);

private:
    GridSpace2D grid;
    bool _initialCalculation = false;
    std::vector<double> u_xi;
    std::vector<DeltaGrid> msrmGridList;
};

#endif // PROBLEM2H_IBVP_H

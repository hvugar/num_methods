#ifndef ISYSTEM_LINEAR_ODE_NONLOCALCONTIONSM_H
#define ISYSTEM_LINEAR_ODE_NONLOCALCONTIONSM_H

#include <vector2d.h>
#include <matrix2d.h>
#include <grid/grid.h>
#include <ode/cauchyp.h>
#include <vector>
#include "islodenlcsv.h"

class MINIMUMSHARED_EXPORT ISystemLinearODENonLocalContionsM : public SystemLinearODE
{
public:
    virtual void calculateForward(DoubleMatrix &x);
    virtual void calculateBackward(DoubleMatrix &x);

    virtual void calculateBackwardCP(DoubleMatrix &x, std::vector<std::vector<DoubleVector>> &m);

    void setLeftSeparatedCondition(const ISystemLinearODENonLocalContionsV::Condition &lscs);
    void setRightSeparatedCondition(const ISystemLinearODENonLocalContionsV::Condition &lscs);
    void addNonSeparatedCondition(const ISystemLinearODENonLocalContionsV::Condition &nsc);
    const std::vector<ISystemLinearODENonLocalContionsV::Condition>& nonSeparatedConditions() const;
    void setBetta(const DoubleMatrix &betta);
    void setSystemOrder(unsigned int n);
    unsigned int systemOrder() const;
public:
    virtual double A(double t UNUSED_PARAM, unsigned int k, unsigned int row, unsigned int col) const = 0;
    virtual double B(double t UNUSED_PARAM, unsigned int k, unsigned int row, unsigned int col) const = 0;

    std::vector<ISystemLinearODENonLocalContionsV::Condition> nscs;
    ISystemLinearODENonLocalContionsV::Condition lscs;
    ISystemLinearODENonLocalContionsV::Condition rscs;
    DoubleMatrix betta;
};

#endif // ISYSTEM_LINEAR_ODE_NONLOCALCONTIONSM_H

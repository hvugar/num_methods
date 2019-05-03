#ifndef ISYSTEM_LINEAR_ODE_NONLOCALCONTIONSM_H
#define ISYSTEM_LINEAR_ODE_NONLOCALCONTIONSM_H

#include <vector2d.h>
#include <matrix2d.h>
#include <grid/grid.h>
//#include <ode/cauchyp.h>
#include <ode/nlode1o.h>
#include <vector>
#include "islodenlcsv.h"

//Depreciated

//class MINIMUMSHARED_EXPORT ISystemLinearODENonLocalContionsM : public SystemDifferentialEquation, public ISystemLinearODENonLocalContions
//{
//public:
//    virtual void calculateForward(DoubleMatrix &x);
//    virtual void calculateBackward(DoubleMatrix &x);

//    virtual void calculateBackwardCP(DoubleMatrix &x, std::vector<std::vector<DoubleVector>> &m);

//    void setLeftSeparatedCondition(const Condition &lscs);
//    void setRightSeparatedCondition(const Condition &lscs);
//    void addNonSeparatedCondition(const Condition &nsc);
//    const std::vector<Condition>& nonSeparatedConditions() const;
//    void setBetta(const DoubleMatrix &betta);
//    void setSystemOrder(unsigned int n);
//    unsigned int systemOrder() const;
//public:
//    virtual double A(double t UNUSED_PARAM, unsigned int k, unsigned int row, unsigned int col) const = 0;
//    virtual double B(double t UNUSED_PARAM, unsigned int k, unsigned int row, unsigned int col) const = 0;

//    std::vector<Condition> nscs;
//    Condition lscs;
//    Condition rscs;
//    DoubleMatrix betta;
//};

#endif // ISYSTEM_LINEAR_ODE_NONLOCALCONTIONSM_H

#ifndef ISYSTEM_LINEAR_ODE_NONLOCALCONTIONSV_H
#define ISYSTEM_LINEAR_ODE_NONLOCALCONTIONSV_H

#include <vector2d.h>
#include <matrix2d.h>
#include <vector>
//#include <ode/cauchyp.h>
#include <ode/nlode1o.h>
#include "islodenlcs.h"

using namespace std;

//Depreciated

/**
 * @brief The CauchyProblemNonLocalContions class
 * Численное решение систем дифференциальных уравнений с нелокальными условиями.
 * Numerical solution of systems of linear ordinary differential equations with non-local conditions.
 */

//class MINIMUMSHARED_EXPORT ISystemLinearODENonLocalContionsV : public SystemLinearODE1stOrder, public ISystemLinearODENonLocalContions
//{
//public:
//    void setLeftSeparatedCondition(const Condition &lscs);
//    void setRightSeparatedCondition(const Condition &lscs);
//    void addNonSeparatedCondition(const Condition &nsc);
//    const std::vector<Condition>& nonSeparatedConditions() const;
//    void setBetta(const DoubleVector &betta);
//    void setSystemOrder(unsigned int n);
//    unsigned int systemOrder() const;

//    void calculateForward(DoubleVector &x);
//    void calculateBackwardCP(const DoubleVector &x, std::vector<DoubleVector> &m);

//    void calculateBackward(DoubleVector &x);
//    void calculateForwardCP(const DoubleVector &x, std::vector<DoubleVector> &m);

//private:
//    /** @brief nlcs Неразделенные многоточечные условия заданные в интервале. */
//    std::vector<Condition> nscs;
//    /** @brief nlcs Разделенные условия заданные на левом конце интервала. */
//    Condition lscs;
//    /** @brief nlcs Разделенные условия заданные на правом конце интервала. */
//    Condition rscs;
//    /** @brief betta */
//    DoubleVector betta;

//private:
//    void calculateIntervalF(unsigned int s, unsigned int r);
//    void calculateIntervalB(unsigned int s, unsigned int r);
////    void calculateDiffEquation(const Condition &sc, const Condition &ec, double h, const DoubleVector &x, DoubleVector &rx,
////                               unsigned int minN, unsigned int maxN);

//public:
//    virtual double A(double t UNUSED_PARAM, unsigned int k, unsigned int row, unsigned int col) const = 0;
//    virtual double B(double t UNUSED_PARAM, unsigned int k, unsigned int row) const = 0;
//};

#endif // ISYSTEM_LINEAR_ODE_NONLOCALCONTIONSV_H

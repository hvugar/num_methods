#ifndef SYSTEM_LINEARODE_NONLOCALCONTIONS_H
#define SYSTEM_LINEARODE_NONLOCALCONTIONS_H

#include <vector2d.h>
#include <matrix2d.h>
#include <vector>
#include <grid/cauchyp.h>

using namespace std;

/**
 * @brief The CauchyProblemNonLocalContions class
 * Численное решение систем дифференциальных уравнений с нелокальными условиями.
 * Numerical solution of systems of linear ordinary differential equations with non-local conditions.
 */

class MINIMUMSHARED_EXPORT ISystemLinearODENonLocalContions : public SystemLinearODE1stOrder
{
public:
    enum ConditionType
    {
        SeparatedLeft = 0,
        SeparatedRight = 1,
        NonSeparated = 2
    };

    struct Condition
    {
        ConditionType type;
        double time;
        unsigned int nmbr;
        DoubleMatrix alpha;
    };

    ISystemLinearODENonLocalContions(const ODEGrid &grid);

    void setLeftSeparatedCondition(const Condition &lscs);
    void setRightSeparatedCondition(const Condition &lscs);
    void addNonSeparatedCondition(const Condition &nsc);
    const std::vector<Condition>& nonSeparatedConditions() const;
    void setBetta(const DoubleVector &betta);
    void setSystemOrder(unsigned int n);
    unsigned int systemOrder() const;

    void calculateBackwardCP(const DoubleVector &x, DoubleMatrix &m);
    void calculateForwardCP(const DoubleVector &x, DoubleMatrix &m);

private:
    /**
     * @brief nlcs Неразделенные многоточечные условия заданные в интервале.
     */
    std::vector<Condition> nscs;
    /**
     * @brief nlcs Разделенные условия заданные на левом конце интервала.
     */
    Condition lscs;
    /**
     * @brief nlcs Разделенные условия заданные на правом конце интервала.
     */
    Condition rscs;
    /**
     * @brief betta
     */
    DoubleVector betta;
    /**
     * @brief n0 Число неразделенных многоточечных условий заданных в интервале.
     */
    unsigned int n0;
    /**
     * @brief n1 Число разделенных условий заданных на левом конце интервала.
     */
    unsigned int n1;
    /**
     * @brief n1 Число разделенных условий заданных на правом конце интервала.
     */
    unsigned int n2;

    unsigned int n;
    unsigned int L;

public:
    void calculateForward(DoubleVector &x);
    void calculateBackward(DoubleVector &x);

private:
    void calculateIntervalF(unsigned int s, unsigned int r);
    void calculateIntervalB(unsigned int s, unsigned int r);

public:
    virtual double A(double t UNUSED_PARAM, unsigned int k, unsigned int row, unsigned int col) const = 0;
    virtual double B(double t UNUSED_PARAM, unsigned int k, unsigned int row) const = 0;
};

#endif // CAUCHYPROBLEMNONLOCALCONTIONS_H

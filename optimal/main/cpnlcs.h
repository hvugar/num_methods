#ifndef CAUCHYPROBLEMNONLOCALCONTIONS_H
#define CAUCHYPROBLEMNONLOCALCONTIONS_H

#include <vector2d.h>
#include <matrix2d.h>
#include <vector>
#include <grid/cauchyp.h>

#define SAMPLE_3

using namespace std;

/**
 * @brief The CauchyProblemNonLocalContions class
 * Численное решение систем дифференциальных уравнений с нелокальными условиями.
 * Numerical solution of systems of linear ordinary differential equations with non-local conditions.
 */

class CauchyProblemNonLocalContions : public SystemLinearODE1stOrder
{
public:
    static void Main(int agrc, char *argv[]);

    enum ConditionType
    {
        SeparatedLeft = 0,
        SeparatedRight = 1,
        NonSeparated = 2
    } ;

    struct Condition
    {
        ConditionType type;
        double time;
        unsigned int nmbr;
        DoubleMatrix alpha;
    };

    CauchyProblemNonLocalContions(const ODEGrid &grid);

    virtual void initialize();

    /**
     * @brief n0 Число неразделенных многоточечных условий заданных в интервале.
     */
    unsigned int n0;
    /**
     * @brief nlcs Неразделенные многоточечные условия заданные в интервале.
     */
    std::vector<Condition> nscs;
    /**
     * @brief n1 Число разделенных условий заданных на левом конце интервала.
     */
    unsigned int n1;
    /**
     * @brief nlcs Разделенные условия заданные на левом конце интервала.
     */
    Condition lscs;
    /**
     * @brief n1 Число разделенных условий заданных на правом конце интервала.
     */
    unsigned int n2;
    /**
     * @brief nlcs Разделенные условия заданные на правом конце интервала.
     */
    Condition rscs;
    /**
     * @brief betta
     */
    DoubleVector betta;

    unsigned int n;
    unsigned int L;

    double x(unsigned int k, int i) const;

    void calculateForward(DoubleVector &x);
    void calculateBackward(DoubleVector &x);

    void calculateIntervalF(unsigned int s, unsigned int r);
    void calculateIntervalB(unsigned int s, unsigned int r);

public:
    virtual double A(double t UNUSED_PARAM, unsigned int k, unsigned int row, unsigned int col) const;
    virtual double B(double t UNUSED_PARAM, unsigned int k, unsigned int row) const;
};

class CauchyProblemM1stOrderA : public CauchyProblemM1stOrder
{
public:
    CauchyProblemM1stOrderA(CauchyProblemNonLocalContions &parent, const ODEGrid& grid) : CauchyProblemM1stOrder(grid), p(parent) {}

protected:
    virtual double f(double t, const DoubleVector &x, unsigned int k, unsigned int i) const;
    double S0(double t, const DoubleVector &x, unsigned int k) const;
private:
    CauchyProblemNonLocalContions &p;
};

class CauchyProblemM1stOrderB : public CauchyProblemM1stOrder
{
public:
    CauchyProblemM1stOrderB(CauchyProblemNonLocalContions &parent, const ODEGrid& grid) : CauchyProblemM1stOrder(grid), p(parent) {}
protected:
    virtual double f(double t, const DoubleVector &x, unsigned int k, unsigned int i) const;
private:
    CauchyProblemNonLocalContions &p;
};

#endif // CAUCHYPROBLEMNONLOCALCONTIONS_H

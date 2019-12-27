#ifndef FIRTS_ORDER_LINEAR_ODE_H
#define FIRTS_ORDER_LINEAR_ODE_H

#include "diffequ.h"
#include "../grid/ibvp.h"

struct MINIMUMSHARED_EXPORT NonLocalCondition
{
    NonLocalCondition();
    NonLocalCondition(unsigned int i, const PointNodeODE &node, const DoubleMatrix &m);
    virtual ~NonLocalCondition();

    unsigned int i;
    PointNodeODE n;
    DoubleMatrix m;
};

struct Condition
{
    double time;
    DoubleMatrix mtrx;
    unsigned int nmbr;
    unsigned int index;
};

/**
 * @brief Линейное дифференциальное уравнение первого порядка с переменными коэффициентами
 * The Linear ODE 1st order in canonical (normal) form y'(x) = A(x)y(x) + B(x);
 */
class MINIMUMSHARED_EXPORT FirstOrderLinearODE :
        virtual public LinearODE,
        virtual public InitialValueProblemODE,
        virtual public BoundaryValueProblemODE
{
public:
    enum class AccuracyStep
    {
        Step_2 = 2,
        Step_4 = 4,
        Step_6 = 6
    };

    //FirstOrderLinearODE();
    virtual ~FirstOrderLinearODE();

    /**
     * @brief solve transfer of conditions
     * @param C
     * @param d
     * @param x
     * @param k
     * @param M
     * @param direction
     */
    void transferOfCondition(const std::vector<NonLocalCondition> &C, const DoubleVector &d, std::vector<DoubleVector> &x, unsigned int k) const;
    void transferOfCondition1(const std::vector<NonLocalCondition> &C, const DoubleVector &d, std::vector<DoubleVector> &x, unsigned int k) const;

    /**
     * @brief solveInitialValueProblem
     * @param rv
     */
    void solveInitialValueProblem(DoubleVector &rv) const;

    /**
     * @brief solveInitialValueProblem
     * @param rv
     */
    void solveInitialValueProblem(std::vector<DoubleVector> &rv, ODESolverMethod method = ODESolverMethod::EULER) const;

protected:

    /**
     * @brief A  A nxn dimensional matrix-function
     * @param node
     * @param row <= n
     * @param col <= n
     * @return
     */
    virtual auto A(const PointNodeODE &node, unsigned int row = 1, unsigned int col = 1) const -> double = 0;

    /**
     * @brief B n dimensional vector-function
     * @param node
     * @param row
     * @return
     */
    virtual auto B(const PointNodeODE &node, unsigned int row = 1) const -> double = 0;

protected:
    /**
     * @brief initial
     * @param condition
     * @param row
     * @return
     */
    virtual auto initial(InitialCondition condition, unsigned int row = 1) const -> double = 0;

    /**
     * @brief boundary
     * @param node
     * @param condition
     * @param row
     * @return
     */
    virtual auto boundary(const PointNodeODE &node, BoundaryConditionODE &condition, unsigned int row = 1) const -> double = 0;

protected:
    virtual auto count() const -> unsigned int = 0;

private:

    auto discritize(const std::vector<NonLocalCondition> &co, std::vector<NonLocalCondition> &cn, unsigned int k = 4) const -> void;

    auto solveInitialValueProblemRK2(std::vector<DoubleVector> &rv) const -> void;
    auto solveInitialValueProblemRK4(std::vector<DoubleVector> &rv) const -> void;
    auto solveInitialValueProblemRK6(std::vector<DoubleVector> &rv) const -> void;
    auto solveInitialValueProblemEuler(std::vector<DoubleVector> &rv) const -> void;
    auto solveInitialValueProblemEulerMod(std::vector<DoubleVector> &rv) const -> void;
};


#endif // FIRTS_ORDER_LINEAR_ODE_H

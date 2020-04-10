#ifndef FIRTS_ORDER_LINEAR_ODE_H
#define FIRTS_ORDER_LINEAR_ODE_H

#include "diffequ.h"
#include "../grid/ibvp.h"

struct MINIMUMSHARED_EXPORT NonLocalCondition
{
    NonLocalCondition();
    NonLocalCondition(unsigned int i, const PointNodeODE &node, const DoubleMatrix &m);
    virtual ~NonLocalCondition();

    /**
     * @brief i index of non-local condition
     */
    unsigned int i;
    /**
     * @brief n point-node on grid
     */
    PointNodeODE n;
    /**
     * @brief m non-local condition matrix
     */
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
class MINIMUMSHARED_EXPORT IFirstOrderLinearODE :
        virtual public LinearODE,
        virtual public InitialValueProblemODE
{
public:
    enum class AccuracyStep
    {
        Step_2 = 2,
        Step_4 = 4,
        Step_6 = 6
    };

    PUBLIC_CONSTRUCTORS_VIRTUAL_DESTRUCTOR(IFirstOrderLinearODE);

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
    void transferOfConditionN(const std::vector<NonLocalCondition> &C, const DoubleVector &d, std::vector<DoubleVector> &x, unsigned int k, unsigned int schema) const;
    void transferOfConditionM(const std::vector<NonLocalCondition> &C, const DoubleVector &d, std::vector<DoubleVector> &x, unsigned int k) const;
    void transferOfConditionS(const std::vector<NonLocalCondition> &C, const DoubleVector &d, std::vector<DoubleVector> &x, unsigned int k, unsigned int schema) const;
    void transferOfConditionP(const std::vector<NonLocalCondition> &C, const DoubleVector &d, std::vector<DoubleVector> &x, unsigned int k, unsigned int schema) const;

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

    /**
     * @brief initial
     * @param condition
     * @param row
     * @return
     */
    virtual auto initial(InitialCondition condition, unsigned int row = 1) const -> double = 0;

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

class MINIMUMSHARED_EXPORT IFirstOrderLinearODEFBVP :
        virtual public LinearODE,
        virtual public FinalValueProblemODE
{
    PUBLIC_CONSTRUCTORS_VIRTUAL_DESTRUCTOR(IFirstOrderLinearODEFBVP);

public:
    void solveFinalValueProblem(std::vector<DoubleVector> &rv, ODESolverMethod method = ODESolverMethod::EULER) const;

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

    /**
     * @brief initial
     * @param condition
     * @param row
     * @return
     */
    virtual auto final(FinalCondition condition, unsigned int row = 1) const -> double = 0;

private:

    auto solveFinalValueProblemEuler(std::vector<DoubleVector> &rv) const -> void;
    auto solveFinalValueProblemEulerMod(std::vector<DoubleVector> &rv) const -> void;
    auto solveFinalValueProblemRK2(std::vector<DoubleVector> &rv) const -> void;
    auto solveFinalValueProblemRK4(std::vector<DoubleVector> &rv) const -> void;
    auto solveFinalValueProblemRK6(std::vector<DoubleVector> &rv) const -> void;
};


#endif // FIRTS_ORDER_LINEAR_ODE_H

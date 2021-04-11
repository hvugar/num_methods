#ifndef FIRTS_ORDER_LINEAR_ODE_H
#define FIRTS_ORDER_LINEAR_ODE_H

#include "diffequ.h"

struct MINIMUMSHARED_EXPORT NonLocalCondition
{
    NonLocalCondition();
    NonLocalCondition(size_t i, const PointNodeODE &node, const DoubleMatrix &m);
    virtual ~NonLocalCondition();

    /**
     * @brief i index of non-local condition
     */
    size_t i;
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
class MINIMUMSHARED_EXPORT IFirstOrderLinearODEIVP : virtual public IFirstOrderLinearODE,
                                                     virtual public InitialValueProblemODE
{
public:
    PUBLIC_CONSTRUCTORS_VIRTUAL_DESTRUCTOR(IFirstOrderLinearODEIVP);

    void solveInitialValueProblem(DoubleVector &rv) const;
    void solveInitialValueProblem(ODESolverMethod method = ODESolverMethod::EULER) const;
    void solveInitialValueProblem(std::vector<DoubleVector> &rv, ODESolverMethod method = ODESolverMethod::EULER) const;

    void start(DoubleVector &x, PointNodeODE &n) const;
    void next(const DoubleVector &x0, const PointNodeODE &n0, DoubleVector &x, PointNodeODE &n, ODESolverMethod method = ODESolverMethod::EULER) const;

    void transferOfCondition(const std::vector<NonLocalCondition> &C, const DoubleVector &d, std::vector<DoubleVector> &x, unsigned int k) const;
    void transferOfConditionN(const std::vector<NonLocalCondition> &C, const DoubleVector &d, std::vector<DoubleVector> &x, unsigned int k, unsigned int schema) const;
    void transferOfConditionM(const std::vector<NonLocalCondition> &C, const DoubleVector &d, std::vector<DoubleVector> &x, unsigned int k) const;
    void transferOfConditionS(const std::vector<NonLocalCondition> &C, const DoubleVector &d, std::vector<DoubleVector> &x, unsigned int k, unsigned int schema) const;
    void transferOfConditionP(const std::vector<NonLocalCondition> &C, const DoubleVector &d, std::vector<DoubleVector> &x, unsigned int k, unsigned int schema) const;
private:

    auto discritize(const std::vector<NonLocalCondition> &co, std::vector<NonLocalCondition> &cn, unsigned int k = 4) const -> void;

    auto solveInitialValueProblemRK2(std::vector<DoubleVector> &rv) const -> void;
    auto solveInitialValueProblemRK4(std::vector<DoubleVector> &rv) const -> void;
    auto solveInitialValueProblemRK6(std::vector<DoubleVector> &rv) const -> void;
    auto solveInitialValueProblemEuler(std::vector<DoubleVector> &rv) const -> void;
    auto solveInitialValueProblemEulerMod(std::vector<DoubleVector> &rv) const -> void;

    auto solveInitialValueProblemRK2() const -> void;
    auto solveInitialValueProblemRK4() const -> void;
    auto solveInitialValueProblemRK6() const -> void;
    auto solveInitialValueProblemEuler() const -> void;
    auto solveInitialValueProblemEulerMod() const -> void;

};

/**
 * @brief Линейное дифференциальное уравнение первого порядка с переменными коэффициентами
 * The Linear ODE 1st order in canonical (normal) form y'(x) = A(x)y(x) + B(x);
 */
class MINIMUMSHARED_EXPORT IFirstOrderLinearODEFVP : virtual public IFirstOrderLinearODE,
                                                     virtual public FinalValueProblemODE
{
public:
    PUBLIC_CONSTRUCTORS_VIRTUAL_DESTRUCTOR(IFirstOrderLinearODEFVP);

    void solveFinalValueProblem(DoubleVector &rv) const;
    void solveFinalValueProblem(ODESolverMethod method = ODESolverMethod::EULER) const;
    void solveFinalValueProblem(std::vector<DoubleVector> &rv, ODESolverMethod method = ODESolverMethod::EULER) const;

    void start(DoubleVector &x, PointNodeODE &n) const;
    void next(const DoubleVector &x0, const PointNodeODE &n0, DoubleVector &x, PointNodeODE &n, ODESolverMethod method = ODESolverMethod::EULER) const;

private:

    auto solveFinalValueProblemRK2(std::vector<DoubleVector> &rv) const -> void;
    auto solveFinalValueProblemRK4(std::vector<DoubleVector> &rv) const -> void;
    auto solveFinalValueProblemRK6(std::vector<DoubleVector> &rv) const -> void;
    auto solveFinalValueProblemEuler(std::vector<DoubleVector> &rv) const -> void;
    auto solveFinalValueProblemEulerMod(std::vector<DoubleVector> &rv) const -> void;

    auto solveFinalValueProblemRK2() const -> void;
    auto solveFinalValueProblemRK4() const -> void;
    auto solveFinalValueProblemRK6() const -> void;
    auto solveFinalValueProblemEuler() const -> void;
    auto solveFinalValueProblemEulerMod() const -> void;
};


#endif // FIRTS_ORDER_LINEAR_ODE_H

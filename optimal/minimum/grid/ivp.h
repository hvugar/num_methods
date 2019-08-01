#ifndef INITIAL_VALUE_PROBLEM_H
#define INITIAL_VALUE_PROBLEM_H

#include "grid.h"

enum class InitialCondition
{
    InitialValue = 0,
    FirstDerivative = 1,
    SecondDerivative = 2
};

enum class FinalCondition
{
    FinalValue = 0,
    FinalFirstDerivative = 1,
    FinalSecondDerivative = 2
};


class MINIMUMSHARED_EXPORT InitialConditionODE
{
public:
    InitialCondition initialConditionType;
    double value;
};

class MINIMUMSHARED_EXPORT InitialConditionPDE {};

/**
 * @brief The InitialValueProblem class
 */
class MINIMUMSHARED_EXPORT InitialValueProblem
{
public:
    virtual ~InitialValueProblem();
};

class MINIMUMSHARED_EXPORT InitialValueProblemODE : public InitialValueProblem
{
public:
    virtual ~InitialValueProblemODE();
protected:
    virtual double initial(InitialCondition condition, unsigned int row = 1) const = 0;
};

class MINIMUMSHARED_EXPORT InitialValueProblemPDE : public InitialValueProblem
{
public:
    virtual ~InitialValueProblemPDE();
protected:
    virtual double initial(const SpaceNodePDE &sn, InitialCondition condition) const = 0;
};

/**
 * @brief The FinalValueProblem class
 */
class MINIMUMSHARED_EXPORT FinalValueProblem {};

class MINIMUMSHARED_EXPORT FinalValueProblemODE : public FinalValueProblem
{
public:
    virtual ~FinalValueProblemODE() = 0;
protected:
    virtual double final(FinalCondition condition, unsigned int row = 1) const = 0;
};

class MINIMUMSHARED_EXPORT FinalValueProblemPDE : public FinalValueProblem
{
public:
    virtual ~FinalValueProblemPDE() = 0;
protected:
    virtual double final(const SpaceNodePDE &sn, FinalCondition condition) const = 0;
};

#endif // INITIAL_VALUE_PROBLEM_H

#ifndef INITIAL_VALUE_PROBLEM_H
#define INITIAL_VALUE_PROBLEM_H

#include "grid.h"

enum class InitialCondition
{
    InitialValue = 0,
    FirstDerivative = 1//,
    //SecondDerivative = 2
};

class MINIMUMSHARED_EXPORT InitialConditionODE
{
public:
    InitialCondition initialConditionType;
    double value;
};

class MINIMUMSHARED_EXPORT InitialValueProblem {};

class MINIMUMSHARED_EXPORT InitialValueProblemODE : public InitialValueProblem
{
public:
    virtual ~InitialValueProblemODE();

protected:
    virtual auto initial(InitialCondition condition, unsigned int row = 1) const -> double = 0;
};

class MINIMUMSHARED_EXPORT InitialValueProblemPDE : public InitialValueProblem {};

#endif // INITIAL_VALUE_PROBLEM_H

#ifndef INITIAL_VALUE_PROBLEM_H
#define INITIAL_VALUE_PROBLEM_H

#include "grid.h"

enum class InitialCondition
{
    InitialValue = 0,
    InitialFirstDerivative = 1,
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
    PUBLIC_CONSTRUCTORS_VIRTUAL_DESTRUCTOR(InitialConditionODE);
public:
    InitialCondition initialConditionType;
    double value;
};

class MINIMUMSHARED_EXPORT InitialConditionPDE
{
    PUBLIC_CONSTRUCTORS_VIRTUAL_DESTRUCTOR(InitialConditionPDE);
};

class MINIMUMSHARED_EXPORT InitialValueProblem
{
    PUBLIC_CONSTRUCTORS_VIRTUAL_DESTRUCTOR(InitialValueProblem);
};

class MINIMUMSHARED_EXPORT FinalValueProblem
{
    PUBLIC_CONSTRUCTORS_VIRTUAL_DESTRUCTOR(FinalValueProblem);
};

class MINIMUMSHARED_EXPORT InitialValueProblemODE : public InitialValueProblem
{
    PUBLIC_CONSTRUCTORS_VIRTUAL_DESTRUCTOR(InitialValueProblemODE);

protected:
    virtual double initial(InitialCondition condition, unsigned int row = 1) const = 0;

    virtual void iterationInfo(double y, const PointNodeODE &node);
    virtual void iterationInfo(const DoubleVector &v, const PointNodeODE &node) const;
};

class MINIMUMSHARED_EXPORT FinalValueProblemODE : public FinalValueProblem
{
    PUBLIC_CONSTRUCTORS_VIRTUAL_DESTRUCTOR(FinalValueProblemODE);
protected:
    virtual double final(FinalCondition condition, unsigned int row = 1) const = 0;

    virtual void iterationInfo(double y, const PointNodeODE &node);
    virtual void iterationInfo(const DoubleVector &v, const PointNodeODE &node) const;
};

class MINIMUMSHARED_EXPORT InitialValueProblemPDE : public InitialValueProblem
{
    PUBLIC_CONSTRUCTORS_VIRTUAL_DESTRUCTOR(InitialValueProblemPDE);
protected:
    virtual double initial(const SpaceNodePDE &sn, InitialCondition condition) const = 0;
};

class MINIMUMSHARED_EXPORT FinalValueProblemPDE : public FinalValueProblem
{
    PUBLIC_CONSTRUCTORS_VIRTUAL_DESTRUCTOR(FinalValueProblemPDE);
protected:
    virtual double final(const SpaceNodePDE &sn, FinalCondition condition) const = 0;
};

#endif // INITIAL_VALUE_PROBLEM_H

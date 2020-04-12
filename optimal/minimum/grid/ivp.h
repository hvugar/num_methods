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

class MINIMUMSHARED_EXPORT InitialValueProblem
{
protected:
    PUBLIC_CONSTRUCTORS_VIRTUAL_DESTRUCTOR(InitialValueProblem);
};

class MINIMUMSHARED_EXPORT FinalValueProblem
{
protected:
    PUBLIC_CONSTRUCTORS_VIRTUAL_DESTRUCTOR(FinalValueProblem);
};

class MINIMUMSHARED_EXPORT InitialValueProblemODE : public InitialValueProblem
{
protected:
    PUBLIC_CONSTRUCTORS_VIRTUAL_DESTRUCTOR(InitialValueProblemODE);
    virtual auto initial(InitialCondition condition, unsigned int row = 1) const -> double = 0;
    virtual auto iterationInfo(double y, const PointNodeODE &node) const -> void;
    virtual auto iterationInfo(const DoubleVector &v, const PointNodeODE &node) const -> void;
};

class MINIMUMSHARED_EXPORT InitialValueProblemPDE : public InitialValueProblem
{
protected:
    PUBLIC_CONSTRUCTORS_VIRTUAL_DESTRUCTOR(InitialValueProblemPDE);
    virtual double initial(const SpaceNodePDE &sn, InitialCondition condition) const = 0;
};

class MINIMUMSHARED_EXPORT FinalValueProblemODE : public FinalValueProblem
{
protected:
    PUBLIC_CONSTRUCTORS_VIRTUAL_DESTRUCTOR(FinalValueProblemODE);
    virtual auto final(FinalCondition condition, unsigned int row = 1) const -> double = 0;
    virtual auto iterationInfo(double y, const PointNodeODE &node) const -> void;
    virtual auto iterationInfo(const DoubleVector &v, const PointNodeODE &node) const -> void;
};

class MINIMUMSHARED_EXPORT FinalValueProblemPDE : public FinalValueProblem
{
protected:
    PUBLIC_CONSTRUCTORS_VIRTUAL_DESTRUCTOR(FinalValueProblemPDE);
    virtual double final(const SpaceNodePDE &sn, FinalCondition condition) const = 0;
};

#endif // INITIAL_VALUE_PROBLEM_H

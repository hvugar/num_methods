#include "ivp.h"

//PUBLIC_CONSTRUCTORS_VIRTUAL_DESTRUCTOR_IMPL(InitialConditionODE);
//PUBLIC_CONSTRUCTORS_VIRTUAL_DESTRUCTOR_IMPL(InitialConditionPDE);
//PUBLIC_CONSTRUCTORS_VIRTUAL_DESTRUCTOR_IMPL(InitialValueProblem);

//**********************************************************************************************//

InitialConditionODE::InitialConditionODE() {}

InitialConditionODE::InitialConditionODE(const InitialConditionODE &) {}

InitialConditionODE& InitialConditionODE::operator=(const InitialConditionODE &other)
{
    if (this == &other) { return *this; }
    this->initialConditionType = other.initialConditionType;
    this->value = other.value;
    return *this;
}

InitialConditionODE::~InitialConditionODE() {}

//**********************************************************************************************//

InitialConditionPDE::InitialConditionPDE() {}

InitialConditionPDE::InitialConditionPDE(const InitialConditionPDE &) {}

InitialConditionPDE & InitialConditionPDE::operator=(const InitialConditionPDE &other)
{
    if (this == &other) { return *this; }
    return *this;
}

InitialConditionPDE::~InitialConditionPDE() {}

//**********************************************************************************************//

InitialValueProblem::InitialValueProblem() {}

InitialValueProblem::InitialValueProblem(const InitialValueProblem &) {}

InitialValueProblem& InitialValueProblem::operator=(const InitialValueProblem &other)
{
    if (this == &other) { return *this; }
    return *this;
}

InitialValueProblem::~InitialValueProblem() {}

//**********************************************************************************************//

FinalValueProblem::FinalValueProblem() {}

FinalValueProblem::FinalValueProblem(const FinalValueProblem &) {}

FinalValueProblem& FinalValueProblem::operator=(const FinalValueProblem &other)
{
    if (this == &other) { return *this; }
    return *this;
}

FinalValueProblem::~FinalValueProblem() {}

//**********************************************************************************************//

InitialValueProblemODE::InitialValueProblemODE() {}

InitialValueProblemODE::InitialValueProblemODE(const InitialValueProblemODE &) {}

InitialValueProblemODE& InitialValueProblemODE::operator=(const InitialValueProblemODE &other)
{
    if (this == &other) { return *this; }
    return *this;
}

InitialValueProblemODE::~InitialValueProblemODE() {}

void InitialValueProblemODE::iterationInfo(double y, const PointNodeODE &node) {}

void InitialValueProblemODE::iterationInfo(const DoubleVector &v, const PointNodeODE &node) const {}

//**********************************************************************************************//

InitialValueProblemPDE::InitialValueProblemPDE() {}

InitialValueProblemPDE::InitialValueProblemPDE(const InitialValueProblemPDE &) {}

InitialValueProblemPDE & InitialValueProblemPDE::operator=(const InitialValueProblemPDE &other)
{
    if (this == &other) { return *this; }
    return *this;
}

InitialValueProblemPDE::~InitialValueProblemPDE() {}

//**********************************************************************************************//

FinalValueProblemODE::FinalValueProblemODE() {}

FinalValueProblemODE::FinalValueProblemODE(const FinalValueProblemODE &) {}

FinalValueProblemODE& FinalValueProblemODE::operator=(const FinalValueProblemODE &other)
{
    if (this == &other) { return *this; }
    return *this;
}

FinalValueProblemODE::~FinalValueProblemODE() {}

void FinalValueProblemODE::iterationInfo(double, const PointNodeODE &) {}

void FinalValueProblemODE::iterationInfo(const DoubleVector &, const PointNodeODE &) const {}

//**********************************************************************************************//

FinalValueProblemPDE::FinalValueProblemPDE() {}

FinalValueProblemPDE::FinalValueProblemPDE(const FinalValueProblemPDE &) {}

FinalValueProblemPDE& FinalValueProblemPDE::operator =(const FinalValueProblemPDE &other)
{
    if (this == &other) { return *this; }
    return *this;
}

FinalValueProblemPDE::~FinalValueProblemPDE() {}

//**********************************************************************************************//

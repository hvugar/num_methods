#include "ivp.h"

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

InitialValueProblemODE::InitialValueProblemODE() {}

InitialValueProblemODE::InitialValueProblemODE(const InitialValueProblemODE &) {}

InitialValueProblemODE& InitialValueProblemODE::operator=(const InitialValueProblemODE &other)
{
    if (this == &other) { return *this; }
    return *this;
}

InitialValueProblemODE::~InitialValueProblemODE() {}

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

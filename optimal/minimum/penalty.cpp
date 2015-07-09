#include "penalty.h"

PenaltyMethod::PenaltyMethod()
{
}

PenaltyMethod::~PenaltyMethod()
{
}

const PRnFunctionList& PenaltyMethod::h() const
{
    return m_h;
}

const PRnFunctionList& PenaltyMethod::g() const
{
    return m_g;
}

void PenaltyMethod::calculate()
{

}


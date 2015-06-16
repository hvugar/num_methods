#include "penalty.h"

PenaltyMethod::PenaltyMethod()
{
}

PenaltyMethod::~PenaltyMethod()
{
}

std::vector<RnFunction*>& PenaltyMethod::h() const
{
    return m_h;
}

std::vector<RnFunction*>& PenaltyMethod::g() const
{
    return m_g;
}

void PenaltyMethod::calculate()
{

}


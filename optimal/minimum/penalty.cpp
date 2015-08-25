#include "penalty.h"

PenaltyMethod::PenaltyMethod()
{
}

PenaltyMethod::~PenaltyMethod()
{
}

const std::vector<RnFunction>& PenaltyMethod::h() const
{
    return m_h;
}

const std::vector<RnFunction>& PenaltyMethod::g() const
{
    return m_g;
}

void PenaltyMethod::calculate()
{

}

void PenaltyMethod::P(double x)
{

}


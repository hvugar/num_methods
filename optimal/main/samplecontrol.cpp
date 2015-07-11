#include "samplecontrol.h"
#include "cfunction1.h"
#include "cfunction2.h"

SampleControl::SampleControl() : SteepestDescentGradient()
{}

SampleControl::~SampleControl()
{
}

void SampleControl::calculateGradient()
{
    CFunction1* func = dynamic_cast<CFunction1*>(function());

    func->gradientJ(grad_step, m_g, m_x);

    printf("J[%d] = %.10f\n", iterationCount, func->fx(m_x));
}

void SampleControl::print()
{
//    SteepestDescentGradient::print();
}

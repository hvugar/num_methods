#include "samplecontrol.h"
#include "cfunction1.h"
#include "cfunction2.h"

SampleControl::SampleControl() : SteepestDescentGradient()
{}

SampleControl::~SampleControl()
{
}

void SampleControl::calcGradient()
{
    CFunction1* func = dynamic_cast<CFunction1*>(function());

    func->gradientJ(grad_step, mg, mx);

    printf("J[%d] = %.10f\n", k, func->fx(mx));
}

void SampleControl::print()
{
//    SteepestDescentGradient::print();
}

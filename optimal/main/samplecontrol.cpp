#include "samplecontrol.h"
#include "cfunction1.h"
#include "cfunction2.h"

SampleControl::SampleControl() : ConjugateGradient()
{}

SampleControl::~SampleControl()
{
}

void SampleControl::calcGradient()
{
    CFunction2* func = dynamic_cast<CFunction2*>(f());

    func->gradientJ(grad_step, mg, mx);

    printf("J[%d] = %.10f\n", k, func->fx(mx));
}

void SampleControl::print()
{}

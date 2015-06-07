#include "samplecontrol.h"
#include "cfunction.h"

SampleControl::SampleControl() : ConjugateGradient()
{}

SampleControl::~SampleControl()
{
}

void SampleControl::calcGradient()
{
    CFunction* func = dynamic_cast<CFunction*>(f());

    func->gradientJ(grad_step, mg, mx);

    //printX("gr", mg);
    printf("J[%d] = %.10f\n", mcount, func->fx(mx));
}

void SampleControl::print()
{}

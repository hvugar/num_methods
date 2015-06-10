#ifndef SAMPLECONTROL_H
#define SAMPLECONTROL_H

#include <sdgradient.h>
#include <cjtgradient.h>

class SampleControl : public ConjugateGradient
{
public:
    SampleControl();
    ~SampleControl();

    virtual void calcGradient();
    virtual void print();
};

#endif // SAMPLECONTROL_H

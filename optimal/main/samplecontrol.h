#ifndef SAMPLECONTROL_H
#define SAMPLECONTROL_H

#include <sdgradient.h>
#include <cjtgradient.h>

class SampleControl : public SteepestDescentGradient
{
public:
    SampleControl();
    ~SampleControl();

    virtual void calculateGradient();
    virtual void print();
};

#endif // SAMPLECONTROL_H

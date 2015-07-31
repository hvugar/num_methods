#ifndef SAMPLECONTROL_H
#define SAMPLECONTROL_H

#include <gradient_sd.h>
#include <gradient_cjt.h>

class SampleControl : public ConjugateGradient
{
public:
    SampleControl();
    ~SampleControl();

    virtual void calculateGradient();
    virtual void print();

    static void Main();
};



#endif // SAMPLECONTROL_H

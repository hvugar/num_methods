#ifndef NEWTONHEATEQUATIONEX1_H
#define NEWTONHEATEQUATIONEX1_H

#include <grid/pibvp.h>

#define SAMPLE_2

class MINIMUMSHARED_EXPORT NewtonHeatEquationEx1 : public NewtonHeatEquation
{
public:
    static void Main(int argc, char* argv[]);
    double U(const SpaceNodePDE &sn,const TimeNodePDE &tn) const;

protected:
    virtual double initial(const SpaceNodePDE &sn) const;
    virtual double boundary(const SpaceNodePDE &sn, const TimeNodePDE &tn) const;
    virtual double f(const SpaceNodePDE &sn, const TimeNodePDE &tn) const;
    virtual double a(const SpaceNodePDE &sn, const TimeNodePDE &tn) const;

    virtual double theta0(const TimeNodePDE &tn) const;
    virtual double theta1(const TimeNodePDE &tn) const;
    virtual double theta2(const TimeNodePDE &tn) const;
};

#endif // NEWTONHEATEQUATIONEX1_H

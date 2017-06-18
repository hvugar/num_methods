#ifndef NEWTONHEATEQUATIONEX1_H
#define NEWTONHEATEQUATIONEX1_H

#include <grid/nhpibvp.h>

#define SAMPLE_2

class MINIMUMSHARED_EXPORT NewtonHeatEquationEx1 : public NewtonHeatEquation
{
public:
    static void Main(int argc, char* argv[]);
    double U(const SpaceNode &sn,const TimeNode &tn) const;

protected:
    virtual double initial(const SpaceNode &sn) const;
    virtual double boundary(const SpaceNode &sn, const TimeNode &tn, BoundaryType boundary = Unused) const;
    virtual double f(const SpaceNode &sn, const TimeNode &tn) const;
    virtual double a(const SpaceNode &sn, const TimeNode &tn) const;

    virtual double theta0(const TimeNode &tn) const;
    virtual double theta1(const TimeNode &tn) const;
    virtual double theta2(const TimeNode &tn) const;
};

#endif // NEWTONHEATEQUATIONEX1_H

#ifndef IBACKWARDLOADEDHEATEQAUATION_H
#define IBACKWARDLOADEDHEATEQAUATION_H

#include <grid/pibvp.h>

class IBackwardLoadedHeatEqauation : public IParabolicIBVP
{
public:
    IBackwardLoadedHeatEqauation();

protected:
    virtual double initial(const SpaceNodePDE &sn) const = 0;
    virtual double boundary(const SpaceNodePDE &sn, const TimeNodePDE &tn, BoundaryType boundary = Unused) const = 0;
    virtual double f(const SpaceNodePDE &sn, const TimeNodePDE &tn) const = 0;

    //virtual double a(const SpaceNode &sn, const TimeNode &tn) const;
    virtual double g(const TimeNodePDE &tn) const = 0;
    virtual double h(const TimeNodePDE &tn) const = 0;
};

#endif // IBACKWARDLOADEDHEATEQAUATION_H

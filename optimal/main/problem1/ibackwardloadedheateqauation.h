#ifndef IBACKWARDLOADEDHEATEQAUATION_H
#define IBACKWARDLOADEDHEATEQAUATION_H

#include <grid/pibvp.h>

class IBackwardLoadedHeatEqauation : public IParabolicIBVP
{
public:
    IBackwardLoadedHeatEqauation();

protected:
    virtual double initial(const SpaceNode &sn) const = 0;
    virtual double boundary(const SpaceNode &sn, const TimeNode &tn, BoundaryType boundary = Unused) const = 0;
    virtual double f(const SpaceNode &sn, const TimeNode &tn) const = 0;

    //virtual double a(const SpaceNode &sn, const TimeNode &tn) const;
    virtual double g(const TimeNode &tn) const = 0;
    virtual double h(const TimeNode &tn) const = 0;
};

#endif // IBACKWARDLOADEDHEATEQAUATION_H

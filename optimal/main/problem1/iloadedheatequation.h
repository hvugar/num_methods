#ifndef ILOADEDHEATEQUATION_H
#define ILOADEDHEATEQUATION_H

#include <grid/pibvp.h>

class ILoadedHeatEquation : public IParabolicIBVP
{
public:
    struct Parameter {
        double k;
        double z;
        double e;
        unsigned int xi;
    } ;

    double a;
    double lambda0;
    double lambda1;
    double lambda2;
    double theta;

    void calculateM1(DoubleVector &u);
    void calculateM2(DoubleVector &u);

    unsigned int L;
    Parameter *params;

    virtual void layerInfo(const DoubleVector &, unsigned int) const {}
    virtual void layerInfo(const DoubleMatrix &, unsigned int) const {}

protected:
    virtual double initial(const SpaceNode &sn) const = 0;
    virtual double boundary(const SpaceNode &sn, const TimeNode &tn, BoundaryType boundary = Unused) const = 0;
    virtual double f(const SpaceNode &sn, const TimeNode &tn) const = 0;

    //virtual double a(const SpaceNode &sn, const TimeNode &tn) const;
    virtual double g(const TimeNode &tn) const = 0;
    virtual double h(const TimeNode &tn) const = 0;
};

#endif // ILOADEDHEATEQUATION_H

#ifndef HEATEQUATIONIBVP_1_H
#define HEATEQUATIONIBVP_1_H

#include <grid/hpibvp.h>
#include <printer.h>

class HeatEquationIBVP1 : public InitialBoundaryValueProblemPDE //: public HeatEquationIBVP
{
public:
    HeatEquationIBVP1();
    virtual ~HeatEquationIBVP1();

//protected:
//    virtual double initial(const SpaceNodePDE &) const;
    virtual double boundary(const SpaceNodePDE &, const TimeNodePDE &, BoundaryType boundary = Unused) const;
//    virtual double f(const SpaceNodePDE &, const TimeNodePDE &) const;

public:
    virtual void layerInfo(const DoubleVector &u, unsigned int ln);
//    virtual void layerInfo(const DoubleMatrix &u, unsigned int ln);

    void gridMethod1(DoubleVector &u, double a = 1.0);
};

#endif // HEATEQUATIONIBVP_1_H

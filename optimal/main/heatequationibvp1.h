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


class HeatEquationIBVP2D1 : public HeatEquationIBVP2D
{
public:
    HeatEquationIBVP2D1();
    virtual ~HeatEquationIBVP2D1();

protected:
    virtual double initial(const SpaceNodePDE &sn) const;
    virtual double boundary(const SpaceNodePDE &sn, const TimeNodePDE &tn, BoundaryType boundary = Unused) const;
    virtual double f(const SpaceNodePDE &sn, const TimeNodePDE &tn) const;

    virtual double env0(const SpaceNodePDE &sn, const TimeNodePDE &tn) const;
    virtual double env1(const SpaceNodePDE &sn, const TimeNodePDE &tn) const;
public:
    double a = 1.0;
    double alpha = 0.0;
    double lambda = 0.001;
};

#endif // HEATEQUATIONIBVP_1_H

#ifndef ISYSTEMLINEARODENONLOCALCONTIONSV2_H
#define ISYSTEMLINEARODENONLOCALCONTIONSV2_H

#include "islodenlcs.h"

class MINIMUMSHARED_EXPORT ISystemLinearODENonLocalContionsV2 : public ISystemLinearODENonLocalContions, public SystemDifferentialEquation
{
public:
    virtual double A(TimeNodePDE node, unsigned int row = 0, unsigned int col = 0) const = 0;
    virtual double B(TimeNodePDE node, unsigned int s, unsigned int row = 0, unsigned int col = 0) const = 0;
    virtual double C(TimeNodePDE node, unsigned int row = 0) const = 0;

    //virtual void calculateForward();
    void calculateBackward();

    void addCondition(const Condition &nlsc);
    void addLoadPoint(const LoadPoint &lpnt);
    void setRightSize(const DoubleVector &gamma);

    std::vector<Condition> nlscs;
    std::vector<LoadPoint> lpnts;
    DoubleVector gamma;
protected:
    void calculateCauchyProblem(const Condition &sc, const Condition &ec, const DoubleVector &x0, std::vector<DoubleVector> &rx, double h);
};

#endif // ISYSTEMLINEARODENONLOCALCONTIONSV2_H

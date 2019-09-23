#ifndef NEWTONHEATEQUATIONEX1_H
#define NEWTONHEATEQUATIONEX1_H

#include "border_global.h"

#define SAMPLE_2

class BORDERSHARED_EXPORT NewtonHeatEquationEx1 : public NewtonHeatEquation
{
public:
    static void Main(int argc, char* argv[]);
    double U(const SpaceNodePDE &sn,const TimeNodePDE &tn) const;

protected:
    virtual double initial(const SpaceNodePDE &sn, InitialCondition condition) const;
    virtual double boundary(const SpaceNodePDE &sn, const TimeNodePDE &tn, BoundaryConditionPDE &condition) const;
    virtual double f(const SpaceNodePDE &sn, const TimeNodePDE &tn) const;
    virtual double a(const SpaceNodePDE &sn, const TimeNodePDE &tn) const;

    virtual double theta0(const TimeNodePDE &tn) const;
    virtual double theta1(const TimeNodePDE &tn) const;
    virtual double theta2(const TimeNodePDE &tn) const;

protected:
    virtual const Dimension& timeDimension() const { return _timeDimension; }
    virtual const Dimension& spaceDimensionX() const { return _spaceDimensionX; }
    virtual const Dimension& spaceDimensionY() const { return _spaceDimensionY; }
    virtual const Dimension& spaceDimensionZ() const { return _spaceDimensionZ; }

    void setTimeDimension(const Dimension &timeDimension) { this->_timeDimension = timeDimension; }
    void setSpaceDimensionX(const Dimension &spaceDimensionX) { this->_spaceDimensionX = spaceDimensionX; }
    void setSpaceDimensionY(const Dimension &spaceDimensionY) { this->_spaceDimensionY = spaceDimensionY; }

private:
    Dimension _timeDimension;
    Dimension _spaceDimensionX;
    Dimension _spaceDimensionY;
    Dimension _spaceDimensionZ;
};

#endif // NEWTONHEATEQUATIONEX1_H

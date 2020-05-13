#ifndef HEATEQUATIONIBVP_H
#define HEATEQUATIONIBVP_H

#include "border_global.h"

class BORDERSHARED_EXPORT HeatEquationIBVP : public IHeatEquationIBVP
{
public:
    static void Main(int argc, char *argv[]);

    HeatEquationIBVP();

protected:
    virtual double initial(const SpaceNodePDE &sn, InitialCondition condition) const;
    virtual double boundary(const SpaceNodePDE &sn, const TimeNodePDE &tn, BoundaryConditionPDE &condition) const;
    virtual double f(const SpaceNodePDE &sn, const TimeNodePDE &tn) const;

    virtual void layerInfo(const DoubleVector&, const TimeNodePDE&) const;
    virtual void layerInfo(const DoubleMatrix&, const TimeNodePDE&) const;

protected:
    double weight() const { return 0.5; }

protected:
    virtual Dimension timeDimension() const { return _timeDimension; }
    virtual Dimension spaceDimensionX() const { return _spaceDimensionX; }
    virtual Dimension spaceDimensionY() const { return _spaceDimensionY; }
    virtual Dimension spaceDimensionZ() const { return _spaceDimensionZ; }

    void setTimeDimension(const Dimension &timeDimension) { this->_timeDimension = timeDimension; }
    void setSpaceDimensionX(const Dimension &spaceDimensionX) { this->_spaceDimensionX = spaceDimensionX; }
    void setSpaceDimensionY(const Dimension &spaceDimensionY) { this->_spaceDimensionY = spaceDimensionY; }

private:
    Dimension _timeDimension;
    Dimension _spaceDimensionX;
    Dimension _spaceDimensionY;
    Dimension _spaceDimensionZ;
    DeltaGrid2D deltaGrid;
};

class BORDERSHARED_EXPORT HeatEquationFBVP : public IHeatEquationFBVP
{
public:
    static void Main(int argc, char *argv[]);

    HeatEquationFBVP();

protected:
    virtual double final(const SpaceNodePDE &sn, FinalCondition condition) const;
    virtual double boundary(const SpaceNodePDE &sn, const TimeNodePDE &tn, BoundaryConditionPDE &condition) const;
    virtual double f(const SpaceNodePDE &sn, const TimeNodePDE &tn) const;

    virtual void layerInfo(const DoubleVector&, const TimeNodePDE&) const;
    virtual void layerInfo(const DoubleMatrix&, const TimeNodePDE&) const;

protected:
    double weight() const { return 0.5; }

protected:
    virtual Dimension timeDimension() const { return _timeDimension; }
    virtual Dimension spaceDimensionX() const { return _spaceDimensionX; }
    virtual Dimension spaceDimensionY() const { return _spaceDimensionY; }
    virtual Dimension spaceDimensionZ() const { return _spaceDimensionZ; }

    void setTimeDimension(const Dimension &timeDimension) { this->_timeDimension = timeDimension; }
    void setSpaceDimensionX(const Dimension &spaceDimensionX) { this->_spaceDimensionX = spaceDimensionX; }
    void setSpaceDimensionY(const Dimension &spaceDimensionY) { this->_spaceDimensionY = spaceDimensionY; }

private:
    Dimension _timeDimension;
    Dimension _spaceDimensionX;
    Dimension _spaceDimensionY;
    Dimension _spaceDimensionZ;
};

#endif // HEATEQUATIONIBVP_H

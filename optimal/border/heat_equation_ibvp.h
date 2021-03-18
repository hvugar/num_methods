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
    double weight() const { return 1.0; }

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

public:
    double _initialTemperature = 1.0;
    double _environmentTemperature = 2.0;
    double _lambda = 0.01;
    double _lambda0 = 0.001;
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

class BORDERSHARED_EXPORT LoadedHeatEquationIBVP : public ILoadedHeatEquationIBVP
{
public:
    static void Main(int argc, char *argv[]);

protected:
    virtual double initial(const SpaceNodePDE &sn, InitialCondition condition) const;
    virtual double boundary(const SpaceNodePDE &sn, const TimeNodePDE &tn, BoundaryConditionPDE &condition) const;
    virtual double f(const SpaceNodePDE &sn, const TimeNodePDE &tn) const;

    virtual void layerInfo(const DoubleVector&, const TimeNodePDE&) const;
    virtual void layerInfo(const DoubleMatrix&, const TimeNodePDE&) const;

protected:
    double weight() const { return 0.5; }

protected:
    virtual Dimension timeDimension() const { return Dimension(0.1, 0, 10); }
    virtual Dimension spaceDimensionX() const { return Dimension(0.1, 20, 30); }
    virtual Dimension spaceDimensionY() const { return Dimension(0.1, 30, 40); }
    virtual Dimension spaceDimensionZ() const { return Dimension(0.1, 0, 10); }
};

class BORDERSHARED_EXPORT LoadedHeatEquationFBVP : public ILoadedHeatEquationFBVP
{
public:
    static void Main(int argc, char *argv[]);

protected:
    virtual double final(const SpaceNodePDE &sn, FinalCondition condition) const;
    virtual double boundary(const SpaceNodePDE &sn, const TimeNodePDE &tn, BoundaryConditionPDE &condition) const;
    virtual double f(const SpaceNodePDE &sn, const TimeNodePDE &tn) const;

    virtual void layerInfo(const DoubleVector&, const TimeNodePDE&) const;
    virtual void layerInfo(const DoubleMatrix&, const TimeNodePDE&) const;

protected:
    double weight() const { return 0.5; }

protected:
    virtual Dimension timeDimension() const { return Dimension(0.01, 0, 100); }
    virtual Dimension spaceDimensionX() const { return Dimension(0.01, 0, 100); }
    virtual Dimension spaceDimensionY() const { return Dimension(0.01, 0, 100); }
    virtual Dimension spaceDimensionZ() const { return Dimension(0.01, 0, 100); }
};

#endif // HEATEQUATIONIBVP_H

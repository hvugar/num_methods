#ifndef SOLVER3_H
#define SOLVER3_H

#include "global.h"

namespace p3p3
{

class Functional;

class PROBLEM3P_SHARED_EXPORT HeatEquationIBVP : public IHeatEquationIBVP
{
public:
    HeatEquationIBVP(Functional*);
    HeatEquationIBVP(const HeatEquationIBVP &);
    HeatEquationIBVP & operator =(const HeatEquationIBVP &);
    virtual ~HeatEquationIBVP() override;

protected:
    virtual auto initial(const SpaceNodePDE &sn, InitialCondition ic) const -> double override;
    virtual auto boundary(const SpaceNodePDE &sn, const TimeNodePDE &tn, BoundaryConditionPDE &bc) const -> double override;
    virtual auto f(const SpaceNodePDE &sn, const TimeNodePDE &tn) const -> double override;

    virtual auto layerInfo(const DoubleMatrix &u, const TimeNodePDE &tn) const -> void override;
    virtual auto timeDimension() const -> Dimension override;
    virtual auto spaceDimensionX() const -> Dimension override;
    virtual auto spaceDimensionY() const -> Dimension override;
    virtual auto spaceDimensionZ() const -> Dimension override;

private:
    double _initial_temperature = 0.0;
    double _enviroment_temperature = 0.5;

    Functional* _functional;
};


class PROBLEM3P_SHARED_EXPORT HeatEquationFBVP : virtual public IHeatEquationFBVP
{
public:
    HeatEquationFBVP(Functional *function);
    HeatEquationFBVP(const HeatEquationFBVP &);
    HeatEquationFBVP & operator =(const HeatEquationFBVP &);
    virtual ~HeatEquationFBVP() override;

protected:
    virtual auto final(const SpaceNodePDE &sn, FinalCondition condition) const -> double override;
    virtual auto boundary(const SpaceNodePDE &sn, const TimeNodePDE &tn, BoundaryConditionPDE &condition) const -> double override;
    virtual auto f(const SpaceNodePDE &sn, const TimeNodePDE &tn) const -> double override;

public:
    virtual auto layerInfo(const DoubleMatrix &, const TimeNodePDE &) const -> void override;
    virtual auto timeDimension() const -> Dimension override;
    virtual auto spaceDimensionX() const -> Dimension override;
    virtual auto spaceDimensionY() const -> Dimension override;
    virtual auto spaceDimensionZ() const -> Dimension override;

public:
    size_t i;
    Functional* _functional;
};

class PROBLEM3P_SHARED_EXPORT Functional : public RnFunction, public IGradient
{
public:
    static void Main(int argc, char** argv);

    Functional();

    virtual double fx(const DoubleVector &x) const;
    auto integral(const DoubleMatrix &u) const -> double;
    virtual void gradient(const DoubleVector &x, DoubleVector &g) const;

    auto timeDimension() const -> Dimension { return Dimension(0.01, 0, 100) /*Dimension(0.0000005, 0, 20000000)*/; }

    auto spaceDimensionX() const -> Dimension { return Dimension(0.01, 0, 100); }

    auto spaceDimensionY() const -> Dimension { return Dimension(0.01, 0, 100); }

    auto spaceDimensionZ() const -> Dimension { return Dimension(0.01, 0, 100); }

    size_t heat_source_number = 2;
    size_t measure_point_number = 4;
    double _lambda1 = 0.01;

    DoubleMatrix U;
    DoubleMatrix V;

    DoubleVector *pp = nullptr;
    DoubleVector *px;
    DoubleVector *py;


    SpacePoint tr(const TimeNodePDE &tn, size_t i) const;
    double v(const TimeNodePDE &tn, size_t i) const;
    double q(const TimeNodePDE &tn, size_t i) const;

    HeatEquationIBVP *ih;
    HeatEquationFBVP *fh;

    DoubleVector x;

    void setVector(const DoubleVector &x) const;

    const double R[2] = {0.40, 0.20};

    void frw_saveToImage(const DoubleMatrix &u, const TimeNodePDE &tn) const;
};

};

#endif // SOLVER3_H

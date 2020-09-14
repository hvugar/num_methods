#ifndef SOLVER2_H
#define SOLVER2_H

#include "global.h"

namespace p3p0
{

class Solver2;

class Common
{
public:
    Common() {}
    Common(const Common &) {}
    Common & operator =(const Common &) { return *this; }
    virtual ~Common() {}

    virtual double lambda1() const { return _lambda1; }
    virtual double theta() const { return _theta; }
    virtual double mu(const SpaceNodePDE &/*sn*/) const { return 1.0; }

    auto convert(const DoubleVector &x, size_t size, size_t length) -> void;
    auto convert(size_t size, size_t length, DoubleVector &x) const -> void;

    double _lambda1 = +0.001;
    double _theta = +2.0;
    size_t i;

    DoubleVector V;
    DoubleVector U;

    double **q = nullptr;
};

class PROBLEM3P_SHARED_EXPORT HeatEquationIBVP : virtual public IHeatEquationIBVP
{
public:
    HeatEquationIBVP();
    HeatEquationIBVP(const HeatEquationIBVP &);
    HeatEquationIBVP & operator =(const HeatEquationIBVP &);
    virtual ~HeatEquationIBVP() override;

public:
    auto q(size_t i, const TimeNodePDE &tn) const -> double;
    auto v(size_t i, const PointNodeODE &tn, SpacePoint &vl) const -> void;
    auto z(size_t i, const TimeNodePDE &tn) const -> double;

protected:
    virtual auto initial(const SpaceNodePDE &sn, InitialCondition ic) const -> double override;
    virtual auto boundary(const SpaceNodePDE &sn, const TimeNodePDE &tn, BoundaryConditionPDE &bc) const -> double override;
    virtual auto f(const SpaceNodePDE &sn, const TimeNodePDE &tn) const -> double override;

    virtual auto layerInfo(const DoubleVector &u, const TimeNodePDE &tn) const -> void override;
    virtual auto timeDimension() const -> Dimension override;
    virtual auto spaceDimensionX() const -> Dimension override;
    virtual auto spaceDimensionY() const -> Dimension override;
    virtual auto spaceDimensionZ() const -> Dimension override;

public:
    Solver2 *s;
};

class PROBLEM3P_SHARED_EXPORT HeatEquationFBVP : virtual public IHeatEquationFBVP
{
public:
    HeatEquationFBVP();
    HeatEquationFBVP(const HeatEquationFBVP &);
    HeatEquationFBVP & operator =(const HeatEquationFBVP &);
    virtual ~HeatEquationFBVP() override;

protected:
    virtual auto final(const SpaceNodePDE &sn, FinalCondition ic) const -> double override;
    virtual auto boundary(const SpaceNodePDE &sn, const TimeNodePDE &tn, BoundaryConditionPDE &bc) const -> double override;
    virtual auto f(const SpaceNodePDE &sn, const TimeNodePDE &tn) const -> double override;

    virtual auto layerInfo(const DoubleVector &u, const TimeNodePDE &tn) const -> void override;
    virtual auto timeDimension() const -> Dimension override;
    virtual auto spaceDimensionX() const -> Dimension override;
    virtual auto spaceDimensionY() const -> Dimension override;
    virtual auto spaceDimensionZ() const -> Dimension override;

public:
    Solver2 *s;
};

class PROBLEM3P_SHARED_EXPORT Solver2 : virtual public RnFunction,
                                        virtual public IGradient,
                                        virtual public Common
{
public:
    static void Main(int argc, char** argv);

public:
    Solver2();

    inline auto timeDimension() const -> Dimension { return Dimension(0.001, 0, 1000); }
    inline auto spaceDimensionX() const -> Dimension { return Dimension(0.01, 0, 100); }

protected:
    virtual auto gradient(const DoubleVector &x, DoubleVector &g) const -> void;
    virtual auto fx(const DoubleVector &x) const -> double;
    auto integral(const DoubleMatrix &) const -> double;

    HeatEquationIBVP forward;
    HeatEquationFBVP backward;
};

}

#endif // SOLVER2_H

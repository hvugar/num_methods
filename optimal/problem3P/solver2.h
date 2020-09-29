#ifndef SOLVER2_H
#define SOLVER2_H

#include "global.h"

#define OPTIMIZE_Q
//#define OPTIMIZE_Q_FB

namespace p3p0
{

class CommonParameter
{
public:
    CommonParameter() { mp = nullptr; }
    CommonParameter(const CommonParameter &) {}
    CommonParameter & operator =(const CommonParameter &) { return *this; }
    virtual ~CommonParameter();

    inline auto lambda1() const -> double { return _lambda1; }
    inline auto theta() const -> double { return _theta; }
    inline auto mu(const SpaceNodePDE &/*sn*/) const -> double { return 1.0; }
    inline auto timeDimension() const -> Dimension { return _timeDimension; }
    inline auto spaceDimensionX() const -> Dimension { return _spaceDimensionX; }
    virtual auto setTimeDimension(const Dimension &timeDimension) -> void;
    virtual auto setSpaceDimensionX(const Dimension &spaceDimensionX) -> void;
    virtual auto setControlSize(size_t heatSourceNumber, size_t measrPointNumber) -> void;

    auto convertFromVector(const DoubleVector &x) -> void;
    auto convertToVector(DoubleVector &x) const -> void;
    auto qNorm1(double t) const -> DoubleVector;
    inline auto sqr(double x) const -> double { return x*x; }

public:
    auto q(const TimeNodePDE &tn) const -> DoubleVector;
    auto z(const TimeNodePDE &tn) const -> DoubleVector;

    double _lambda1 = +0.001;
    double _theta = +2.0;
    double initialTemperature = 0.0;

    DoubleVector V;
    DoubleVector U;

    DoubleVector *mq = nullptr;
    DoubleVector *mp = nullptr;
    DoubleVector *mz = nullptr;

    size_t heatSourceNumber = 2;

#ifdef OPTIMIZE_Q_FB
    DoubleMatrix alpha;
    DoubleMatrix betta;
    DoubleMatrix nomnU;
    DoubleVector measurePoints;
    size_t measrPointNumber = 4;
#endif

    double epsilon = 0.0001;
    double no_norm = 0.0000;

    unsigned int _w = 12;
    unsigned int _p = 6;

protected:
    Dimension _timeDimension;
    Dimension _spaceDimensionX;
};

class PROBLEM3P_SHARED_EXPORT HeatEquationIBVP : virtual public IHeatEquationIBVP
{
public:
    HeatEquationIBVP();
    HeatEquationIBVP(const HeatEquationIBVP &);
    HeatEquationIBVP & operator =(const HeatEquationIBVP &);
    virtual ~HeatEquationIBVP() override;

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
    CommonParameter *common;
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

    virtual auto layerInfo(const DoubleVector &psi, const TimeNodePDE &tn) const -> void override;
    virtual auto timeDimension() const -> Dimension override;
    virtual auto spaceDimensionX() const -> Dimension override;
    virtual auto spaceDimensionY() const -> Dimension override;
    virtual auto spaceDimensionZ() const -> Dimension override;

public:
    CommonParameter *common;
};

class PROBLEM3P_SHARED_EXPORT Functional : virtual public RnFunction, virtual public IGradient,
        virtual public IPrinter, virtual public CommonParameter
{
public:
    static void Main(int argc, char** argv);

    Functional(double thermalDiffusivity, double thermalConvection, double thermalConductivity);

protected:

    virtual auto fx(const DoubleVector &x) const -> double;
    auto integral(const DoubleVector &) const -> double;
    auto norm(const DoubleVector &) const -> double;

    virtual auto gradient(const DoubleVector &x, DoubleVector &g) const -> void;

    virtual auto print(unsigned int iteration, const DoubleVector &x, const DoubleVector &g, double f, double alpha, GradientBasedMethod::MethodResult result) const -> void;

    HeatEquationIBVP forward;
    HeatEquationFBVP backward;
};

}

#endif // SOLVER2_H

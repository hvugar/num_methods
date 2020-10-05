#ifndef SOLVER2_H
#define SOLVER2_H

#include "global.h"

#define OPTIMIZE_Y
//#define OPTIMIZE_ETA

namespace p3p0
{

class CommonParameter
{
public:
    CommonParameter() { mp = nullptr; }
    CommonParameter(const CommonParameter &) {}
    CommonParameter & operator =(const CommonParameter &) { return *this; }
    virtual ~CommonParameter();

    inline auto lambda0() const -> double { return _lambda0; }
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

    auto g0(size_t i, size_t ln) const -> double;
    auto gi(size_t i, size_t ln) const -> double;
    auto gp(size_t i, size_t ln) const -> double;

    double _lambda0 = +0.0001;
    double _lambda1 = +0.001;
    double _theta = +2.0;
    double initialTemperature = 0.0;

    DoubleVector V;
    DoubleVector U;

    DoubleVector *mq = nullptr;
    DoubleVector *mp = nullptr;
    DoubleVector *mz = nullptr;

    size_t heatSourceNumber = 2;

#ifdef OPTIMIZE_Y
    DoubleMatrix alpha;
    DoubleMatrix betta;
    DoubleMatrix omega;
    DoubleVector mPnts;
    size_t measrPointNumber = 4;
    DoubleVector *uv = nullptr;
    DoubleVector *ud = nullptr;

    DoubleMatrix alphaN;
    DoubleMatrix bettaN;
    DoubleMatrix omegaN;
    DoubleVector mPntsN;
#endif

    double epsilon = 0.0000;
    double no_norm = 1.0000;

    DoubleVector *qMin = nullptr;
    DoubleVector *qMax = nullptr;

    unsigned int _w = 8;
    unsigned int _p = 4;

    double R = 0.0;

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
        virtual public IPrinter, virtual public CommonParameter, virtual public IProjection
{
public:
    static void Main(int argc, char** argv);

    Functional(double thermalDiffusivity, double thermalConductivity, double thermalConvection);

protected:

    virtual auto fx(const DoubleVector &x) const -> double;
    auto integral(const DoubleVector &x) const -> double;
    auto norm(const DoubleVector &x) const -> double;
    auto penalty(const DoubleVector &x) const -> double;

    virtual auto gradient(const DoubleVector &x, DoubleVector &g) const -> void;

    virtual auto print(unsigned int iteration, const DoubleVector &x, const DoubleVector &g, double f, double alpha, GradientBasedMethod::MethodResult result) const -> void;

    virtual void project(DoubleVector &x, size_t index);

    HeatEquationIBVP forward;
    HeatEquationFBVP backward;
};

}

#endif // SOLVER2_H

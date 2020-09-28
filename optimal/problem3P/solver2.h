#ifndef SOLVER2_H
#define SOLVER2_H

#include "global.h"

#define OPTIMIZE_Q

namespace p3p0
{

class Solver2;

class Common
{
public:
    Common() {}
    Common(const Common &) {}
    Common & operator =(const Common &) { return *this; }
    virtual ~Common();

    virtual double lambda1() const { return _lambda1; }
    virtual double theta() const { return _theta; }
    virtual double mu(const SpaceNodePDE &/*sn*/) const { return 1.0; }

    auto convert1(const DoubleVector &x, size_t size, size_t length) -> void;
    auto convert2(size_t size, size_t length, DoubleVector &x) const -> void;

    double _lambda1 = +0.0001;
    double _theta = +2.0;
    size_t i;

    DoubleVector V;
    DoubleVector U;

    DoubleVector *q = nullptr;
    DoubleVector *p = nullptr;
    size_t heatSourceNumber = 2;
    size_t measrPointNumber = 4;

    DoubleMatrix alpha;
    DoubleMatrix betta;
    DoubleMatrix nomnU;
    DoubleVector measurePoints;
};

class PROBLEM3P_SHARED_EXPORT HeatEquationIBVP : virtual public IHeatEquationIBVP
{
public:
    HeatEquationIBVP();
    HeatEquationIBVP(const HeatEquationIBVP &);
    HeatEquationIBVP & operator =(const HeatEquationIBVP &);
    virtual ~HeatEquationIBVP() override;

public:
    auto q(const TimeNodePDE &tn) const -> DoubleVector;
    auto z(const TimeNodePDE &tn) const -> DoubleVector;
    auto v(size_t i, const PointNodeODE &tn, SpacePoint &vl) const -> void;

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

    virtual auto layerInfo(const DoubleVector &psi, const TimeNodePDE &tn) const -> void override;
    virtual auto timeDimension() const -> Dimension override;
    virtual auto spaceDimensionX() const -> Dimension override;
    virtual auto spaceDimensionY() const -> Dimension override;
    virtual auto spaceDimensionZ() const -> Dimension override;

public:
    Solver2 *s;
    HeatEquationIBVP *h;
};

class PROBLEM3P_SHARED_EXPORT Solver2 : virtual public RnFunction, virtual public IGradient, virtual public IPrinter, virtual public Common
{
public:
    static void Main(int argc, char** argv);

public:
    Solver2();

    inline auto timeDimension() const -> Dimension { return _timeDimension; }
    inline auto spaceDimensionX() const -> Dimension { return _spaceDimensionX; }

    auto qNorm1(double t) const -> DoubleVector;

protected:

    virtual auto fx(const DoubleVector &x) const -> double;
    auto integral(const DoubleVector &) const -> double;
    auto norm(const DoubleVector &) const -> double;

    virtual auto gradient(const DoubleVector &x, DoubleVector &g) const -> void;

    virtual auto print(unsigned int iteration, const DoubleVector &x, const DoubleVector &g, double f, double alpha, GradientBasedMethod::MethodResult result) const -> void;

    HeatEquationIBVP forward;
    HeatEquationFBVP backward;

private:
    Dimension _timeDimension;
    Dimension _spaceDimensionX;
    double epsilon = 0.01;
    double no_norm = 0.00;

    unsigned int _w = 12;
    unsigned int _p = 6;
};

}

/********************************************************************************************************************************/

namespace p3p2
{

class CommonParameters
{
public:
    virtual ~CommonParameters();

    auto q(const TimeNodePDE &tn) const -> DoubleVector;
    auto z(const TimeNodePDE &tn) const -> DoubleVector;

    virtual double lambda1() const { return m_lambda1; }
    virtual double theta() const { return m_theta; }
    virtual double mu(const SpaceNodePDE &/*sn*/) const { return 1.0; }

    auto timeDimension() const -> const Dimension& { return _timeDimension; }
    auto spaceDimensionX() const -> const Dimension& { return _spaceDimensionX; }

    virtual void setTimeDimension(const Dimension &timeDimension)
    {
        _timeDimension = timeDimension;
        const auto time_size = _timeDimension.size();

        if (mq != nullptr)
        {
            for (size_t ln=0; ln<time_size; ln++) { mq[ln].clear(); }
            delete [] mq;
            mq = nullptr;
        }
        mq = new DoubleVector[time_size];
        for (size_t ln=0; ln<time_size; ln++) { mq[ln].resize(heatSourceNumber); }

        if (mz != nullptr)
        {
            for (size_t ln=0; ln<time_size; ln++) { mz[ln].clear(); }
            delete [] mz;
            mz = nullptr;
        }
        mz = new DoubleVector[time_size];
        for (size_t ln=0; ln<time_size; ln++) { mz[ln].resize(heatSourceNumber); }

        if (mp != nullptr)
        {
            for (size_t ln=0; ln<time_size; ln++) { mp[ln].clear(); }
            delete [] mp;
            mp = nullptr;
        }
        mp = new DoubleVector[time_size];
        for (size_t ln=0; ln<time_size; ln++) { mp[ln].resize(heatSourceNumber); }

    }

    virtual void setSpaceDimensionX(const Dimension &spaceDimensionX)
    {
        _spaceDimensionX = spaceDimensionX;
        const auto spaceX_size = _spaceDimensionX.size();

        U.resize(spaceX_size, 0.0);
        V.resize(spaceX_size, 0.0);
    }

    auto convert1(const DoubleVector &x, size_t size, size_t length) -> void;
    auto convert2(size_t size, size_t length, DoubleVector &x) const -> void;

    auto qNorm1(double t) const -> DoubleVector;
    inline auto sqr(double x) const -> double { return x*x; }

public:
    Dimension _timeDimension;
    Dimension _spaceDimensionX;

    DoubleVector V;
    DoubleVector U;

    size_t heatSourceNumber = 2;
    size_t measrPointNumber = 4;

    DoubleMatrix alpha;
    DoubleMatrix betta;
    DoubleMatrix nomnU;
    double* measurePoints;

    DoubleVector *mq = nullptr;
    DoubleVector *mz = nullptr;
    DoubleVector *mp = nullptr;

    double m_a = 0.01;
    double m_lambda1 = +0.0001;
    double m_theta = +2.0;
    double m_b = 0.0;

    double epsilon = 0.01;
    double no_norm = 0.00;

    const unsigned int _w = 12;
    const unsigned int _p = 6;
};

class HeatEquationIBVP1 : virtual public IHeatEquationIBVP, virtual public CommonParameters
{
public:
    virtual ~HeatEquationIBVP1() override;

    auto calculate_forward(const DoubleVector &x) const -> void;

protected:
    virtual auto initial(const SpaceNodePDE &sn, InitialCondition ic) const -> double override;
    virtual auto boundary(const SpaceNodePDE &sn, const TimeNodePDE &tn, BoundaryConditionPDE &bc) const -> double override;
    virtual auto f(const SpaceNodePDE &sn, const TimeNodePDE &tn) const -> double override;

    virtual auto layerInfo(const DoubleVector &u, const TimeNodePDE &tn) const -> void override;

protected:
    virtual auto timeDimension() const -> Dimension override;
    virtual auto spaceDimensionX() const -> Dimension override;
    virtual auto spaceDimensionY() const -> Dimension override;
    virtual auto spaceDimensionZ() const -> Dimension override;
};

class HeatEquationFBVP1 : virtual public IHeatEquationFBVP, virtual public HeatEquationIBVP1
{
public:
    virtual ~HeatEquationFBVP1() override;

    auto calculate_backward(const DoubleVector &x) const -> void;

protected:
    virtual auto final(const SpaceNodePDE &sn, FinalCondition ic) const -> double override;
    virtual auto boundary(const SpaceNodePDE &sn, const TimeNodePDE &tn, BoundaryConditionPDE &bc) const -> double override;
    virtual auto f(const SpaceNodePDE &sn, const TimeNodePDE &tn) const -> double override;

    virtual auto layerInfo(const DoubleVector &psi, const TimeNodePDE &tn) const -> void override;
    virtual auto timeDimension() const -> Dimension override;
    virtual auto spaceDimensionX() const -> Dimension override;
    virtual auto spaceDimensionY() const -> Dimension override;
    virtual auto spaceDimensionZ() const -> Dimension override;
};

class PROBLEM3P_SHARED_EXPORT Functional : virtual public RnFunction, virtual public IGradient,
        virtual public IPrinter, virtual public HeatEquationFBVP1
{
public:
    static void Main(int argc, char** argv);

public:
    Functional(double thermalDiffusivity, double thermalConvection, double thermalConductivity);
    virtual ~Functional() override;

    virtual auto fx(const DoubleVector &x) const -> double override;
    auto integral(const DoubleVector &) const -> double;
    auto norm(const DoubleVector &) const -> double;

    virtual auto gradient(const DoubleVector &x, DoubleVector &g) const -> void override;

    virtual auto print(unsigned int iteration, const DoubleVector &x, const DoubleVector &g, double f, double alpha,
                       GradientBasedMethod::MethodResult result) const -> void override;
};

}

#endif // SOLVER2_H

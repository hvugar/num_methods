#ifndef SOLVER3_H
#define SOLVER3_H

#include "global.h"

namespace p3p3
{

class Functional;

class PROBLEM3P_SHARED_EXPORT LoadedHeatEquationIBVP : public ILoadedHeatEquationIBVP
{
public:
    LoadedHeatEquationIBVP(Functional*);
    LoadedHeatEquationIBVP(const LoadedHeatEquationIBVP &);
    LoadedHeatEquationIBVP & operator =(const LoadedHeatEquationIBVP &);
    virtual ~LoadedHeatEquationIBVP() override;

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

    double _fx;
};


class PROBLEM3P_SHARED_EXPORT LoadedHeatEquationFBVP : virtual public ILoadedHeatEquationFBVP
{
public:
    LoadedHeatEquationFBVP(Functional *function);
    LoadedHeatEquationFBVP(const LoadedHeatEquationFBVP &);
    LoadedHeatEquationFBVP & operator =(const LoadedHeatEquationFBVP &);
    virtual ~LoadedHeatEquationFBVP() override;

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
    virtual auto integral(const DoubleMatrix &u) const -> double;
    virtual void gradient(const DoubleVector &x, DoubleVector &g) const;

    auto timeDimension() const -> Dimension { return Dimension(0.01, 0, 1000) /*Dimension(0.0000005, 0, 20000000)*/; }

    auto spaceDimensionX() const -> Dimension { return Dimension(0.01, 0, 100); }

    auto spaceDimensionY() const -> Dimension { return Dimension(0.01, 0, 100); }

    auto spaceDimensionZ() const -> Dimension { return Dimension(0.01, 0, 100); }

    size_t heating_point_number = 2;
    size_t measure_point_number = 4;
    double _lambda1 = 0.01;

    DoubleMatrix k;
    DoubleMatrix z;
    DoubleVector s;
    std::vector<SpacePoint> measure_point;

    DoubleMatrix U;
    DoubleMatrix V;

    DoubleVector *pp = nullptr;
    DoubleVector *px = nullptr;
    DoubleVector *py = nullptr;

    DoubleVector *uu = nullptr;
    DoubleVector *ux = nullptr;
    DoubleVector *uy = nullptr;


    SpacePoint tr(const TimeNodePDE &tn, size_t i) const;
    double v(const TimeNodePDE &tn, size_t i) const;
    double q(const TimeNodePDE &tn, size_t i) const;

    LoadedHeatEquationIBVP *ih;
    LoadedHeatEquationFBVP *fh;

    void setVector(const DoubleVector &x) const;

    const double R[2] = {0.40, 0.20};

    void frw_saveToImage(const DoubleMatrix &u, const TimeNodePDE &tn) const;
};

};

#endif // SOLVER3_H

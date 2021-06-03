#ifndef SOLVER6_H
#define SOLVER6_H

#include "global.h"

#define P3P6_TIME_STEP 0.01
#define P3P6_TIME_MAX 1000
#define P3P6_TIME_SIZE P3P6_TIME_MAX+1
#define P3P6_DIMX_STEP 0.01
#define P3P6_DIMX_MAX 100
#define P3P6_NUM_GRAD_STEP 0.01

#define P3P6_OPTIMIZE_Q
//#define P3P6_OPTIMIZE_Y
#define P3P6_CALCULATE_GRAD

namespace p3p6
{

class Functional;

class Shared
{
public:
    auto q(const TimeNodePDE &tn, size_t i) const -> double;
    auto s(const TimeNodePDE &tn, size_t i) const -> double;

public:
    DoubleVector U = DoubleVector(P3P6_TIME_SIZE, 10.0);
    DoubleVector V = DoubleVector(P3P6_TIME_SIZE, 10.0);

    DoubleVector uj_v;
    DoubleVector pi_v;

    DoubleVector qi_v;
    DoubleVector gj_v;

#ifdef P3P6_OPTIMIZE_Y
    DoubleMatrix k;
    DoubleMatrix z;
#endif

protected:
    Dimension _timeDimension   = Dimension(P3P6_TIME_STEP, 0, P3P6_TIME_MAX);
    Dimension _spaceDimensionX = Dimension(P3P6_DIMX_STEP, 0, P3P6_DIMX_MAX);

public:
    const size_t heating_source_number = 2;
    const size_t meausere_point_number = 4;
    SpacePoint *measurePoint;
    double lambda1 = 0.01;
    double initial_temperature = 0.0;
    double envrmnt_temperature = 5.0;

    size_t drawImages = 0;
};

class PROBLEM3P_SHARED_EXPORT LoadedHeatEquationIBVP : public IHeatEquationIBVP, public virtual Shared /*public ILoadedHeatEquationIBVP*/
{
public:
    LoadedHeatEquationIBVP(double a = 1.0, double lambda0 = 0.0, double lambda1 = 0.0);
    LoadedHeatEquationIBVP(const LoadedHeatEquationIBVP &);
    LoadedHeatEquationIBVP & operator =(const LoadedHeatEquationIBVP &);
    virtual ~LoadedHeatEquationIBVP() override;

protected:
    virtual auto initial(const SpaceNodePDE &sn, InitialCondition ic) const -> double override;
    virtual auto boundary(const SpaceNodePDE &sn, const TimeNodePDE &tn, BoundaryConditionPDE &bc) const -> double override;
    virtual auto f(const SpaceNodePDE &sn, const TimeNodePDE &tn) const -> double override;

    virtual auto layerInfo(const DoubleVector &u, const TimeNodePDE &tn) const -> void override;
    virtual auto timeDimension() const -> Dimension override;
    virtual auto spaceDimensionX() const -> Dimension override;
    virtual auto spaceDimensionY() const -> Dimension override;
    virtual auto spaceDimensionZ() const -> Dimension override;

    auto saveImage(const DoubleVector &u, const TimeNodePDE &tn) const -> void;
};

class PROBLEM3P_SHARED_EXPORT LoadedHeatEquationFBVP : public IHeatEquationFBVP, public virtual Shared /*public ILoadedHeatEquationFBVP*/
{
public:
    LoadedHeatEquationFBVP(double a = -1.0, double lambda0 = 0.0, double lambda1 = 0.0);
    LoadedHeatEquationFBVP(const LoadedHeatEquationFBVP &);
    LoadedHeatEquationFBVP & operator =(const LoadedHeatEquationFBVP &);
    virtual ~LoadedHeatEquationFBVP() override;

protected:
    virtual auto final(const SpaceNodePDE &sn, FinalCondition fc) const -> double override;
    virtual auto boundary(const SpaceNodePDE &sn, const TimeNodePDE &tn, BoundaryConditionPDE &bc) const -> double override;
    virtual auto f(const SpaceNodePDE &sn, const TimeNodePDE &tn) const -> double override;

    virtual auto layerInfo(const DoubleVector &u, const TimeNodePDE &tn) const -> void override;
    virtual auto timeDimension() const -> Dimension override;
    virtual auto spaceDimensionX() const -> Dimension override;
    virtual auto spaceDimensionY() const -> Dimension override;
    virtual auto spaceDimensionZ() const -> Dimension override;

    auto saveImage(const DoubleVector &u, const TimeNodePDE &tn) const -> void;
};

class PROBLEM3P_SHARED_EXPORT Functional :
        public LoadedHeatEquationIBVP,
        public LoadedHeatEquationFBVP,
        public RnFunction,
        public IGradient,
        public IProjection,
        public IPrinter
{
public:
    static void Main(int argc, char** argv);

    Functional(double diffusivity = 1.0, double conductivity = 0.0, double convection = 0.0, double lambda = 0.0);

    virtual auto fx(const DoubleVector &x) const -> double override;
    virtual auto integral(const DoubleMatrix &U) const -> double;
    virtual auto gradient(const DoubleVector &x, DoubleVector &g) const -> void override;

    virtual auto project(DoubleVector &x, size_t index) -> void override;
    virtual auto print(unsigned int iteration, const DoubleVector &x, const DoubleVector &g, double f, double alpha, GradientBasedMethod::MethodResult result) const -> void override;

    GradientBasedMethod *gm;
    double lastFx = 0.0;

    size_t functionCount = 0;
    bool optimizeK = true;
    bool optimizeZ = true;
};

}



#endif // SOLVER6_H

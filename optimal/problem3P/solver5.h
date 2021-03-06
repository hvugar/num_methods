#ifndef SOLVER5_H
#define SOLVER5_H

#include "global.h"

#undef TIME_STEP
#define TIME_STEP 0.01
#undef TIME_MAX
#define TIME_MAX 600
#define DIMX_STEP 0.01
#define DIMX_MAX 100
#define DIMY_STEP 0.01
#define DIMY_MAX 100
#define NUM_GRAD_STEP 0.01

//#define OPTIMIZE_Q
#define OPTIMIZE_Y
#define CALCULATE_GRAD

namespace p3p5
{

class Functional;

class Shared
{
public:
    auto q(const TimeNodePDE &tn, size_t i) const -> double;
    auto s(const TimeNodePDE &tn, size_t i) const -> SpacePoint;

public:
    DoubleMatrix U = DoubleMatrix(DIMY_MAX+1, DIMX_MAX+1, 10.0);
    DoubleMatrix V = DoubleMatrix(DIMY_MAX+1, DIMX_MAX+1, 10.0);

    DoubleVector uj_v;
    DoubleVector pi_v;

    DoubleVector qi_v;
    DoubleVector gj_v;

#ifdef OPTIMIZE_Y
    DoubleMatrix k;
    DoubleMatrix z;
#endif

protected:
    Dimension _timeDimension   = Dimension(TIME_STEP, 0, TIME_MAX);
    Dimension _spaceDimensionX = Dimension(DIMX_STEP, 0, DIMX_MAX);
    Dimension _spaceDimensionY = Dimension(DIMY_STEP, 0, DIMY_MAX);
    Dimension _spaceDimensionZ = Dimension(DIMX_STEP, 0, DIMX_MAX);

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

    virtual auto layerInfo(const DoubleMatrix &u, const TimeNodePDE &tn) const -> void override;
    virtual auto timeDimension() const -> Dimension override;
    virtual auto spaceDimensionX() const -> Dimension override;
    virtual auto spaceDimensionY() const -> Dimension override;
    virtual auto spaceDimensionZ() const -> Dimension override;

    auto saveImage(const DoubleMatrix &u, const TimeNodePDE &tn) const -> void;
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

    virtual auto layerInfo(const DoubleMatrix &u, const TimeNodePDE &tn) const -> void override;
    virtual auto timeDimension() const -> Dimension override;
    virtual auto spaceDimensionX() const -> Dimension override;
    virtual auto spaceDimensionY() const -> Dimension override;
    virtual auto spaceDimensionZ() const -> Dimension override;

    auto saveImage(const DoubleMatrix &u, const TimeNodePDE &tn) const -> void;
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



#endif // SOLVER5_H

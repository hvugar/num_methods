#ifndef SOLVER2_H
#define SOLVER2_H

#include "global.h"

#define OPTIMIZE_Y
//#define OPTIMIZE_BETTA
#define OPTIMIZE_ETA

namespace p3p0
{

class CommonParameter
{
public:
    CommonParameter() { }
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
    //auto convertToVector(DoubleVector &x) const -> void;
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

    //DoubleVector *mq = nullptr;
    //DoubleVector *mp = nullptr;
    //DoubleVector *mz = nullptr;
    DoubleMatrix mq;
    DoubleMatrix mp;
    DoubleMatrix mz;

    size_t heatSourceNumber = 2;


#ifdef OPTIMIZE_Y
    DoubleMatrix alpha;
    DoubleMatrix betta;
    DoubleMatrix omega;
    DoubleVector mPnts;
    size_t measrPointNumber = 4;
    DoubleMatrix uv;
    DoubleMatrix ud;
//    DoubleVector *uv = nullptr;
//    DoubleVector *ud = nullptr;

    DoubleMatrix alphaN;
    DoubleMatrix bettaN;
    DoubleMatrix omegaN;
    DoubleVector mPntsN;
#endif

    double R       = 1.0;
    double epsilon = 0.10000;
    double no_norm = 1.00000;

    DoubleMatrix qMin;
    DoubleMatrix qMax;

//    DoubleVector *qMin = nullptr;
//    DoubleVector *qMax = nullptr;

    unsigned int _w = 7;
    unsigned int _p = 4;

protected:
    Dimension _timeDimension;
    Dimension _spaceDimensionX;

public:

//    const double VCTR_1[28] = {
//        -0.4384,-0.2379, 1.4993, -0.1458, 1.2402, 0.6780, -1.4138, -0.3162,
//        0.0000, 0.0000, 0.0000,  0.0000, 0.0000, 0.0000,  0.0000,  0.0000,
//        +4.7275, 4.0221, 5.6837, +3.3758,+3.1954, 4.5372, +4.7883, 	3.9435,
//        0.2444, 0.4913, 0.6988,  0.8745 };

//    const double NORM_1[28] = {
//        0.5112, -2.5011, 1.1734, -0.6380, 0.8051, -1.5831, 0.2565, -1.4589,
//        0.0000,  0.0000, 0.0000,  0.0000, 0.0000,  0.0000, 0.0000,  0.0000,
//        5.0249,  5.2474, 5.4654,  5.5379, 4.9167,  4.9639, 4.9973,  5.0597,
//        0.1078,  0.2612, 0.7232 , 0.9125}; /*0.00085968*/

//    const double RESULT_1[28] = {
//        0.5112, -2.5011, 1.1734, -0.6380, 0.8051, -1.5831, 0.2565, -1.4589,
//        0.0000,  0.0000, 0.0000,  0.0000, 0.0000,  0.0000, 0.0000,  0.0000,
//        5.0249,  5.2474, 5.4654,  5.5379, 4.9167,  4.9639, 4.9973,  5.0597,
//        0.1078,  0.2612, 0.7232 , 0.9125 }; /*0.00085968*/

    /*fx:29.08423687 int:  6.93129 nrm: 13.62076 penalty:  8.53218*/
    const double VCTR_1[28] = { // 1.0 1.0 1.0 0.5 0.0001 R *= 2.0 epsilon *= 0.5
        0.8308, 0.7776, 0.4642,-3.8820,-0.7526, 1.6451,-0.0113,- 1.2315,
        0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000,  0.0000,
        5.0665, 4.9530, 5.0094, 5.0119, 5.1485, 4.9594, 5.1899,  5.0338,
        0.2439, 0.3608, 0.5999, 0.7879 };  /*0.00110691*/

    const double NORM_1[28] = {
        0.0527,1.0949,-1.4853,-1.6131,-1.1310,1.3967,-1.5548,-0.2611,
        0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000,
        5.3490, 4.9490, 5.2758, 4.8753, 4.7894, 5.0547, 5.2184, 5.3647,
        0.1634,0.3895,0.6836,0.8595 };  /*0.00110691*/

    const double RESULT_1[28] = { /*0.0008164817*/
        0.0538, 1.1115, -1.4909, -1.6223, -1.1340, 1.3753, -1.5611, -0.2472,
        0.0000, 0.0000,  0.0000,  0.0000,  0.0000, 0.0000,  0.0000,  0.0000,
        5.3496, 4.9623,  5.2580,  4.8559,  4.7911, 5.0526,  5.2208,  5.3651,
        0.1618, 0.3861,  0.7037,  0.8543 };
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

};

#endif // SOLVER2_H

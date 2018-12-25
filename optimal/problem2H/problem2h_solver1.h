#ifndef PROBLEM2H_SOLVER1_H
#define PROBLEM2H_SOLVER1_H

#include "problem2h_common.h"
#include <grid/hibvp.h>

class PROBLEM2HSHARED_EXPORT Problem2HNDirichlet1 : public RnFunction, public IGradient, public InitialBoundaryValueProblemPDE, public IProjection, public IPrinter, public IVectorNormalizer
{
public:
    static void Main(int argc, char* argv[]);
    static void checkGradient1(const Problem2HNDirichlet1 &prob);
    static void checkGradient2(const Problem2HNDirichlet1 &prob);

    Problem2HNDirichlet1();
    virtual ~Problem2HNDirichlet1();

    virtual auto fx(const DoubleVector &x) const -> double;
    virtual auto gradient(const DoubleVector &, DoubleVector &) const -> void;

    virtual auto norm(const DoubleVector &v) const -> double;
    virtual auto normalize(DoubleVector &v) const -> void;

    /** Integral part of functional */
    double integral(const std::vector<DoubleMatrix> &u) const;
    double integralU(const DoubleMatrix &u) const;
    double norm(const EquationParameterH &eprm, const OptimizeParameterH &oprm, const OptimizeParameterH &r_prm) const;

    /** Penalty part of functional */
    double penalty(const spif_vectorH &info, const OptimizeParameterH &o_prm) const;
    double gpi(unsigned int i, unsigned int layer, const spif_vectorH &u_info, const OptimizeParameterH &o_prm) const;
    double g0i(unsigned int i, unsigned int layer, const spif_vectorH &u_info, const OptimizeParameterH &o_prm) const;
    double sign(double x) const;

public:
    void PrmToVector(const OptimizeParameterH &prm, DoubleVector &x) const;
    void VectorToPrm(const DoubleVector &x, OptimizeParameterH &prm) const;

    virtual double boundary(const SpaceNodePDE &, const TimeNodePDE &) const { return NAN; }

    virtual void print(unsigned int iteration, const DoubleVector &x, const DoubleVector &g, double f, double alpha, GradientMethod::MethodResult result) const;

    virtual void project(DoubleVector &x, unsigned int index);
    virtual void project(DoubleVector &x) const;

    virtual void projectControlPoints(DoubleVector &x, unsigned int index) const;
    virtual void projectMeasurePoints(DoubleVector &x, unsigned int index) const;

    void solveForwardIBVP(std::vector<DoubleMatrix> &u, spif_vectorH &u_info, bool use) const
    {
        solveForwardIBVP2(u, u_info, use);
    }

    void solveBackwardIBVP(const std::vector<DoubleMatrix> &u, spif_vectorH &p_info, bool use, const spif_vectorH &u_info) const
    {
        solveBackwardIBVP2(u, p_info, use, u_info);
    }
    //forward -------------------------------------
    void solveForwardIBVP2(std::vector<DoubleMatrix> &u, spif_vectorH &u_info, bool use) const;
    void solveBackwardIBVP2(const std::vector<DoubleMatrix> &u, spif_vectorH &p_info, bool use, const spif_vectorH &u_info) const;
    //backward ------------------------------------
    //void solveForwardIBVP1(std::vector<DoubleMatrix> &u, spif_vectorH &u_info, bool use) const;
    //void solveBackwardIBVP1(const std::vector<DoubleMatrix> &u, spif_vectorH &p_info, bool use, const spif_vectorH &u_info) const;
    //void f_initialLayers(DoubleMatrix &u00, DoubleMatrix &u10, spif_vectorH &u_info, bool use, unsigned int N, unsigned int M,
    //                     double hx, double hy, double ht, double aa__hxhx, double aa__hyhy, double lambda,
    //                     const std::vector<ExtendedSpacePointH> &qExtSpacePoints,
    //                     const std::vector<ExtendedSpacePointH> &msnExtSpacePoints,
    //                     const std::vector<ExtendedSpacePointH> &cntExtSpacePoints) const;
    //void b_initialLayers(DoubleMatrix &p00, DoubleMatrix &p10, spif_vectorH &p_info, bool use, const spif_vectorH &u_info, unsigned int N, unsigned int M,
    //                     double hx, double hy, double ht, double aa__hxhx, double aa__hyhy, double lambda,
    //                     const std::vector<ExtendedSpacePointH> &cntExtSpacePoints,
    //                     const std::vector<ExtendedSpacePointH> &msnExtSpacePoints) const;
    double f_initial1(const SpaceNodePDE &sn) const;
    double b_initial1(const SpaceNodePDE &sn) const;
    double f_initial2(const SpaceNodePDE &sn) const;
    double b_initial2(const SpaceNodePDE &sn) const;
    double f_boundary(const SpaceNodePDE &sn, const TimeNodePDE &tn) const;
    double b_boundary(const SpaceNodePDE &sn, const TimeNodePDE &tn) const;

    //void f_findRowsCols(GridH &grid, uint_vectorH &rows0, uint_vectorH &rows1, uint_vectorH &rows2, uint_vectorH &cols0, uint_vectorH &cols1, uint_vectorH &cols2, unsigned int N, unsigned int M,
    //                    const std::vector<ExtendedSpacePointH> &cntExtSpacePoints, const std::vector<ExtendedSpacePointH> &msnExtSpacePoints) const;
    //void b_findRowsCols(GridH &grid, uint_vectorH &rows0, uint_vectorH &rows1, uint_vectorH &rows2, uint_vectorH &cols0, uint_vectorH &cols1, uint_vectorH &cols2, unsigned int N, unsigned int M,
    //                    const std::vector<ExtendedSpacePointH> &msnExtSpacePoints, const std::vector<ExtendedSpacePointH> &cntExtSpacePoints) const;
    auto f_prepareInfo(unsigned int No, const std::vector<SpacePoint> &points, spif_vectorH &u_info, unsigned int LLD) const -> void;
    auto b_prepareInfo(unsigned int Nc, const std::vector<SpacePoint> &points, spif_vectorH &p_info, unsigned int LLD) const -> void;
    void f_borderLayer(DoubleMatrix &u, DoubleMatrix &uh, unsigned int ln) const;
    auto f_add2Info(const DoubleMatrix &u, spif_vectorH &u_info, unsigned int ln, double hx, double hy, const std::vector<ExtendedSpacePointH> &extMsmnts, int method = 4) const -> void;
    auto b_add2Info(const DoubleMatrix &p, spif_vectorH &p_info, unsigned int ln, double hx, double hy, const std::vector<ExtendedSpacePointH> &extCntrls, int method = 4) const -> void;
    void f_layerInfo(const DoubleMatrix &u, unsigned int ln) const;
    void b_layerInfo(const DoubleMatrix &p, unsigned int ln) const;
    auto b_characteristic(const DoubleMatrix &u, unsigned int n, unsigned int m) const -> double;

    // common -----------------------------------
    auto distributeTimeDelta(double t, double ht, unsigned int ln, const SpaceNodePDE &sn, const std::vector<ExtendedSpacePointH> &qPointNodes) const -> double;
    auto newDistributeDeltaGaussPulse(const std::vector<SpacePoint> &thetas, std::vector<ExtendedSpacePointH> &extThetas, const Dimension &dimX, const Dimension &dimY) const -> void;
    auto newDistributeDeltaGaussCntrl(const std::vector<SpacePoint> &cntrls, std::vector<ExtendedSpacePointH> &extCntrls, const Dimension &dimX, const Dimension &dimY) const -> void;
    auto newDistributeDeltaGaussMsmnt(const std::vector<SpacePoint> &msmnts, std::vector<ExtendedSpacePointH> &extMsmnts, const Dimension &dimX, const Dimension &dimY) const -> void;

    auto v(unsigned int i, OptimizeParameterH o_prm, EquationParameterH e_prm, const spif_vectorH &u_info, unsigned int ln) const -> double;
    auto mu(unsigned int, unsigned int) const -> double { return 1.0; }
public:
    EquationParameterH mEquParameter;
    OptimizeParameterH mOptParameter;
    OptimizeParameterH mRegParameter;
    double r;

    DoubleVector vmin;
    DoubleVector vmax;
    double regEpsilon;
    double noise = 0.0;

    unsigned int LD;

    bool optimizeK;
    bool optimizeZ;
    bool optimizeC;
    bool optimizeO;

    bool usePenalty;
    bool useNormal;

    GradientMethod *gm;
    //    EquationParameterHE mParameter;
    std::vector<DoubleMatrix> vu;
};

#endif // PROBLEM2H_SOLVER1_H

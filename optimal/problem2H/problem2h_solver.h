#ifndef PROBLEM2HN_H
#define PROBLEM2HN_H

#include "problem2h_common.h"

class PROBLEM2HSHARED_EXPORT Problem2HNDirichlet : public RnFunction,
        public IGradient, public InitialBoundaryValueProblemPDE, public IProjection, public IPrinter,
        public IVectorNormalizer
{
public:
    static void Main(int argc, char* argv[]);
    static void checkGradient1(const Problem2HNDirichlet &prob);
    static void checkGradient2(const Problem2HNDirichlet &prob);

    Problem2HNDirichlet();
    virtual ~Problem2HNDirichlet();

    virtual auto fx(const DoubleVector &x) const -> double;
    virtual auto gradient(const DoubleVector &, DoubleVector &) const -> void;

    virtual auto fx_norm(const DoubleVector &v) const -> double;
    virtual auto normalize(DoubleVector &v) const -> void;

    /** Integral part of functional */
    double integral(const std::vector<DoubleMatrix> &u) const;
    double integralU(const DoubleMatrix &u) const;
    double fx_norm(const EquationParameter &eprm, const OptimizeParameter &oprm, const OptimizeParameter &r_prm) const;

    /** Penalty part of functional */
    double penalty(const spif_vector &info, const OptimizeParameter &o_prm) const;
    double gpi(unsigned int i, unsigned int layer, const spif_vector &u_info, const OptimizeParameter &o_prm) const;
    double g0i(unsigned int i, unsigned int layer, const spif_vector &u_info, const OptimizeParameter &o_prm) const;
    double sign(double x) const;

public:
    void PrmToVector(const OptimizeParameter &prm, DoubleVector &x) const;
    void VectorToPrm(const DoubleVector &x, OptimizeParameter &prm) const;

    virtual double boundary(const SpaceNodePDE &, const TimeNodePDE &, BoundaryType) const { return NAN; }

    virtual void print(unsigned int iteration, const DoubleVector &x, const DoubleVector &g, double f, double alpha, GradientMethod::MethodResult result) const;

    virtual void project(DoubleVector &x, unsigned int index);
    virtual void project(DoubleVector &x) const;

    virtual void projectControlPoints(DoubleVector &x, unsigned int index) const;
    virtual void projectMeasurePoints(DoubleVector &x, unsigned int index) const;

    //forward -------------------------------------
    //backward ------------------------------------
    void solveForwardIBVP(std::vector<DoubleMatrix> &u, spif_vector &u_info, bool use) const;
    void solveBackwardIBVP(const std::vector<DoubleMatrix> &u, spif_vector &p_info, bool use, const spif_vector &u_info) const;
    void f_initialLayers(DoubleMatrix &u00, DoubleMatrix &u10, spif_vector &u_info, bool use, unsigned int N, unsigned int M,
                         double hx, double hy, double ht, double aa__hxhx, double aa__hyhy, double lambda,
                         const std::vector<ExtendedSpacePoint> &qExtSpacePoints,
                         const std::vector<ExtendedSpacePoint> &msnExtSpacePoints,
                         const std::vector<ExtendedSpacePoint> &cntExtSpacePoints) const;
    void b_initialLayers(DoubleMatrix &p00, DoubleMatrix &p10, spif_vector &p_info, bool use, const spif_vector &u_info, unsigned int N, unsigned int M,
                         double hx, double hy, double ht, double aa__hxhx, double aa__hyhy, double lambda,
                         const std::vector<ExtendedSpacePoint> &cntExtSpacePoints,
                         const std::vector<ExtendedSpacePoint> &msnExtSpacePoints) const;
    double f_initial1(const SpaceNodePDE &sn) const;
    double b_initial1(const SpaceNodePDE &sn) const;
    double f_initial2(const SpaceNodePDE &sn) const;
    double b_initial2(const SpaceNodePDE &sn) const;
    double f_boundary(const SpaceNodePDE &sn, const TimeNodePDE &tn, BoundaryType boundary = Unused) const;
    double b_boundary(const SpaceNodePDE &sn, const TimeNodePDE &tn, BoundaryType boundary = Unused) const;

    void f_findRowsCols(uint_vector &rows0, uint_vector &rows1, uint_vector &rows2, uint_vector &cols0, uint_vector &cols1, uint_vector &cols2, unsigned int N, unsigned int M,
                        const std::vector<ExtendedSpacePoint> &cntExtSpacePoints, const std::vector<ExtendedSpacePoint> &msnExtSpacePoints) const;
    void b_findRowsCols(uint_vector &rows0, uint_vector &rows1, uint_vector &rows2, uint_vector &cols0, uint_vector &cols1, uint_vector &cols2, unsigned int N, unsigned int M,
                        const std::vector<ExtendedSpacePoint> &msnExtSpacePoints, const std::vector<ExtendedSpacePoint> &cntExtSpacePoints) const;
    void f_prepareInfo(unsigned int No, const std::vector<SpacePoint> &points, spif_vector &u_info, unsigned int L, const Dimension &dimX, const Dimension &dimY) const;
    void b_prepareInfo(unsigned int Nc, const std::vector<SpacePoint> &points, spif_vector &p_info, unsigned int L, const Dimension &dimX, const Dimension &dimY) const;
    void f_borderLayer(DoubleMatrix &u, DoubleMatrix &uh, unsigned int ln) const;
    void f_add2Info(const DoubleMatrix &u, spif_vector &u_info, unsigned int ln, double hx, double hy, const std::vector<ExtendedSpacePoint> &extMsmnts, int method = 4) const;
    void b_add2Info(const DoubleMatrix &p, spif_vector &p_info, unsigned int ln, double hx, double hy, const std::vector<ExtendedSpacePoint> &extCntrls, int method = 4) const;
    void f_layerInfo(const DoubleMatrix &u, unsigned int ln) const;
    void b_layerInfo(const DoubleMatrix &p, unsigned int ln) const;
    auto b_characteristic(const DoubleMatrix &u, unsigned int n, unsigned int m) const -> double;

    // common -----------------------------------
    auto distributeTimeDelta(double t, double ht, unsigned int ln, const SpaceNodePDE &sn, const std::vector<ExtendedSpacePoint> &qPointNodes) const -> double;
    void newDistributeDeltaGaussPulse(const std::vector<SpacePoint> &thetas, std::vector<ExtendedSpacePoint> &extThetas, const Dimension &dimX, const Dimension &dimY) const;
    void newDistributeDeltaGaussCntrl(const std::vector<SpacePoint> &cntrls, std::vector<ExtendedSpacePoint> &extCntrls, const Dimension &dimX, const Dimension &dimY) const;
    void newDistributeDeltaGaussMsmnt(const std::vector<SpacePoint> &msmnts, std::vector<ExtendedSpacePoint> &extMsmnts, const Dimension &dimX, const Dimension &dimY) const;
public:
    EquationParameter mEquParameter;
    OptimizeParameter mOptParameter;
    OptimizeParameter mRegParameter;
    double r;

    DoubleVector vmin;
    DoubleVector vmax;
    double regEpsilon;

    unsigned int LD;

    bool optimizeK;
    bool optimizeZ;
    bool optimizeC;
    bool optimizeO;

    bool usePenalty;
    bool useNormal;

    GradientMethod *gm;
};

#endif // PROBLEM2H_H

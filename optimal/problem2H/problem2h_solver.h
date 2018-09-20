#ifndef PROBLEM2HN_H
#define PROBLEM2HN_H

#include "common.h"

class PROBLEM2HSHARED_EXPORT Problem2HNDirichlet : public RnFunction, public IGradient, public InitialBoundaryValueProblemPDE, public IProjection, public IPrinter
{
public:
    static void Main(int argc, char* argv[]);
    static void checkGradient(const Problem2HNDirichlet &prob);

    Problem2HNDirichlet();
    virtual ~Problem2HNDirichlet();

    virtual double fx(const DoubleVector &x) const;
    virtual void gradient(const DoubleVector &, DoubleVector &) const;

protected:
    inline double mu(double x, double y) const;
    double integral(const std::vector<DoubleMatrix> &u) const;
    double integralU(const DoubleMatrix &u) const;
    double norm(const EquationParameter &eprm, const OptimizeParameter &oprm, const OptimizeParameter &oprm0) const;

    double penalty(const spif_vector &info, const OptimizeParameter &o_prm) const;
    double gpi(unsigned int i, unsigned int layer, const spif_vector &info, const OptimizeParameter &o_prm) const;
    double g0i(unsigned int i, unsigned int layer, const spif_vector &info, const OptimizeParameter &o_prm) const;

public:
    void PrmToVector(const OptimizeParameter &prm, DoubleVector &x) const;
    void VectorToPrm(const DoubleVector &x, OptimizeParameter &prm) const;

    virtual double boundary(const SpaceNodePDE &, const TimeNodePDE &, BoundaryType) const { return NAN; }

    virtual void print(unsigned int iteration, const DoubleVector &x, const DoubleVector &g, double f, double alpha, GradientMethod::MethodResult result) const;

    virtual void project(DoubleVector &x, unsigned int index);
    virtual void project(DoubleVector &x) const;

    virtual void projectControlPoints(DoubleVector &x, unsigned int index) const;
    virtual void projectMeasurePoints(DoubleVector &x, unsigned int index) const;

    void solveForwardIBVP(std::vector<DoubleMatrix> &u, spif_vector &u_info, bool use) const;
    void solveBackwardIBVP(const std::vector<DoubleMatrix> &u, spif_vector &p_info, bool use, const spif_vector &u_info) const;

    void solveForwardIBVP1(std::vector<DoubleMatrix> &u_, spif_vector &u_info, bool use) const;
    void f_initialLayers1(DoubleMatrix &u00, DoubleMatrix &u10, spif_vector &info, bool use, espn_vector &obsPointNodes, espn_vector &cntDeltaNodes, espn_vector &qPointNodes, unsigned int N, unsigned int M, double hx, double hy, double ht, double aa__hxhx, double aa__hyhy, double lambda) const;
    void solveBackwardIBVP1(const std::vector<DoubleMatrix> &u, spif_vector &p_info, bool use, const spif_vector &u_info) const;
    void b_initialLayers1(DoubleMatrix &p00, DoubleMatrix &p10, spif_vector &p_info, bool use, espn_vector &cntPointNodes, espn_vector &obsDeltaNodes, unsigned int N, unsigned int M, double hx, double hy, double ht, double aa__hxhx, double aa__hyhy, double lambda) const;

    //forward -------------------------------------
    //void solveForwardIBVP0(std::vector<DoubleMatrix> &u_, spif_vector &u_info, bool use) const;

    double f_initial1(const SpaceNodePDE &sn) const;
    double f_initial2(const SpaceNodePDE &sn) const;
    double f_boundary(const SpaceNodePDE &sn, const TimeNodePDE &tn, BoundaryType boundary = Unused) const;

    void f_findRowsCols(uint_vector &rows0, uint_vector &rows1, uint_vector &rows2, uint_vector &cols0, uint_vector &cols1, uint_vector &cols2, espn_vector &obsPointNodes, espn_vector &cntDeltaNodes, unsigned int N, unsigned int M) const;
    void f_prepareInfo(unsigned int No, const std::vector<SpacePoint> &points, spif_vector &u_info, unsigned int L, const Dimension &dimX, const Dimension &dimY) const;
    void f_initialLayers(DoubleMatrix &u00, DoubleMatrix &u05, DoubleMatrix &u10, spif_vector &info, bool use, espn_vector &obsPointNodes, espn_vector &cntDeltaNodes, espn_vector &qPointNodes, unsigned int N, unsigned int M, double hx, double hy, double ht, double aa__hxhx, double aa__hyhy, double lambda) const;
    void f_borderLayer(DoubleMatrix &u, DoubleMatrix &uh, unsigned int ln) const;
    void f_add2Info(const DoubleMatrix &u, spif_vector &u_info, const espn_vector &obsPointNodes, unsigned int ln, double hx, double hy, int method = 4) const;
    void f_layerInfo(const DoubleMatrix &u, unsigned int ln) const;

    // backward -----------------------------------
    //void solveBackwardIBVP0(const std::vector<DoubleMatrix> &u, spif_vector &p_info, bool use, const spif_vector &u_info) const;
    double b_initial1(const SpaceNodePDE &sn) const;
    double b_initial2(const SpaceNodePDE &sn) const;
    double b_boundary(const SpaceNodePDE &sn, const TimeNodePDE &tn, BoundaryType boundary = Unused) const;

    void b_findRowsCols(uint_vector &rows0, uint_vector &rows1, uint_vector &rows2, uint_vector &cols0, uint_vector &cols1, uint_vector &cols2, espn_vector &cntPointNodes, espn_vector &obsDeltaNodes, unsigned int N, unsigned int M) const;
    void b_prepareInfo(unsigned int Nc, const std::vector<SpacePoint> &points, spif_vector &p_info, unsigned int L, const Dimension &dimX, const Dimension &dimY) const;
    void b_initialLayers(DoubleMatrix &p00, DoubleMatrix &p05, DoubleMatrix &p10, spif_vector &p_info, bool use, espn_vector &cntPointNodes, espn_vector &obsDeltaNodes, unsigned int N, unsigned int M, double hx, double hy, double ht, double aa__hxhx, double aa__hyhy, double lambda) const;
    void b_add2Info(const DoubleMatrix &p, spif_vector &p_info, const espn_vector &cntPointNodes, unsigned int ln, double hx, double hy, int method = 4) const;
    void b_layerInfo(const DoubleMatrix &p, unsigned int ln) const;

    // common -----------------------------------
    void distributeDelta0(const SpacePoint &pt, unsigned int id, espn_vector &nodes, const Dimension &dimX, const Dimension &dimY, unsigned int k, int method=4) const;
    void distributeDeltaP(const SpacePoint &pt, unsigned int id, espn_vector &nodes, const Dimension &dimX, const Dimension &dimY, unsigned int k=0) const;
    void distributeDeltaR(const SpacePoint &pt, unsigned int id, espn_vector &nodes, const Dimension &dimX, const Dimension &dimY, unsigned int k=0) const;
    void distributeDeltaG(const SpacePoint &pt, unsigned int id, espn_vector &nodes, const Dimension &dimX, const Dimension &dimY, unsigned int k) const;

public:
    EquationParameter mEquParameter;
    OptimizeParameter mOptParameter;
    OptimizeParameter mRegParameter;
    double r;

    DoubleVector vmin;
    DoubleVector vmax;
    double regEpsilon;

    DoubleMatrix V0;
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
